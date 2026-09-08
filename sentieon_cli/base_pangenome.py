"""
A base class for pangenome pipelines
"""

import copy
import pathlib
import sys
from typing import Iterable, List, Optional, Sequence

from . import command_strings as cmds
from .dag import DAG
from .driver import (
    AlignmentStat,
    BaseAlgo,
    BaseDistributionByCycle,
    CoverageMetrics,
    GCBias,
    InsertSizeMetricAlgo,
    MeanQualityByCycle,
    QualDistribution,
    WgsMetricsAlgo,
)
from .job import Job
from .pipeline import BasePipeline
from .stages.base import StageContext
from .stages.cnv import PangenomeCNVResult, PangenomeCNVStage
from .stages.expansion import ExpansionHunterResult, ExpansionHunterStage
from .stages.metrics import MetricsPaths
from .stages.segdup import SegdupResult, SegdupStage
from .stages.t1k import (
    DEFAULT_T1K_HLA_LOCUS,
    DEFAULT_T1K_KIR_LOCUS,
    T1K_MIN_VERSIONS,
    T1KResult,
    T1KStage,
)
from .util import path_arg, require_versions, sample_sex_arg


class BasePangenome(BasePipeline):
    """A pipeline base class for short reads"""

    params = copy.deepcopy(BasePipeline.params)
    params.update(
        {
            # Required arguments
            "gbz": {
                "help": "The pangenome graph file in GBZ format.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "hapl": {
                "help": "The haplotype file.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "model_bundle": {
                "flags": ["-m", "--model_bundle"],
                "help": "The model bundle file.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "r1_fastq": {
                "nargs": "*",
                "help": "Sample R1 fastq files.",
                "type": path_arg(exists=True, is_file=True),
            },
            "r2_fastq": {
                "nargs": "*",
                "help": "Sample R2 fastq files.",
                "type": path_arg(exists=True, is_file=True),
            },
            # Additional arguments
            "bam_format": {
                "help": (
                    "Use the BAM format instead of CRAM for output aligned "
                    "files."
                ),
                "action": "store_true",
            },
            "dbsnp": {
                "flags": ["-d", "--dbsnp"],
                "help": (
                    "dbSNP vcf file Supplying this file will annotate "
                    "variants with their dbSNP refSNP ID numbers."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "kmer_memory": {
                "help": "Memory limit for KMC in GB.",
                "default": 30,
                "type": int,
            },
            "pcr_free": {
                "help": "Use arguments for PCR-free data processing",
                "action": "store_true",
            },
            "par_bed": {
                "help": (
                    "A BED file of the pseudo-autosomal regions (PAR), used "
                    "for sex-aware CNV calling of male samples. Overrides the "
                    "PAR BED file selected for the reference genome."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "sample_sex": {
                "help": (
                    "The sample sex, used by the sex-aware callers. "
                    "Supplying this argument overrides the sex estimated "
                    "from read coverage."
                ),
                "metavar": "{male,female}",
                "type": sample_sex_arg,
            },
            "segdup_caller": {
                "nargs": "*",
                "help": (
                    "Call variants in difficult segmental duplications with "
                    "segdup-caller. Supply the flag with no arguments to run "
                    "the caller's default gene set. Supply a comma-separated "
                    "list of gene names (e.g. 'CFH,CYP2D6,SMN1') to restrict "
                    "calling to those genes."
                ),
            },
            "expansion_catalog": {
                "help": (
                    "An ExpansionHunter variant catalog. Required for short "
                    "tandem repeat expansion calling."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "t1k_hla_seq": {
                "help": (
                    "The DNA HLA seq FASTA file for T1K. Required for HLA "
                    "calling."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "t1k_hla_coord": {
                "help": (
                    "The DNA HLA coord FASTA file for T1K. Required for HLA "
                    "calling."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "t1k_hla_locus": {
                "default": DEFAULT_T1K_HLA_LOCUS,
                "help": (
                    "Reference interval covering the HLA locus. Reads "
                    "overlapping this region are extracted before being "
                    f"passed to T1K (default: {DEFAULT_T1K_HLA_LOCUS})."
                ),
            },
            "t1k_kir_seq": {
                "help": (
                    "The DNA KIR seq FASTA file for T1K. Required for KIR "
                    "calling."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "t1k_kir_coord": {
                "help": (
                    "The DNA KIR coord FASTA file for T1K. Required for KIR "
                    "calling."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "t1k_kir_locus": {
                "default": DEFAULT_T1K_KIR_LOCUS,
                "help": (
                    "Reference interval covering the KIR locus. Reads "
                    "overlapping this region are extracted before being "
                    f"passed to T1K (default: {DEFAULT_T1K_KIR_LOCUS})."
                ),
            },
        }
    )

    positionals = BasePipeline.positionals

    def __init__(self) -> None:
        super().__init__()
        self.gbz: Optional[pathlib.Path] = None
        self.hapl: Optional[pathlib.Path] = None
        self.model_bundle: Optional[pathlib.Path] = None
        self.r1_fastq: List[pathlib.Path] = []
        self.r2_fastq: List[pathlib.Path] = []
        self.bam_format = False
        self.dbsnp: Optional[pathlib.Path] = None
        self.kmer_memory = 30
        self.pcr_free = False
        self.par_bed: Optional[pathlib.Path] = None
        self.segdup_caller: Optional[List[str]] = None
        self.expansion_catalog: Optional[pathlib.Path] = None
        self.t1k_hla_seq: Optional[pathlib.Path] = None
        self.t1k_hla_coord: Optional[pathlib.Path] = None
        self.t1k_hla_locus: str = DEFAULT_T1K_HLA_LOCUS
        self.t1k_kir_seq: Optional[pathlib.Path] = None
        self.t1k_kir_coord: Optional[pathlib.Path] = None
        self.t1k_kir_locus: str = DEFAULT_T1K_KIR_LOCUS
        self.has_cnv_model = False
        # Stashed by `build_dag` for the second, sex-aware DAG
        self.ploidy_json: Optional[pathlib.Path] = None
        self.sr_alignments: List[pathlib.Path] = []
        self.sr_replace_rg: Optional[List[List[str]]] = None

    def validate_t1k(self) -> None:
        hla_requested = self.t1k_hla_seq is not None or (
            self.t1k_hla_coord is not None
        )
        if hla_requested and not (self.t1k_hla_seq and self.t1k_hla_coord):
            self.logger.error(
                "For HLA calling, both `--t1k_hla_seq` and `--t1k_hla_coord` "
                "must be supplied."
            )
            sys.exit(2)

        kir_requested = self.t1k_kir_seq is not None or (
            self.t1k_kir_coord is not None
        )
        if kir_requested and not (self.t1k_kir_seq and self.t1k_kir_coord):
            self.logger.error(
                "For KIR calling, both `--t1k_kir_seq` and `--t1k_kir_coord` "
                "must be supplied."
            )
            sys.exit(2)

        if not (hla_requested or kir_requested):
            return

        require_versions(T1K_MIN_VERSIONS, skip=self.skip_version_check)

    def output_path(self, suffix: str) -> pathlib.Path:
        """A path next to the output VCF, with `.vcf.gz` replaced"""
        output_vcf = self.required(self.output_vcf, "output_vcf")
        return pathlib.Path(str(output_vcf).replace(".vcf.gz", suffix))

    def add_t1k(
        self,
        dag: DAG,
        ctx: StageContext,
        inputs: List[pathlib.Path],
        upstream: Iterable[Job] = (),
    ) -> List[T1KResult]:
        """Genotype every enabled T1K gene group from `inputs`"""
        deps = set(upstream)
        gene_groups = (
            (
                "hla",
                "hla-wgs",
                self.t1k_hla_seq,
                self.t1k_hla_coord,
                self.t1k_hla_locus,
            ),
            (
                "kir",
                "kir-wgs",
                self.t1k_kir_seq,
                self.t1k_kir_coord,
                self.t1k_kir_locus,
            ),
        )

        results: List[T1KResult] = []
        for tag, preset, gene_seq, gene_coord, locus in gene_groups:
            if not (gene_seq and gene_coord):
                continue
            results.append(
                T1KStage(
                    ctx=ctx,
                    inputs=inputs,
                    gene_seq=gene_seq,
                    gene_coord=gene_coord,
                    locus=locus,
                    preset=preset,
                    tag=tag,
                    output_dir=self.output_path(f"_{tag}"),
                ).add_to(dag, deps)
            )
        return results

    def add_expansion(
        self,
        dag: DAG,
        ctx: StageContext,
        alignment: pathlib.Path,
    ) -> Optional[ExpansionHunterResult]:
        """Genotype repeat expansions, when a catalog was supplied"""
        if self.expansion_catalog is None:
            return None

        return ExpansionHunterStage(
            ctx=ctx,
            alignment=alignment,
            variant_catalog=self.expansion_catalog,
            output_prefix=self.output_path("_expansion"),
            sample_sex=self.sample_sex,
        ).add_to(dag)

    def add_segdup(
        self,
        dag: DAG,
        ctx: StageContext,
        sr_alignment: pathlib.Path,
        *,
        lr_alignment: Optional[pathlib.Path] = None,
        lr_bundle: Optional[pathlib.Path] = None,
        overrides: Sequence[str] = (),
    ) -> Optional[SegdupResult]:
        """Call segmental duplications, when `--segdup_caller` was set"""
        if self.segdup_caller is None:
            return None

        return SegdupStage(
            ctx=ctx,
            sr_alignment=sr_alignment,
            sr_bundle=self.required(self.model_bundle, "model_bundle"),
            output_dir=self.output_path("_segdups"),
            input_vcf=ctx.output_vcf,
            sample_sex=self.sample_sex,
            genes=",".join(self.segdup_caller) or None,
            lr_alignment=lr_alignment,
            lr_bundle=lr_bundle,
            overrides=overrides,
        ).add_to(dag)

    def add_pangenome_cnv(
        self,
        dag: DAG,
        ctx: StageContext,
        sv_vcf: pathlib.Path,
        inputs: List[pathlib.Path],
        *,
        replace_rg: Optional[List[List[str]]] = None,
        interval: Optional[pathlib.Path] = None,
    ) -> PangenomeCNVResult:
        """Call CNVs from the short reads and the PangenomeSV output"""
        bundle = self.required(self.model_bundle, "model_bundle")

        return PangenomeCNVStage(
            ctx=ctx,
            inputs=inputs,
            model=bundle.joinpath("cnv.model"),
            sv_vcf=sv_vcf,
            sample_sex=self.sample_sex,
            par_bed=self.cnv_par_bed,
            interval=interval,
            replace_rg=replace_rg,
        ).add_to(dag)

    def build_kmc_job(
        self, kmer_prefix: pathlib.Path, job_threads: int
    ) -> Job:
        """Build KMC k-mer counting jobs"""
        # Create file list for KMC
        file_list = pathlib.Path(str(kmer_prefix) + ".paths")
        all_fastqs = []

        # Add R1 files
        all_fastqs.extend(self.r1_fastq)

        # Add R2 files if present
        if self.r2_fastq:
            all_fastqs.extend(self.r2_fastq)

        # Write file list
        if not self.dry_run:
            with open(file_list, "w") as f:
                for fq in all_fastqs:
                    f.write(f"{fq}\n")

        # Create KMC job
        kmc_job = Job(
            cmds.cmd_kmc(
                kmer_prefix,
                file_list,
                self.tmp_dir,
                memory=self.kmer_memory,
                threads=self.cores,
            ),
            "kmc",
            job_threads,
            task_name="kmer-counting",
        )

        return kmc_job

    def pangenome_metrics_algos(self, paths: MetricsPaths) -> List[BaseAlgo]:
        """The metrics collected from a deduplicated pangenome alignment"""
        return [
            InsertSizeMetricAlgo(paths.insert_size),
            MeanQualityByCycle(paths.mean_qual_by_cycle),
            BaseDistributionByCycle(paths.base_distribution_by_cycle),
            QualDistribution(paths.qual_distribution),
            AlignmentStat(paths.alignment_stat),
            GCBias(paths.gc_bias, summary=paths.gc_bias_summary),
            WgsMetricsAlgo(paths.wgs, include_unpaired="true"),
            CoverageMetrics(paths.coverage),
        ]
