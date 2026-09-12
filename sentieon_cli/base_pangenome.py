"""
A base class for pangenome pipelines
"""

import argparse
import copy
import pathlib
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Set

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
from .pangenome_meta import (
    PangenomeMetadataError,
    PangenomeReference,
    detect_pangenome_reference,
    read_hapl_header,
)
from .pipeline import BasePipeline
from .shard import GRCH38_CONTIGS, detect_reference_build
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

# The pangenome reference of the graphs the pipelines shipped with, used
# when a graph's own metadata cannot be read
DEFAULT_PANGENOME_REF_NAME = "GRCh38"
DEFAULT_PANGENOME_CONTIG_PREFIX = "GRCh38#0#"

# The linear reference build that each pangenome backbone reference
# implies. Contig names alone cannot tell the builds apart, as both name
# their chromosomes chr1..chrM.
PANGENOME_REFERENCE_BUILDS: Dict[str, str] = {
    "GRCh38": "hg38",
    "CHM13": "chm13",
}


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
            # Hidden arguments. The pangenome reference name and contig
            # prefix are read from the `--gbz` file; these override the
            # detected values for a graph the detection cannot handle.
            "pangenome_ref_name": {
                "help": argparse.SUPPRESS,
            },
            "pangenome_contig_prefix": {
                "help": argparse.SUPPRESS,
            },
            "skip_contig_checks": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_pangenome_name_checks": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
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
        # The pangenome reference, resolved from the graph's metadata by
        # `resolve_pangenome_reference`. The two arguments are hidden
        # overrides, so they stay `None` unless the user supplies them.
        self.pangenome_ref_name: Optional[str] = None
        self.pangenome_contig_prefix: Optional[str] = None
        self.pangenome_reference: Optional[PangenomeReference] = None
        self.skip_contig_checks = False
        self.skip_pangenome_name_checks = False
        # Parsed by `validate`; declared here for the shared checks
        self.fai_data: Dict[str, Dict[str, int]] = {}
        self.pop_vcf_contigs: Dict[str, Optional[int]] = {}
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

    def resolve_pangenome_reference(self) -> None:
        """Resolve the pangenome reference name and contig prefix.

        Both values are read from the `--gbz` file's own metadata. The
        hidden `--pangenome_ref_name` and `--pangenome_contig_prefix`
        arguments override the detected values, and disagreeing with the
        graph is an error unless `--skip_pangenome_name_checks` is set.

        The check matters because a wrong reference name is silent: `vg
        haplotypes --set-reference` only warns, writes a sampled graph
        with no reference paths, `vg convert -Q` then writes a GFA with
        no `SN:Z:` rGFA tags, and `pgutil lift` exits 0 with every read
        unmapped.

        Called from `validate` once the reference index is parsed and
        before `validate_bundle`, which selects the model bundle's
        `extract.<name>.model` member by the resolved name.
        """
        gbz = self.required(self.gbz, "gbz")
        supplied_name = self.pangenome_ref_name
        supplied_prefix = self.pangenome_contig_prefix

        detected: Optional[PangenomeReference] = None
        try:
            detected = detect_pangenome_reference(gbz)
        except (PangenomeMetadataError, OSError) as err:
            self._undetected_pangenome_reference(gbz, err)

        if detected is None:
            self.pangenome_ref_name = (
                supplied_name or DEFAULT_PANGENOME_REF_NAME
            )
            self.pangenome_contig_prefix = (
                supplied_prefix or DEFAULT_PANGENOME_CONTIG_PREFIX
            )
        else:
            self._check_pangenome_overrides(detected)
            self.pangenome_reference = detected
            self.pangenome_ref_name = supplied_name or detected.ref_name
            self.pangenome_contig_prefix = (
                supplied_prefix or detected.contig_prefix
            )

        self.logger.info(
            "Using pangenome reference '%s' (contig prefix '%s')",
            self.pangenome_ref_name,
            self.pangenome_contig_prefix,
        )

        if detected is not None:
            self.check_pangenome_reference_fasta(detected)
        self.check_pangenome_hapl(detected)

    def _undetected_pangenome_reference(
        self, gbz: pathlib.Path, err: Exception
    ) -> None:
        """Handle a `--gbz` file whose reference could not be detected.

        Dry runs (and the unit tests behind them) parse placeholder files,
        so they fall back to the historical defaults. A real run needs
        both overrides to continue without the detected values.
        """
        if self.dry_run:
            self.logger.debug(
                "Could not read the pangenome metadata of '%s': %s", gbz, err
            )
            return
        if self.pangenome_ref_name and self.pangenome_contig_prefix:
            self.logger.warning(
                "Could not read the pangenome metadata of '%s': %s. Using "
                "the supplied `--pangenome_ref_name` and "
                "`--pangenome_contig_prefix`.",
                gbz,
                err,
            )
            return
        self.logger.error(
            "Could not read the pangenome metadata of '%s': %s. Supply "
            "both `--pangenome_ref_name` and `--pangenome_contig_prefix` "
            "to run without the detected values.",
            gbz,
            err,
        )
        sys.exit(2)

    def _check_pangenome_overrides(self, detected: PangenomeReference) -> None:
        """Compare the supplied overrides with the detected values"""
        overrides = (
            (
                "--pangenome_ref_name",
                self.pangenome_ref_name,
                detected.ref_name,
            ),
            (
                "--pangenome_contig_prefix",
                self.pangenome_contig_prefix,
                detected.contig_prefix,
            ),
        )
        mismatched = False
        for flag, supplied, found in overrides:
            if supplied is None or supplied == found:
                continue
            mismatched = True
            if self.skip_pangenome_name_checks:
                self.logger.warning(
                    "The supplied `%s` value '%s' does not match the value "
                    "detected in the pangenome graph, '%s'. Using the "
                    "supplied value.",
                    flag,
                    supplied,
                    found,
                )
            else:
                self.logger.error(
                    "The supplied `%s` value '%s' does not match the value "
                    "detected in the pangenome graph, '%s'. Drop the "
                    "argument to use the detected value.",
                    flag,
                    supplied,
                    found,
                )
        if mismatched and not self.skip_pangenome_name_checks:
            sys.exit(2)

    def _pangenome_check_failed(self, msg: str, *args: object) -> None:
        """Fail a pangenome consistency check.

        `--skip_pangenome_name_checks` downgrades every one of these
        checks to a warning.
        """
        if self.skip_pangenome_name_checks:
            self.logger.warning(msg, *args)
            return
        self.logger.error(msg, *args)
        sys.exit(2)

    def check_pangenome_reference_fasta(
        self, detected: PangenomeReference
    ) -> None:
        """Check the reference FASTA against the pangenome backbone.

        Every contig of the backbone reference must be in the reference
        index, and the build the backbone implies must match the build
        detected from the index. Both checks are needed: contig names
        alone cannot tell GRCh38 from CHM13.
        """
        missing = [ctg for ctg in detected.contigs if ctg not in self.fai_data]
        if missing:
            self._pangenome_check_failed(
                "%d contig(s) of the pangenome reference '%s' are missing "
                "from the reference FASTA index: %s",
                len(missing),
                detected.ref_name,
                ", ".join(missing[:10]),
            )

        expected_build = PANGENOME_REFERENCE_BUILDS.get(detected.ref_name)
        fai_build = detect_reference_build(self.fai_data)
        if expected_build is None or fai_build is None:
            self.logger.info(
                "The pangenome reference is '%s' and the reference FASTA "
                "build is '%s'; the two were not compared",
                detected.ref_name,
                fai_build,
            )
            return
        if expected_build != fai_build:
            self._pangenome_check_failed(
                "The pangenome reference '%s' expects the '%s' reference "
                "build, but the `--reference` FASTA is '%s'",
                detected.ref_name,
                expected_build,
                fai_build,
            )

    def check_pangenome_hapl(
        self, detected: Optional[PangenomeReference]
    ) -> None:
        """Check the `--hapl` file against the `--gbz` graph.

        A mismatched top-level chain count is only a warning: the count
        equals the backbone's path count for every graph seen so far, but
        it is a property of the snarl decomposition rather than a
        guarantee.
        """
        hapl = self.required(self.hapl, "hapl")
        try:
            header = read_hapl_header(hapl)
        except (PangenomeMetadataError, OSError) as err:
            if self.dry_run:
                self.logger.debug(
                    "Could not read the haplotype file '%s': %s", hapl, err
                )
                return
            self._pangenome_check_failed("%s", err)
            return

        if detected is None:
            return
        if header.top_level_chains != detected.n_paths:
            self.logger.warning(
                "The `--hapl` file '%s' has %d top-level chains, but the "
                "pangenome reference '%s' has %d paths. The `--gbz` and "
                "`--hapl` files may not be a matching pair.",
                hapl,
                header.top_level_chains,
                detected.ref_name,
                detected.n_paths,
            )

    def validate_grch38_contigs(self) -> None:
        """Check the reference and pop VCF contig lengths against GRCh38.

        The lengths only describe GRCh38; other pangenome references are
        covered by the reference FASTA checks in
        `resolve_pangenome_reference` instead.
        """
        if self.skip_contig_checks:
            return
        if self.pangenome_ref_name != DEFAULT_PANGENOME_REF_NAME:
            self.logger.info(
                "The pangenome reference is '%s', so the GRCh38 "
                "contig-length check is skipped",
                self.pangenome_ref_name,
            )
            return

        # Check the fai file contigs
        mismatch_contigs: Set[str] = set()
        for ctg, length in GRCH38_CONTIGS.items():
            fai_length = self.fai_data.get(ctg, {}).get("length", -1)
            if length != fai_length:
                mismatch_contigs.add(ctg)
        if mismatch_contigs:
            self.logger.error(
                "Reference contigs with unexpected lengths: %s",
                ", ".join(mismatch_contigs),
            )
            sys.exit(2)

        # Check the pop VCF file contigs
        if self.dry_run:
            return
        mismatch_contigs = set()
        for ctg, length in GRCH38_CONTIGS.items():
            if length != self.pop_vcf_contigs.get(ctg, -1):
                mismatch_contigs.add(ctg)
        if mismatch_contigs:
            self.logger.error(
                "Pop VCF contigs with unexpected lengths: %s",
                ", ".join(mismatch_contigs),
            )
            sys.exit(2)

    def ref_name(self) -> str:
        """The resolved pangenome reference name.

        `resolve_pangenome_reference` always sets it; this narrows the
        optional attribute for the consumers of the graph.
        """
        return self.required(self.pangenome_ref_name, "pangenome_ref_name")

    def contig_prefix(self) -> str:
        """The resolved pangenome contig prefix"""
        return self.required(
            self.pangenome_contig_prefix, "pangenome_contig_prefix"
        )

    def build_check_gbz_job(self, sample_gbz: pathlib.Path) -> Job:
        """Check the sampled pangenome before anything consumes it"""
        return Job(
            cmds.cmd_check_pangenome_gbz(
                sample_gbz, self.ref_name(), self.contig_prefix()
            ),
            "check-sample-gbz",
            1,
            task_name="pangenome",
        )

    def build_check_gfa_job(self, sample_gfa: pathlib.Path) -> Job:
        """Check the converted GFA before anything consumes it"""
        return Job(
            cmds.cmd_check_pangenome_gfa(
                sample_gfa, self.ref_name(), self.contig_prefix()
            ),
            "check-sample-gfa",
            1,
            task_name="pangenome",
        )

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
