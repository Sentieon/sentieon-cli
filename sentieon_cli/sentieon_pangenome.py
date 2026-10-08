"""
Sentieon's pangenome alignment and variant calling pipeline
"""

import argparse
import copy
import json
import pathlib
import sys
from typing import Dict, List, Optional, Set, Union

import packaging.version

from importlib.resources import files

from . import command_strings as cmds
from .archive import ar_load
from .base_pangenome import BasePangenome
from .dag import DAG
from .driver import (
    BaseAlgo,
    DNAscope,
    PangenomeSV,
)
from .job import Job
from .logging import get_logger
from .shell_pipeline import Command, Pipeline
from .stages.alignment import (
    BwaExtractStage,
    find_unzip,
)
from .stages.base import StageContext
from .stages.cnv import CNV_MIN_VERSIONS
from .stages.dedup import DedupStage
from .stages.expansion import EXPANSION_MIN_VERSIONS
from .stages.metrics import MetricsPaths, MetricsStage
from .stages.ploidy import PloidyStage
from .stages.segdup import SEGDUP_MIN_VERSIONS
from .stages.small_variants import (
    ApplySpec,
    DNAscopeStage,
    GVCFtyperStage,
    TransferApplyStage,
    TransferSpec,
)
from .stages.transfer import TransferConfig
from .util import (
    __version__,
    check_kmc_patch,
    parse_rg_line,
    path_arg,
    require_versions,
    total_memory,
)
from .shard import (
    determine_shards_from_fai,
    parse_fai,
)

SENT_PANGENOME_MIN_VERSIONS = {
    "kmc": None,
    "sentieon driver": packaging.version.Version("202503.02"),
    "vg": None,
    "bcftools": packaging.version.Version("1.22"),
    "samtools": packaging.version.Version("1.16"),
}

logger = get_logger(__name__)


class SentieonPangenome(BasePangenome):
    """The Sentieon pangenome pipeline"""

    params = copy.deepcopy(BasePangenome.params)
    params.update(
        {
            # Required arguments
            "readgroup": {
                "help": "Readgroup information for the fastq files.",
            },
            "sample_input": {
                "flags": ["-i", "--sample_input"],
                "nargs": "*",
                "help": "sample BAM or CRAM file.",
                "type": path_arg(exists=True, is_file=True),
            },
            "pop_vcf": {
                "flags": ["--pop_vcf"],
                "help": (
                    "A VCF containing annotations for use with DNAModelApply."
                ),
                "type": path_arg(exists=True, is_file=True),
                "required": True,
            },
            # Additional arguments
            "bed": {
                "flags": ["-b", "--bed"],
                "help": (
                    "Region BED file. Supplying this file will limit variant "
                    "calling to the intervals inside the BED file."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "call_svs": {
                "help": "Call structural and copy number variants.",
                "action": "store_true",
            },
            "gvcf": {
                "flags": ["-g", "--gvcf"],
                "help": "Generate a gVCF output file.",
                "action": "store_true",
            },
            "skip_metrics": {
                "help": "Skip metrics collection and multiQC",
                "action": "store_true",
            },
            "skip_multiqc": {
                "help": "Skip multiQC report generation",
                "action": "store_true",
            },
            # Hidden arguments
            "skip_pop_vcf_id_check": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_model_apply": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_small_variants": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
        }
    )

    positionals = BasePangenome.positionals

    def __init__(self) -> None:
        super().__init__()
        self.readgroup: Optional[str] = None
        self.sample_input: List[pathlib.Path] = []
        self.pop_vcf: Optional[pathlib.Path] = None
        self.bed: Optional[pathlib.Path] = None
        self.call_svs = False
        self.gvcf = False
        self.extract_model_name = "extract.model"
        self.skip_metrics = False
        self.skip_multiqc = False
        self.skip_pop_vcf_id_check: bool = False
        self.skip_model_apply = False
        self.skip_small_variants = False

    def _cnv_in_second_dag(self) -> bool:
        """CNV calling runs in the second, sex-aware DAG"""
        return self.call_svs and self.has_cnv_model

    def cnv_combine_preset(self) -> str:
        """Ultima bundles carry a single-end CNVscope model"""
        return "SE" if self.tech.upper() == "ULTIMA" else "PE"

    def _needs_second_dag(self) -> bool:
        """The run has jobs that depend on the estimated sample sex"""
        return bool(
            self.expansion_catalog
            or self.segdup_caller is not None
            or self._cnv_in_second_dag()
        )

    def validate(self) -> None:
        """Validate pipeline inputs"""
        self.validate_ref()
        self.fai_data = parse_fai(pathlib.Path(str(self.reference) + ".fai"))
        self.shards = determine_shards_from_fai(
            self.fai_data, 10 * 1000 * 1000
        )
        # Before `validate_bundle`, which picks the bundle's
        # `extract.<pangenome_ref_name>.model` member by the resolved name
        self.resolve_pangenome_reference()
        self.load_pop_vcf_header(self.pop_vcf)

        self.validate_bundle()
        self.validate_fastq_rg()
        self.validate_output_vcf()
        self.collect_readgroups()

        if not self.sample_input and not self.r1_fastq:
            self.logger.error(
                "Please supply either the `--sample_input` or `--r1_fastq` "
                "and `--readgroups` arguments"
            )
            sys.exit(2)

        if self.sample_input and self.r1_fastq:
            self.logger.error(
                "Supplying both `--r1_fastq` and `--sample_input` is not "
                "supported"
            )
            sys.exit(2)

        if self.r1_fastq:
            self.validate_bwa_index()

        self.validate_segdup()
        self.validate_expansion()
        self.validate_t1k()
        self.validate_cnv()

        require_versions(
            SENT_PANGENOME_MIN_VERSIONS, skip=self.skip_version_check
        )

        if not self.skip_version_check and self.sample_input:
            if not check_kmc_patch("kmc"):
                self.logger.error(
                    "Error: The 'kmc' executable in the PATH does not "
                    "support reading from stdin. Please ensure "
                    "you are using the patched version of KMC from "
                    "https://github.com/Sentieon/KMC/releases."
                )
                sys.exit(2)

        if self.bed is None:
            self.logger.info(
                "A BED file is recommended to avoid variant calling "
                "across decoy and unplaced contigs."
            )

        self.validate_pangenome_contig_lengths()

    def validate_segdup(self) -> None:
        if self.segdup_caller is None:
            return

        if len(self.sample_input) > 1:
            self.logger.error(
                "`--segdup_caller` accepts only a single `--sample_input` "
                "BAM/CRAM file."
            )
            sys.exit(2)

        if self.skip_small_variants:
            self.logger.error(
                "`--segdup_caller` requires the small-variant VCF as "
                "`--input_vcf` and cannot be combined with "
                "`--skip_small_variants`."
            )
            sys.exit(2)

        require_versions(SEGDUP_MIN_VERSIONS, skip=self.skip_version_check)

    def validate_expansion(self) -> None:
        if self.expansion_catalog is None:
            return

        if self.tech.upper() == "ULTIMA":
            self.logger.error(
                "`--expansion_catalog` is not supported with single-end "
                "(Ultima) input."
            )
            sys.exit(2)

        if len(self.sample_input) > 1:
            self.logger.error(
                "`--expansion_catalog` accepts only a single `--sample_input` "
                "BAM/CRAM file."
            )
            sys.exit(2)

        require_versions(EXPANSION_MIN_VERSIONS, skip=self.skip_version_check)

    def validate_cnv(self) -> None:
        """Validate the arguments used for sex-aware CNV calling"""
        # `validate_bundle` has already set `self.has_cnv_model`
        cnv_will_run = self._cnv_in_second_dag()
        self.resolve_cnv_par_bed(self.fai_data, self.par_bed, cnv_will_run)
        if not cnv_will_run:
            return

        # A PAR BED file is required for CNV calling, whatever the sex
        self.validate_cnv_par(True)

        require_versions(CNV_MIN_VERSIONS, skip=self.skip_version_check)

    def validate_bundle(self) -> None:
        self.required(self.pop_vcf, "pop_vcf")
        gbz = self.required(self.gbz, "gbz")
        bundle_info_bytes = ar_load(
            str(self.model_bundle) + "/bundle_info.json"
        )
        if isinstance(bundle_info_bytes, list):
            bundle_info_bytes = b"{}"

        bundle_info = json.loads(bundle_info_bytes.decode())
        try:
            req_version = packaging.version.Version(
                bundle_info["minScriptVersion"]
            )
            bundle_pipeline = bundle_info["pipeline"]
            self.tech: str = bundle_info["platform"].upper()
            bundle_vcf_id = bundle_info["SentieonVcfID"]
            bundle_pangenome = bundle_info.get(
                "pangenome", "hprc-v2.0-mc-grch38.gbz"
            )
        except KeyError:
            self.logger.error(
                "The model bundle does not have the expected attributes"
            )
            sys.exit(2)
        if req_version > packaging.version.Version(__version__):
            self.logger.error(
                "The model bundle requires version %s or later of the "
                "sentieon-cli.",
                req_version,
            )
            sys.exit(2)
        if bundle_pipeline != "Sentieon pangenome":
            self.logger.error("The model bundle is for a different pipeline.")
            sys.exit(2)
        if gbz.name != bundle_pangenome:
            self.logger.warning(
                "The `--gbz` file name is not '%s'. "
                "This model is optimized for the %s pangenome.",
                bundle_pangenome,
                bundle_pangenome,
            )

        bundle_members = set(ar_load(str(self.model_bundle)))

        # Prefer a reference-specific extract model. Fall back to the generic
        # 'extract.model' only for the default 'GRCh38' reference so that
        # existing bundles continue to work.
        extract_candidate = f"extract.{self.pangenome_ref_name}.model"
        if extract_candidate in bundle_members:
            self.extract_model_name = extract_candidate
        elif self.pangenome_ref_name == "GRCh38":
            self.extract_model_name = "extract.model"
        else:
            self.extract_model_name = extract_candidate

        if (
            "dnascope.model" not in bundle_members
            or self.extract_model_name not in bundle_members
            or "minimap2.model" not in bundle_members
        ):
            self.logger.error(
                "Expected model files not found in the model bundle file"
            )
            sys.exit(2)

        self.has_cnv_model = "cnv.model" in bundle_members
        if self.call_svs and not self.has_cnv_model:
            self.logger.warning(
                "The model bundle does not contain a 'cnv.model' file. "
                "CNV calling with CNVscope will be skipped; SV calling with "
                "PangenomeSV will still run."
            )

        if not self.skip_pop_vcf_id_check:
            self.check_pop_vcf_id(bundle_vcf_id)

    def validate_fastq_rg(self) -> None:
        if len(self.r1_fastq) != len(self.r2_fastq):
            self.logger.error(
                "The number of input `--r1_fastq` files does not equal the "
                "number of `--r2_fastq` files"
            )
            sys.exit(2)

        if (len(self.r1_fastq) > 0 and not self.readgroup) or (
            self.readgroup and len(self.r1_fastq) < 1
        ):
            self.logger.error(
                "`--r1_fastq`, `--r2_fastq`, and `--readgroup` are required "
                "with fastq input. `--readgroup` cannot be used with bam/cram "
                "input."
            )
            sys.exit(2)

    def collect_readgroups(self) -> None:
        """Collect readgroup tags"""
        self.bam_readgroups: List[Dict[str, str]] = []
        for aln in self.sample_input:
            aln_rgs = cmds.get_rg_lines(aln, self.dry_run)
            for rg_line in aln_rgs:
                self.bam_readgroups.append(parse_rg_line(rg_line))
                break  # currently just the first RG line for each input

        self.fastq_readgroup: Dict[str, str] = {}
        if not self.readgroup:
            return

        try:
            parsed_rg = parse_rg_line(self.readgroup.replace(r"\t", "\t"))
        except ValueError as e:
            self.logger.error(
                "Invalid --readgroup value '%s': %s", self.readgroup, e
            )
            sys.exit(2)
        if not parsed_rg.get("ID"):
            self.logger.error(
                "Readgroup '%s' does not have a RGID tag",
                self.readgroup,
            )
            sys.exit(2)
        if parsed_rg.get("SM", None) is None:
            self.logger.error(
                "Readgroup '%s' does not have a RGSM tag",
                self.readgroup,
            )
        self.fastq_readgroup = parsed_rg

    def configure(self) -> None:
        """Configure pipeline parameters"""
        pass

    def build_dag(self) -> DAG:
        """Build the first DAG for the Sentieon pangenome pipeline"""
        bundle = self.required(self.model_bundle, "model_bundle")

        ctx = self.stage_context()

        self.logger.info("Building the Sentieon pangenome DAG")
        dag = DAG()

        # Output files
        suffix = "bam" if self.bam_format else "cram"
        out_bwa_aln = pathlib.Path(
            str(ctx.output_vcf).replace(".vcf.gz", f"_bwa_deduped.{suffix}")
        )
        out_mm2_aln = pathlib.Path(
            str(ctx.output_vcf).replace(".vcf.gz", f"_mm2_deduped.{suffix}")
        )
        out_gvcf = pathlib.Path(
            str(ctx.output_vcf).replace(".vcf.gz", ".g.vcf.gz")
        )

        # Intermediate file paths
        bwa_bam = self.tmp_dir.joinpath("sample-bwa.bam")
        ext_fastq = self.tmp_dir.joinpath("sample-bwa.fq.gz")
        kmer_prefix = self.tmp_dir.joinpath("sample.fq")
        kmer_file = pathlib.Path(str(kmer_prefix) + ".kff")
        sample_pangenome = self.tmp_dir.joinpath("sample_pangenome.gbz")
        sample_gfa = self.tmp_dir.joinpath("sample-hap.gfa")
        sample_fasta = self.tmp_dir.joinpath("sample-hap.fa")
        mm2_bam = self.tmp_dir.joinpath("sample-mm2.bam")
        if not self.r1_fastq:
            # with bam/cram input, output the realigned bam/cram
            mm2_bam = pathlib.Path(
                str(ctx.output_vcf).replace(
                    ".vcf.gz", f"_mm2_deduped.{suffix}"
                )
            )
        raw_vcf = self.tmp_dir.joinpath("sample-dnascope.vcf.gz")
        transfer_vcf = self.tmp_dir.joinpath("sample-dnascope_transfer.vcf.gz")

        bwa_lc_dependencies: Set[Job] = set()
        haplotype_dependencies: Set[Job] = set()
        mm2_dependencies: Set[Job] = set()
        dnascope_bams: List[pathlib.Path] = []

        total_mem_gb = total_memory() / (1024.0**3)

        if self.r1_fastq:
            # KMC k-mer counting
            kmc_job = self.build_kmc_job(kmer_prefix, 0)  # run in background
            dag.add_job(kmc_job)
            haplotype_dependencies.add(kmc_job)

            # BWA alignment and extraction
            bwa_result = self.bwa_extract_stage(
                ctx, bwa_bam, ext_fastq
            ).add_to(dag)
            bwa_job = bwa_result.jobs[0]
            mm2_dependencies.add(bwa_job)
            bwa_lc_dependencies.add(bwa_job)
            # Do not run vg-haplotypes with bwa in low-mem environments
            if total_mem_gb < 70:
                haplotype_dependencies.add(bwa_job)
        else:
            dnascope_bams = copy.deepcopy(self.sample_input)
            # ReadWriter cannot write to /dev/stdout directly; pre-create a
            # symlink so the driver writes to a real path that resolves to
            # its stdout (the pipe to pgutil extract).
            rw_bam = self.tmp_dir.joinpath("extract-kmc-rw.bam")
            ln_job = Job(
                Pipeline(Command("ln", "-sf", "/dev/stdout", str(rw_bam))),
                "extract-kmc-symlink",
                1,
                task_name="read-extraction",
            )
            dag.add_job(ln_job)

            extract_kmc_job = Job(
                cmds.cmd_extract_kmc(
                    kmer_prefix,
                    ext_fastq,
                    self.sample_input,
                    ctx.reference,
                    bundle.joinpath(self.extract_model_name),
                    self.tmp_dir,
                    rw_bam,
                    threads=self.cores,
                ),
                "extract-kmc",
                self.cores,
                task_name="read-extraction",
            )
            dag.add_job(extract_kmc_job, {ln_job})
            haplotype_dependencies.add(extract_kmc_job)

        # vg haplotypes - create a sample-specific pangenome
        haplotypes_job = self.build_haplotypes_job(sample_pangenome, kmer_file)
        dag.add_job(haplotypes_job, haplotype_dependencies)

        # Confirm the sampled pangenome kept its reference paths before
        # anything consumes it
        check_gbz_job = self.build_check_gbz_job(sample_pangenome)
        dag.add_job(check_gbz_job, {haplotypes_job})

        # convert the sample pangenome
        gfa_job = self.build_gfa_job(sample_gfa, sample_pangenome)
        fasta_job = self.build_fasta_job(sample_fasta, sample_pangenome)
        dag.add_job(gfa_job, {check_gbz_job})
        dag.add_job(fasta_job, {check_gbz_job})

        # Confirm the GFA carries the reference's rGFA tags, without which
        # `pgutil lift` leaves every read unmapped
        check_gfa_job = self.build_check_gfa_job(sample_gfa)
        dag.add_job(check_gfa_job, {gfa_job})

        # minimap2 alignment of the extracted fastq
        dnascope_dependencies = set()
        mm2_job = self.build_minimap2_lift_job(
            mm2_bam,
            ext_fastq,
            sample_fasta,
            sample_gfa,
        )
        dag.add_job(mm2_job, mm2_dependencies | {check_gfa_job, fasta_job})
        dnascope_dependencies.add(mm2_job)

        # With fastq input, perform dedup and metrics
        sr_alignments: List[pathlib.Path] = []
        cnvscope_deps: Set[Job] = set()
        if self.r1_fastq:
            dnascope_bams.append(out_bwa_aln)
            dnascope_bams.append(out_mm2_aln)

            # Emit Dedup metrics for the primary (bwa) short-read alignment so
            # they land in the metrics directory scanned by MultiQC.
            paths = MetricsPaths.from_output_vcf(ctx.output_vcf)
            dedup_metrics: Optional[pathlib.Path] = None
            if not self.skip_metrics:
                paths.ensure_dir(self.dry_run)
                dedup_metrics = paths.dedup_metrics

            bwa_dedup = DedupStage(
                ctx=ctx,
                tag="bwa",
                inputs=[bwa_bam],
                output=out_bwa_aln,
                score_file=self.tmp_dir.joinpath("sample-bwa-score.txt.gz"),
                dedup_metrics=dedup_metrics,
            ).add_to(dag, bwa_lc_dependencies)
            bwa_dedup_job = bwa_dedup.dedup_job
            dnascope_dependencies.add(bwa_dedup_job)
            mm2_dedup = DedupStage(
                ctx=ctx,
                tag="mm2",
                inputs=[mm2_bam],
                output=out_mm2_aln,
                score_file=self.tmp_dir.joinpath("sample-mm2-score.txt.gz"),
                read_filters=[
                    "IndelLeftAlignReadTransform,rgid="
                    f"{self.fastq_readgroup['ID']}-mm2"
                ],
            ).add_to(dag, {mm2_job})
            dnascope_dependencies.add(mm2_dedup.dedup_job)

            sr_alignments = [out_bwa_aln]
            cnvscope_deps = {bwa_dedup_job}

            if not self.skip_metrics:
                metrics_result = MetricsStage(
                    ctx=ctx,
                    inputs=[out_bwa_aln],
                    algos=self.pangenome_metrics_algos(paths),
                    rehead_metrics=paths.wgs,
                ).add_to(dag, {bwa_dedup_job})
                if not self.skip_multiqc:
                    multiqc_job = self.multiqc()
                    if multiqc_job:
                        dag.add_job(multiqc_job, metrics_result.terminal)
        else:
            dnascope_bams.append(mm2_bam)
            sr_alignments = list(self.sample_input)

        # Stash the short-read alignments for the second DAG
        self.sr_alignments = sr_alignments

        # Estimate the sample ploidy and sex. The JSON output is always
        # written; `--sample_sex` takes precedence for the sex used by
        # the sex-aware callers.
        ploidy_result = PloidyStage(
            ctx=ctx,
            inputs=[self.sr_alignments[0]],
            reference_build=self.reference_build,
        ).add_to(dag, cnvscope_deps)
        self.ploidy_json = ploidy_result.ploidy_json

        # T1K HLA/KIR calling
        self.add_t1k(dag, ctx, self.sr_alignments, cnvscope_deps)

        # DNAscope calling with bwa and mm2 input
        sv_vcf = None
        if self.call_svs:
            sv_vcf = pathlib.Path(
                str(ctx.output_vcf).replace(".vcf.gz", "_sv.vcf.gz")
            )

        if self.skip_small_variants and not self.call_svs:
            return dag

        read_filters: List[str] = []
        if self.tech.upper() == "ULTIMA":
            read_filters.append("UltimaReadFilter")
        pcr_indel_model = "NONE" if self.pcr_free else "CONSERVATIVE"
        model = bundle.joinpath("dnascope.model")
        gfa_file = sample_gfa if self.call_svs else None

        algos: List[BaseAlgo] = []
        if not self.skip_small_variants:
            algos.append(
                DNAscope(
                    raw_vcf,
                    model=model,
                    pcr_indel_model=pcr_indel_model,
                    dbsnp=self.dbsnp,
                    emit_mode="gvcf" if self.gvcf else "variant",
                )
            )
        if sv_vcf and gfa_file:
            algos.append(
                PangenomeSV(
                    sv_vcf,
                    gfa_file=gfa_file,
                    prefix=self.contig_prefix(),
                )
            )
        call = DNAscopeStage(
            ctx=ctx,
            algos=algos,
            inputs=dnascope_bams,
            interval=self.bed,
            read_filter=read_filters,
        ).add_to(dag, dnascope_dependencies)

        if self.skip_small_variants:
            # SV calling only. CNV calling is sex-aware and runs in the
            # second DAG
            return dag

        if self.skip_model_apply and not self.pop_vcf:
            # Nothing post-processes the raw VCF
            return dag

        # When --gvcf is set, the model-apply / transfer outputs are
        # gVCFs; GVCFtyper produces the final VCF at ctx.output_vcf.
        small_variants_out = out_gvcf if self.gvcf else ctx.output_vcf
        snv_apply_vcf = self.tmp_dir.joinpath(
            "sample-snv_apply.g.vcf.gz"
            if self.gvcf
            else "sample-snv_apply.vcf.gz"
        )

        # Transfer annotations from the pop_vcf, then apply the model
        transfer: Optional[TransferSpec] = None
        if self.pop_vcf:
            transfer = TransferSpec(
                config=TransferConfig.from_pipeline(self),
                out_vcf=(
                    transfer_vcf
                    if not self.skip_model_apply
                    else snv_apply_vcf
                ),
            )
        apply_spec: Optional[ApplySpec] = None
        if not self.skip_model_apply:
            apply_spec = ApplySpec(model=model, output=snv_apply_vcf)
        transfer_apply = TransferApplyStage(
            ctx=ctx,
            raw_vcf=raw_vcf,
            transfer=transfer,
            apply=apply_spec,
        ).add_to(dag, call.terminal)

        # Update the overestimated AD/DP of the joint pileup
        ad_update_job = self.build_count_ad_update_job(
            small_variants_out, snv_apply_vcf
        )
        dag.add_job(ad_update_job, transfer_apply.terminal)

        # Genotype the gVCF to also produce a regular VCF at output_vcf
        if self.gvcf:
            GVCFtyperStage(
                ctx=ctx,
                gvcf=out_gvcf,
                output=ctx.output_vcf,
                interval=self.bed,
            ).add_to(dag, {ad_update_job})

        # CNV calling is sex-aware and runs in the second DAG

        return dag

    def bwa_extract_stage(
        self,
        ctx: StageContext,
        sample_bam: pathlib.Path,
        sample_fastq: pathlib.Path,
    ) -> BwaExtractStage:
        """The bwa alignment and read extraction stage.

        The bwa alignment reuses the input readgroup ID with a `-bwa`
        suffix, so it stays distinct from the lifted alignment's.
        """
        bundle = self.required(self.model_bundle, "model_bundle")

        rg = copy.deepcopy(self.fastq_readgroup)
        rg["ID"] = rg["ID"] + "-bwa"
        return BwaExtractStage(
            ctx=ctx,
            output_bam=sample_bam,
            output_fastq=sample_fastq,
            r1_fastq=self.r1_fastq,
            r2_fastq=self.r2_fastq,
            readgroup=(
                "@RG\\t" + "\\t".join([f"{x[0]}:{x[1]}" for x in rg.items()])
            ),
            extract_model=bundle.joinpath(self.extract_model_name),
            bwa_model=bundle.joinpath("bwa.model"),
            unzip=find_unzip(self.logger),
        )

    def build_haplotypes_job(
        self, output_gbz: pathlib.Path, kmer_file: pathlib.Path
    ) -> Job:
        """Build vg haplotypes job"""
        hapl_file = self.required(self.hapl, "hapl")
        gbz_file = self.required(self.gbz, "gbz")

        haplotypes_job = Job(
            cmds.cmd_vg_haplotypes(
                output_gbz,
                kmer_file,
                hapl_file,
                gbz_file,
                threads=self.cores,
                xargs=[
                    "--include-reference",
                    "--diploid-sampling",
                    "--set-reference",
                    self.ref_name(),
                ],
            ),
            "vg-haplotypes",
            self.cores,
            task_name="pangenome",
        )
        return haplotypes_job

    def build_gfa_job(
        self, output_gfa: pathlib.Path, input_gbz: pathlib.Path
    ) -> Job:
        """Build vg convert to GFA job"""
        gfa_job = Job(
            cmds.cmd_vg_convert_gfa(
                output_gfa,
                input_gbz,
                threads=self.cores,
                reference_name=self.ref_name(),
            ),
            "vg-convert-gfa",
            0,
            task_name="pangenome",
        )
        return gfa_job

    def build_fasta_job(
        self, output_fasta: pathlib.Path, input_gbz: pathlib.Path
    ) -> Job:
        """Build vg paths to FASTA job"""
        fasta_job = Job(
            cmds.cmd_vg_paths_fasta(
                output_fasta,
                input_gbz,
            ),
            "vg-paths-fasta",
            0,
            task_name="pangenome",
        )
        return fasta_job

    def build_minimap2_lift_job(
        self,
        mm2_bam: pathlib.Path,
        ext_fastq: pathlib.Path,
        sample_fasta: pathlib.Path,
        sample_gfa: pathlib.Path,
    ) -> Job:
        """Build minimap2 alignment with pgutil lift job"""
        bundle = self.required(self.model_bundle, "model_bundle")
        reference = self.required(self.reference, "reference")

        rg = (
            self.fastq_readgroup
            if self.fastq_readgroup
            else self.bam_readgroups[0]
        )
        rg2 = copy.deepcopy(rg)
        rg2["ID"] = rg2["ID"] + "-mm2"
        rg2["LR"] = "1"

        mm2_model: Union[str, pathlib.Path] = bundle.joinpath("minimap2.model")
        mm2_job = Job(
            cmds.cmd_minimap2_lift(
                mm2_bam,
                sample_fasta,
                ext_fastq,
                sample_gfa,
                reference,
                "@RG\\t" + "\\t".join([f"{x[0]}:{x[1]}" for x in rg2.items()]),
                mm2_model,
                threads=self.cores,
                lift_prefix=self.contig_prefix(),
            ),
            "mm2-lift",
            self.cores,
            task_name="pangenome-alignment",
        )
        return mm2_job

    def build_count_ad_update_job(
        self,
        out_vcf: pathlib.Path,
        in_vcf: pathlib.Path,
    ) -> Job:
        """Update the AD/DP of the joint pileup.

        Reads that are present in both the bwa and the mm2 alignment are
        counted twice, so FORMAT/AD, FORMAT/DP and INFO/DP are inflated.
        The script picks FORMAT/SAD or FORMAT/LAD per sample and derives
        the updated depths from that choice.
        """
        sad_lad_update = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("sad_lad_update.py"))
        )
        return Job(
            cmds.cmd_pyexec_sad_lad_update(
                out_vcf,
                in_vcf,
                sad_lad_update,
                self.cores,
            ),
            "sad-lad-update",
            self.cores,
            task_name="ad-update",
        )

    def build_second_dag(self) -> Optional[DAG]:
        """Build the second DAG for sex-aware downstream tools"""
        if not self._needs_second_dag():
            return None

        assert self.ploidy_json is not None
        self.get_sex(self.ploidy_json)

        self.logger.info("Building the second pangenome DAG")
        dag = DAG()
        ctx = self.stage_context()

        # CNV calling with CNVscope, using the sample sex
        if self._cnv_in_second_dag():
            self.add_pangenome_cnv(
                dag,
                ctx,
                self.output_path("_sv.vcf.gz"),
                self.sr_alignments,
                interval=self.bed,
            )

        if self.sr_alignments:
            self.add_expansion(dag, ctx, self.sr_alignments[0])

            # SegDup calling consumes the small-variant VCF and the
            # inferred sex. segdup-caller's default `main.min_map_qual`
            # of 45 is too strict for Ultima alignments.
            self.add_segdup(
                dag,
                ctx,
                self.sr_alignments[0],
                overrides=(
                    ["main.min_map_qual=30"]
                    if self.tech.upper() == "ULTIMA"
                    else ()
                ),
            )

        return dag
