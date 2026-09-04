"""
Functionality for the DNAscope hybrid pipeline
"""

import argparse
import copy
import json
import pathlib
import sys
from typing import Dict, List, Optional, Set, Tuple

import packaging.version

from importlib.resources import files

from .logging import get_logger
from .archive import ar_load
from . import command_strings as cmds
from .dag import DAG
from .driver import (
    Driver,
    DNAscope,
    HybridStage1,
    HybridStage2,
    HybridStage3,
)
from .job import Job
from .pipeline import BasePipeline
from .shell_pipeline import Command, Pipeline
from .stages.alignment import (
    BWA_FASTQ_MIN_VERSIONS,
    MINIMAP2_REALIGN_MIN_VERSIONS,
    AlignResult,
    BwaFastqStage,
    Minimap2RealignStage,
    find_unzip,
)
from .stages.base import StageContext, driver_job, rm_job
from .stages.cnv import CNV_MIN_VERSIONS, CNVscopeStage
from .stages.metrics import MOSDEPTH_MIN_VERSIONS, MosdepthStage
from .stages.ploidy import PloidyStage
from .stages.preprocessing import (
    ShortReadPreprocessingStage,
    SrPreprocessingResult,
)
from .stages.small_variants import (
    ApplySpec,
    DNAscopeStage,
    TransferApplyStage,
    TransferSpec,
)
from .stages.sv import LONGREADSV_MIN_VERSIONS, LongReadSVStage
from .stages.transfer import TransferConfig
from .util import (
    __version__,
    library_preloaded,
    parse_rg_line,
    path_arg,
    require_versions,
    sample_sex_arg,
    set_bwt_max_mem,
    split_alignment,
    vcf_id,
    versions_available,
)
from .shard import (
    determine_shards_from_fai,
    parse_fai,
    vcf_contigs,
)

logger = get_logger(__name__)


CALLING_MIN_VERSIONS = {
    # 202503.04 writes unsorted BAM to `--hap_bam`, supporting stdout
    "sentieon driver": packaging.version.Version("202503.04"),
    "bedtools": None,
    "bcftools": packaging.version.Version("1.22"),
    "samtools": packaging.version.Version("1.16"),
}

MIN_BUNDLE_VERSION = {
    "ONT": packaging.version.Version("1.2"),
}


class RgInfo:
    """A container class for short and long-read readgroups"""

    def __init__(
        self,
        lr_aln_readgroups: List[List[Dict[str, str]]],
        sr_aln_readgroups: List[List[Dict[str, str]]],
        lr_read_filter: Optional[str],
        sr_read_filter: Optional[str],
        hybrid_set_rg: bool,
        hybrid_rg_sm: str,
        shortread_tech: str,
    ):
        # readgroup handling for long-reads
        self.lr_rg_read_filter: List[str] = []
        self.replace_rg_args: Tuple[List[List[str]], List[List[str]]] = (
            [],
            [],
        )
        for aln_rgs in lr_aln_readgroups:
            self.replace_rg_args[0].append([])
            for rg_line_d in aln_rgs:
                id = rg_line_d.get("ID")
                new_sm = hybrid_rg_sm if hybrid_set_rg else rg_line_d.get("SM")
                self.replace_rg_args[0][-1].append(
                    f"{id}=ID:{id}\\tSM:{new_sm}\\tLR:1"
                )
                if lr_read_filter:
                    self.lr_rg_read_filter.append(
                        f"{lr_read_filter},rgid={id}"
                    )

        # readgroup handling for short-reads
        self.sr_rg_read_filter: List[str] = []
        self.ultima_read_filter: List[str] = []
        for aln_rgs in sr_aln_readgroups:
            if hybrid_set_rg:
                self.replace_rg_args[1].append([])
            for rg_line_d in aln_rgs:
                id = rg_line_d.get("ID")
                new_sm = hybrid_rg_sm if hybrid_set_rg else rg_line_d.get("SM")
                if hybrid_set_rg:
                    self.replace_rg_args[1][-1].append(
                        f"{id}=ID:{id}\\tSM:{new_sm}"
                    )
                if sr_read_filter:
                    self.sr_rg_read_filter.append(
                        f"{sr_read_filter},rgid={id}"
                    )
                if shortread_tech.upper() == "ULTIMA":
                    self.ultima_read_filter.append(
                        f"UltimaReadFilter,rgid={id}"
                    )


class DNAscopeHybridPipeline(BasePipeline):
    """The DNAscope Hybrid pipeline"""

    params = copy.deepcopy(BasePipeline.params)
    params.update(
        {
            # Required arguments
            "lr_aln": {
                "nargs": "*",
                "help": "Long-read BAM or CRAM files.",
                "type": path_arg(exists=True, is_file=True),
                "required": True,
            },
            "model_bundle": {
                "flags": ["-m", "--model_bundle"],
                "help": "The model bundle file.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "sr_aln": {
                "nargs": "*",
                "help": "Short-read BAM or CRAM files",
                "type": path_arg(exists=True, is_file=True),
            },
            "sr_r1_fastq": {
                "nargs": "*",
                "help": "Short-read R1 fastq files",
                "type": path_arg(exists=True, is_file=True),
            },
            "sr_r2_fastq": {
                "nargs": "*",
                "help": "Short-read R2 fastq files",
                "type": path_arg(exists=True, is_file=True),
            },
            "sr_readgroups": {
                "nargs": "*",
                "help": "Readgroup information for the short-read fastq files",
            },
            # Additional arguments
            "bam_format": {
                "help": (
                    "Use the BAM format instead of CRAM for output aligned "
                    "files."
                ),
                "action": "store_true",
            },
            "bed": {
                "flags": ["-b", "--bed"],
                "help": (
                    "Region BED file. Supplying this file will limit variant "
                    "calling to the intervals inside the BED file."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "dbsnp": {
                "flags": ["-d", "--dbsnp"],
                "help": (
                    "dbSNP vcf file Supplying this file will annotate "
                    "variants with their dbSNP refSNP ID numbers."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "gvcf": {
                "flags": ["-g", "--gvcf"],
                "help": (
                    "Generate a gVCF output file along with the VCF."
                    " (default generates only the VCF)"
                ),
                "action": "store_true",
            },
            "lr_align_input": {
                "help": (
                    "Align the input long-read BAM/CRAM/uBAM file to the "
                    "reference genome"
                ),
                "action": "store_true",
            },
            "lr_input_ref": {
                "help": (
                    "Used to decode the input long-read alignment file. "
                    "Required if the input file is in the CRAM/uCRAM formats."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "par_bed": {
                "help": (
                    "A BED file of the pseudo-autosomal regions (PAR), used "
                    "for sex-aware CNV calling of male samples. Overrides the "
                    "PAR BED file selected for the reference genome."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "pop_vcf": {
                "flags": ["--pop_vcf"],
                "help": (
                    "A VCF containing annotations for use with DNAModelApply."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "sample_sex": {
                "help": (
                    "The sample sex, used for sex-aware CNV calling. "
                    "Supplying this argument overrides the sex estimated "
                    "from read coverage."
                ),
                "metavar": "{male,female}",
                "type": sample_sex_arg,
            },
            "rgsm": {
                "help": (
                    "Overwrite the SM tag of the input readgroups for "
                    "compatibility"
                )
            },
            "skip_cnv": {
                "help": "Skip CNV calling.",
                "action": "store_true",
            },
            "skip_metrics": {
                "help": "Skip all metrics collection and multiQC",
                "action": "store_true",
            },
            "skip_mosdepth": {
                "help": "Skip QC with mosdepth.",
                "action": "store_true",
            },
            "skip_multiqc": {
                "help": "Skip multiQC report generation.",
                "action": "store_true",
            },
            "skip_svs": {
                "help": "Skip SV calling",
                "action": "store_true",
            },
            "sr_duplicate_marking": {
                "help": "Options for duplicate marking.",
                "choices": ["markdup", "rmdup", "none"],
                "default": "markdup",
            },
            # Hidden arguments
            "bwa_args": {
                # help="Extra arguments for sentieon bwa",
                "help": argparse.SUPPRESS,
                "default": "",
            },
            "bwa_k": {
                # help="The '-K' argument in bwa",
                "help": argparse.SUPPRESS,
                "default": 100000000,
            },
            "bwt_max_mem": {
                # Manually set `bwt_max_mem`
                "help": argparse.SUPPRESS,
            },
            "lr_fastq_taglist": {
                # help="A comma-separated list of tags to retain. Defaults to "
                # "'%(default)s' and the 'RG' tag is required",
                "help": argparse.SUPPRESS,
                "default": "*",
            },
            "lr_read_filter": {
                "help": argparse.SUPPRESS,
            },
            "minimap2_args": {
                # help="Extra arguments for sentieon minimap2",
                "help": argparse.SUPPRESS,
                "default": "-YL",
            },
            "no_split_alignment": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_model_apply": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_pop_vcf_id_check": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "sr_read_filter": {
                "help": argparse.SUPPRESS,
            },
            "util_sort_args": {
                # help="Extra arguments for sentieon util sort",
                "help": argparse.SUPPRESS,
                "default": "--cram_write_options version=3.0,compressor=rans",
            },
        }
    )

    def __init__(self) -> None:
        super().__init__()
        # Required arguments
        self.lr_aln: List[pathlib.Path] = []
        self.model_bundle: Optional[pathlib.Path] = None
        self.sr_aln: List[pathlib.Path] = []
        self.sr_r1_fastq: List[pathlib.Path] = []
        self.sr_r2_fastq: List[pathlib.Path] = []
        self.sr_readgroups: List[str] = []
        # Additional arguments
        self.bam_format = False
        self.bed: Optional[pathlib.Path] = None
        self.dbsnp: Optional[pathlib.Path] = None
        self.gvcf = False
        self.lr_align_input = False
        self.lr_input_ref: Optional[pathlib.Path] = None
        self.par_bed: Optional[pathlib.Path] = None
        self.pop_vcf: Optional[pathlib.Path] = None
        self.rgsm: Optional[str] = None
        self.skip_cnv = False
        self.skip_metrics = False
        self.skip_mosdepth = False
        self.skip_multiqc = False
        self.skip_svs = False
        self.sr_duplicate_marking = "markdup"
        # Hidden arguments
        self.bwa_args = ""
        self.bwa_k = 100000000
        self.bwt_max_mem: Optional[str] = None
        self.lr_fastq_taglist = "*"
        self.lr_read_filter: Optional[str] = None
        self.minimap2_args = "-YL"
        self.no_split_alignment = False
        self.skip_model_apply = False
        self.skip_pop_vcf_id_check = False
        self.sr_read_filter: Optional[str] = None
        self.util_sort_args = (
            "--cram_write_options version=3.0,compressor=rans"
        )
        # Stashed by `build_dag` for the second, sex-aware DAG
        self.ploidy_json: Optional[pathlib.Path] = None
        self.cnv_sr_aln: List[pathlib.Path] = []
        self.cnv_replace_rg: Optional[List[List[str]]] = None

    def validate(self) -> None:
        self.fai_data = parse_fai(pathlib.Path(str(self.reference) + ".fai"))
        self.pop_vcf_contigs: Dict[str, Optional[int]] = {}
        if self.pop_vcf:
            self.pop_vcf_contigs = vcf_contigs(self.pop_vcf, self.dry_run)
            self.logger.debug("VCF contigs are: %s", self.pop_vcf_contigs)
        self.shards = determine_shards_from_fai(
            self.fai_data, 10 * 1000 * 1000
        )

        self.validate_bundle()
        self.collect_readgroups()
        self.validate_readgroups()
        # `validate_bundle` may have set `skip_cnv`
        self.validate_cnv()

        self.validate_output_vcf()
        if not self.sr_aln and not self.sr_r1_fastq:
            self.logger.error("Please supply a short-read input file")
            sys.exit(2)

        self.skip_multiqc = True if self.skip_metrics else self.skip_multiqc
        self.skip_mosdepth = True if self.skip_metrics else self.skip_mosdepth

        if not library_preloaded("libjemalloc.so"):
            self.logger.warning(
                "jemalloc is recommended, but is not preloaded. See "
                "https://support.sentieon.com/appnotes/jemalloc/"
            )

        if self.bed is None:
            self.logger.info(
                "A BED file is recommended to avoid variant calling "
                "across decoy and unplaced contigs."
            )

        if len(self.sr_r1_fastq) != len(self.sr_readgroups):
            self.logger.error(
                "The number of readgroups does not equal the number of fastq "
                "files"
            )
            sys.exit(2)

        if self.sr_r1_fastq:
            self.validate_bwa_index()

    def validate_cnv(self) -> None:
        """Validate the arguments used for sex-aware CNV calling"""
        self.resolve_cnv_par_bed(
            self.fai_data, self.par_bed, not self.skip_cnv
        )
        if self.skip_cnv:
            return

        # A PAR BED file is required for CNV calling, whatever the sex
        self.validate_cnv_par(True)

        # CNV calling now runs in the second DAG, so the version of the
        # driver is checked before the first DAG starts
        require_versions(CNV_MIN_VERSIONS, skip=self.skip_version_check)

    def validate_bundle(self) -> None:
        bundle_info_bytes = ar_load(
            str(self.model_bundle) + "/bundle_info.json"
        )
        if isinstance(bundle_info_bytes, list):
            bundle_info_bytes = b"{}"

        bundle_info = json.loads(bundle_info_bytes.decode())
        self.longread_tech = bundle_info.get("longReadPlatform")
        self.shortread_tech = bundle_info.get("shortReadPlatform")
        bundle_vcf_id = bundle_info.get("SentieonVcfID")

        if bundle_vcf_id:
            if not self.pop_vcf:
                self.logger.error(
                    "The model bundle requires a population VCF file. Please "
                    "supply the `--pop_vcf` argument."
                )
                sys.exit(2)
            if not self.skip_pop_vcf_id_check and not self.dry_run:
                pop_vcf_id = vcf_id(self.pop_vcf)
                if bundle_vcf_id != pop_vcf_id:
                    self.logger.error(
                        "The ID of the `--pop_vcf` does not match the model "
                        "bundle"
                    )
                    sys.exit(2)

        if not self.longread_tech or not self.shortread_tech:
            self.logger.error(
                "The bundle file does not have the expected attributes. "
                "Please check that you using the latest bundle version."
            )
            sys.exit(2)
        if not self.shortread_tech:
            self.shortread_tech = "Illumina"
        req_version = packaging.version.Version(
            bundle_info.get("minScriptVersion", __version__)
        )
        if req_version > packaging.version.Version(__version__):
            self.logger.error(
                "The model bundle requires version %s or later of the "
                "sentieon-cli.",
                req_version,
            )
            sys.exit(2)
        if bundle_info.get("pipeline", "DNAscope Hybrid") != "DNAscope Hybrid":
            self.logger.error("The model bundle is for a different pipeline.")
            sys.exit(2)
        bundle_version = packaging.version.Version(
            bundle_info.get("bundleVersion", "1.0")
        )
        if self.longread_tech.upper() == "ONT":
            if bundle_version < MIN_BUNDLE_VERSION["ONT"]:
                self.logger.error(
                    "The model bundle is for an older version of the "
                    "pipeline. Please update to the latest model bundle."
                )
                sys.exit(2)

        bundle_members = set(ar_load(str(self.model_bundle)))
        if "longreadsv.model" not in bundle_members:
            self.logger.info("No LongReadSV model found. Skipping SV calling")
            self.skip_svs = True
        if "cnv.model" not in bundle_members:
            self.logger.info("No CNVscope model found. Skipping CNV calling")
            self.skip_cnv = True
        if "bwa.model" not in bundle_members and len(self.sr_r1_fastq) > 0:
            self.logger.error(
                "Alignment with bwa is not supported with this model bundle"
            )
            sys.exit(2)

    def collect_readgroups(self) -> None:
        # Collect the readgroup tags from all input files
        self.all_readgroups: Tuple[
            List[List[Dict[str, str]]],  # long-read alignments
            List[List[Dict[str, str]]],  # short-read alignments
            List[List[Dict[str, str]]],  # short-read fastq
        ] = ([], [], [])
        for i, aln_list in enumerate((self.lr_aln, self.sr_aln)):
            for aln in aln_list:
                self.all_readgroups[i].append([])
                aln_rgs = cmds.get_rg_lines(aln, self.dry_run)
                for rg_line in aln_rgs:
                    self.all_readgroups[i][-1].append(parse_rg_line(rg_line))
        for rg in self.sr_readgroups:
            self.all_readgroups[2].append([])
            try:
                parsed = parse_rg_line(rg.replace(r"\t", "\t"))
            except ValueError as e:
                self.logger.error(
                    "Invalid --sr_readgroups value '%s': %s", rg, e
                )
                sys.exit(2)
            self.all_readgroups[2][-1].append(parsed)

    def validate_readgroups(self) -> None:
        # Confirm that all readgroups have the same RGSM
        self.hybrid_rg_sm = ""
        rg_sm_tag = None
        for rg_aln_list in self.all_readgroups:
            for rg_list in rg_aln_list:
                for aln_rg in rg_list:
                    sm = aln_rg.get("SM")
                    if not sm:
                        self.logger.error(
                            "Found a readgroup without a SM tag: %s",
                            str(aln_rg),
                        )
                        sys.exit(2)
                    if not aln_rg.get("ID"):
                        self.logger.error(
                            "Found a readgroup without an ID tag: %s",
                            str(aln_rg),
                        )

                    if rg_sm_tag is None:
                        rg_sm_tag = sm
                        self.hybrid_rg_sm = sm
                    if self.dry_run:
                        continue
                    elif rg_sm_tag != sm and not self.rgsm:
                        self.logger.error(
                            "Input readgroup '%s' has a different RG-SM tag"
                            " from previously seen alignment files.\n"
                            "found='%s' expected='%s'. Please set the `--rgsm`"
                            " argument to override the SM tag in the input"
                            " files",
                            str(aln_rg),
                            sm,
                            rg_sm_tag,
                        )
                        sys.exit(2)
        self.hybrid_set_rg = True if self.rgsm else False
        if self.rgsm:
            self.hybrid_rg_sm = self.rgsm

    def configure(self) -> None:
        self.configure_readgroups()
        self.numa_nodes: List[str] = []
        n_alignment_jobs = 1
        if not self.no_split_alignment:
            self.numa_nodes = split_alignment(self.cores)
        n_alignment_jobs = max(1, len(self.numa_nodes))

        set_bwt_max_mem(0, n_alignment_jobs, override=self.bwt_max_mem)

    def configure_readgroups(self) -> None:
        self.lr_aln_readgroups = self.all_readgroups[0]
        self.sr_aln_readgroups: List[List[Dict[str, str]]] = []
        if self.sr_duplicate_marking == "none":
            # Without dedup, retain the original rg structure
            self.sr_aln_readgroups = copy.deepcopy(self.all_readgroups[1])
            self.sr_aln_readgroups.extend(
                copy.deepcopy(self.all_readgroups[2])
            )
        else:
            # Flatten
            self.sr_aln_readgroups.append([])
            for rg_aln_list in self.all_readgroups[1:]:
                for rg_list in rg_aln_list:
                    for aln_rg in rg_list:
                        self.sr_aln_readgroups[-1].append(aln_rg)

    def build_dag(self) -> DAG:
        """Build the DAG for the pipeline"""
        self.logger.info("Building the DAG")
        dag = DAG()

        self.required(self.model_bundle, "model_bundle")
        ctx = self.stage_context()

        rg_info = RgInfo(
            self.lr_aln_readgroups,
            self.sr_aln_readgroups,
            self.lr_read_filter,
            self.sr_read_filter,
            self.hybrid_set_rg,
            self.hybrid_rg_sm,
            self.shortread_tech,
        )

        # Short-read alignment
        sr_aln = self.sr_aln
        fq_result = self.add_sr_fastq_alignment(dag, ctx)
        align_fastq_jobs = set(fq_result.jobs)
        sr_aln.extend(fq_result.outputs)

        # Short-read dedup
        preprocessing = self.add_sr_preprocessing(
            dag, ctx, sr_aln, set(align_fastq_jobs)
        )
        sr_aln = preprocessing.deduped
        dedup_job = preprocessing.dedup_job
        # The fastq alignment has cleanup paths only when it is an
        # intermediate, which is exactly when duplicate marking -- and so
        # `dedup_job` -- follows it
        if dedup_job and fq_result.cleanup_paths:
            dag.add_job(
                rm_job(fq_result.cleanup_paths, "rm-fq-aln"), {dedup_job}
            )

        # Long-read alignment
        realign_jobs: Set[Job] = set()
        lr_aln = self.lr_aln
        if self.lr_align_input:
            realign_result = self.add_lr_realignment(dag, ctx)
            lr_aln = realign_result.outputs
            realign_jobs = set(realign_result.jobs)

        if not self.skip_mosdepth:
            self.add_mosdepth(dag, ctx, lr_aln, realign_jobs)

        if not self.skip_multiqc:
            multiqc_job = self.multiqc()
            if multiqc_job:
                dag.add_job(multiqc_job, set(preprocessing.qc_jobs))

        if not self.skip_svs:
            self.add_sv_calling(
                dag,
                ctx,
                lr_aln,
                rg_info.replace_rg_args[0],
                realign_jobs,
            )

        sr_preprocessing_jobs: Set[Job] = set()
        sr_preprocessing_jobs.update(align_fastq_jobs)
        if dedup_job:
            sr_preprocessing_jobs.add(dedup_job)
        # Estimate the sample ploidy and sex. The JSON output is always
        # written; `--sample_sex` takes precedence for the sex used by
        # sex-aware CNV calling.
        ploidy_result = PloidyStage(
            ctx=ctx,
            inputs=sr_aln,
            reference_build=self.reference_build,
        ).add_to(dag, sr_preprocessing_jobs)
        self.ploidy_json = ploidy_result.ploidy_json

        if not self.skip_cnv:
            # CNV calling is sex-aware and runs in the second DAG
            self.cnv_sr_aln = sr_aln
            self.cnv_replace_rg = (
                rg_info.replace_rg_args[1]
                if rg_info.replace_rg_args[1]
                else None
            )

        self.add_variant_calling(
            dag,
            ctx,
            sr_aln,
            lr_aln,
            rg_info,
            realign_jobs | sr_preprocessing_jobs,
        )

        return dag

    def add_sr_fastq_alignment(
        self, dag: DAG, ctx: StageContext
    ) -> AlignResult:
        """Align the short-read fastq files with bwa"""
        bundle = self.required(self.model_bundle, "model_bundle")

        if not self.sr_r1_fastq and not self.sr_readgroups:
            return AlignResult(
                jobs=[], terminal=set(), outputs=[], cleanup_paths=[]
            )

        require_versions(BWA_FASTQ_MIN_VERSIONS, skip=self.skip_version_check)

        return BwaFastqStage(
            ctx=ctx,
            r1_fastq=self.sr_r1_fastq,
            r2_fastq=self.sr_r2_fastq,
            readgroups=self.sr_readgroups,
            model_bundle=bundle,
            numa_nodes=self.numa_nodes,
            bam_format=self.bam_format,
            duplicate_marking=self.sr_duplicate_marking,
            unzip=find_unzip(self.logger),
            bwa_args=self.bwa_args,
            bwa_k=str(self.bwa_k),
            util_sort_args=self.util_sort_args,
        ).add_to(dag)

    def add_sr_preprocessing(
        self,
        dag: DAG,
        ctx: StageContext,
        sr_aln: List[pathlib.Path],
        upstream: Set[Job],
    ) -> SrPreprocessingResult:
        """Mark duplicates and collect metrics from the short reads"""
        return ShortReadPreprocessingStage(
            ctx=ctx,
            inputs=sr_aln,
            duplicate_marking=self.sr_duplicate_marking,
            consensus=False,
            skip_metrics=self.skip_metrics,
            assay="WGS",
            bed=self.bed,
            bam_format=self.bam_format,
        ).add_to(dag, upstream)

    def add_lr_realignment(self, dag: DAG, ctx: StageContext) -> AlignResult:
        """Re-align the long-read input files with minimap2"""
        bundle = self.required(self.model_bundle, "model_bundle")
        require_versions(
            MINIMAP2_REALIGN_MIN_VERSIONS, skip=self.skip_version_check
        )

        return Minimap2RealignStage(
            ctx=ctx,
            inputs=self.lr_aln,
            model_bundle=bundle,
            sample_name=ctx.output_vcf.name.replace(".vcf.gz", ""),
            bam_format=self.bam_format,
            input_ref=self.lr_input_ref,
            fastq_taglist=self.lr_fastq_taglist,
            minimap2_args=self.minimap2_args,
            util_sort_args=self.util_sort_args,
        ).add_to(dag)

    def add_mosdepth(
        self,
        dag: DAG,
        ctx: StageContext,
        inputs: List[pathlib.Path],
        upstream: Set[Job],
    ) -> None:
        """Run mosdepth for QC, when it is available"""
        if not versions_available(
            MOSDEPTH_MIN_VERSIONS, skip=self.skip_version_check
        ):
            self.logger.warning(
                "Skipping mosdepth. mosdepth version %s or later not found",
                MOSDEPTH_MIN_VERSIONS["mosdepth"],
            )
            return

        MosdepthStage(ctx=ctx, inputs=inputs).add_to(dag, upstream)

    def add_sv_calling(
        self,
        dag: DAG,
        ctx: StageContext,
        inputs: List[pathlib.Path],
        replace_rg: Optional[List[List[str]]],
        upstream: Set[Job],
    ) -> None:
        """Call SVs from the long reads with Sentieon LongReadSV"""
        bundle = self.required(self.model_bundle, "model_bundle")
        require_versions(LONGREADSV_MIN_VERSIONS, skip=self.skip_version_check)

        LongReadSVStage(
            ctx=ctx,
            inputs=inputs,
            model=bundle.joinpath("longreadsv.model"),
            interval=self.bed,
            replace_rg=replace_rg,
        ).add_to(dag, upstream)

    def build_second_dag(self) -> Optional[DAG]:
        """Build the second DAG, for sex-aware CNV calling"""
        if self.skip_cnv:
            return None
        bundle = self.required(self.model_bundle, "model_bundle")
        assert self.ploidy_json is not None
        ctx = self.stage_context()

        self.get_sex(self.ploidy_json)

        self.logger.info("Building the DNAscope Hybrid CNV DAG")
        dag = DAG()
        CNVscopeStage(
            ctx=ctx,
            inputs=self.cnv_sr_aln,
            model=bundle.joinpath("cnv.model"),
            cnvscope_vcf=ctx.tmp_dir.joinpath("cnvscope.vcf.gz"),
            cnv_vcf=pathlib.Path(
                str(ctx.output_vcf).replace(".vcf.gz", ".cnv.vcf.gz")
            ),
            sample_sex=self.sample_sex,
            par_bed=self.cnv_par_bed,
            interval=self.bed,
            replace_rg=self.cnv_replace_rg,
        ).add_to(dag)
        return dag

    def transfer_config(self) -> TransferConfig:
        """The inputs the annotation-transfer stage needs.

        Only valid once `--pop_vcf` is known to be set; the calling code
        checks it first.
        """
        assert self.pop_vcf is not None
        return TransferConfig(
            pop_vcf=self.pop_vcf,
            shards=self.shards,
            pop_vcf_contigs=self.pop_vcf_contigs,
            fai_data=self.fai_data,
        )

    def add_cleanup(
        self,
        dag: DAG,
        paths: List[pathlib.Path],
        name: str,
        upstream: Set[Job],
    ) -> None:
        """Remove intermediate files once `upstream` has finished.

        A no-op under `--retain_tmpdir`, which keeps everything the run
        wrote in the temporary directory.
        """
        if self.retain_tmpdir:
            return
        dag.add_job(rm_job(paths, name), upstream)

    def add_variant_calling(
        self,
        dag: DAG,
        ctx: StageContext,
        sr_aln: List[pathlib.Path],
        lr_aln: List[pathlib.Path],
        rg_info: RgInfo,
        upstream: Set[Job],
    ) -> None:
        """
        Call SNVs and indels using the DNAscope hybrid pipeline

        A first calling pass over the input alignments picks the regions
        that need realignment; the three hybrid-realignment stages realign
        the short reads there; a second pass recalls those regions, and the
        two call sets are merged, annotated and filtered.
        """
        bundle = self.required(self.model_bundle, "model_bundle")

        ref_fai = pathlib.Path(str(ctx.reference) + ".fai")

        require_versions(CALLING_MIN_VERSIONS, skip=self.skip_version_check)

        hybrid_model = bundle.joinpath("hybrid.model")

        # First pass - combined variant calling
        vcf_suffix = ".g.vcf.gz" if self.gvcf else ".vcf.gz"
        combined_vcf = self.tmp_dir.joinpath("initial" + vcf_suffix)
        emit_mode = "gvcf" if self.gvcf else None
        call = DNAscopeStage(
            ctx=ctx,
            algos=[
                DNAscope(
                    combined_vcf,
                    dbsnp=self.dbsnp,
                    model=hybrid_model,
                    pcr_indel_model="none",
                    emit_mode=emit_mode,
                )
            ],
            inputs=lr_aln + sr_aln,
            interval=self.bed,
            read_filter=rg_info.ultima_read_filter
            + rg_info.lr_rg_read_filter
            + rg_info.sr_rg_read_filter,
            replace_rg=rg_info.replace_rg_args[0] + rg_info.replace_rg_args[1],
            name="dnascope-1",
        ).add_to(dag, upstream)

        # Region selection
        hybrid_select = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("hybrid_select.py"))
        ).resolve()
        selected_bed = self.tmp_dir.joinpath("selected.bed")
        select_job = Job(
            cmds.cmd_pyexec_hybrid_select(
                out_bed=selected_bed,
                vcf=combined_vcf,
                ref_fai=ref_fai,
                hybrid_select=hybrid_select,
                threads=self.cores,
            ),
            "hybrid-select",
            0,
            task_name="region-selection",
        )
        dag.add_job(select_job, call.terminal)

        mapq0_bed = self.tmp_dir.joinpath("hybrid_mapq0.bed")
        mapq0_job = driver_job(
            ctx,
            [
                HybridStage2(
                    model=bundle.joinpath("HybridStage2_region.model"),
                    all_bed=mapq0_bed,
                )
            ],
            name="mapq0-bed",
            task_name="region-selection",
            inputs=lr_aln + sr_aln,
            read_filter=rg_info.lr_rg_read_filter + rg_info.sr_rg_read_filter,
            replace_rg=rg_info.replace_rg_args[0] + rg_info.replace_rg_args[1],
        )
        dag.add_job(mapq0_job, upstream)

        mapq0_slop_bed = self.tmp_dir.joinpath("hybrid_mapq0.ex1000.bed")
        mapq0_slop_job = Job(
            cmds.cmd_bedtools_slop(
                mapq0_bed,
                mapq0_slop_bed,
                ref_fai,
                1000,
            ),
            "mapq0-bed-slop",
            0,
            task_name="region-selection",
        )
        dag.add_job(mapq0_slop_job, {mapq0_job})

        diff_bed = self.tmp_dir.joinpath("merged_diff.bed")
        cat_merge_job = Job(
            cmds.cmd_bedtools_cat_sort_merge(
                out_bed=diff_bed,
                in_bed=[selected_bed, mapq0_slop_bed],
                ref_fai=ref_fai,
            ),
            "concat-merge-bed",
            0,
            task_name="region-selection",
        )
        dag.add_job(cat_merge_job, {mapq0_slop_job, select_job})
        self.add_cleanup(
            dag,
            [selected_bed, mapq0_slop_bed],
            "rm-tmp1",
            {cat_merge_job},
        )

        stage1_ins_fa = self.tmp_dir.joinpath("stage1_ins.fa")
        stage1_ins_bed = self.tmp_dir.joinpath("stage1_ins.bed")
        ins_driver = Driver(
            reference=ctx.reference,
            thread_count=self.cores,
            replace_rg=rg_info.replace_rg_args[0],
            input=lr_aln,
            read_filter=rg_info.lr_rg_read_filter,
        )
        ins_driver.add_algo(
            HybridStage1(
                "-",
                model=bundle.joinpath("HybridStage1_ins.model"),
                fa_file=stage1_ins_fa,
                bed_file=stage1_ins_bed,
            )
        )

        stage1_hap_bam = self.tmp_dir.joinpath("stage1_hap.bam")
        stage1_hap_bed = self.tmp_dir.joinpath("stage1_hap.bed")
        stage1_hap_vcf = self.tmp_dir.joinpath("stage1_hap.vcf")

        # The algo writes its fastq output to a fifo read by the bwa job
        # and its unsorted haplotype alignments to stdout
        stage1_fifo = self.tmp_dir.joinpath("stage1_hap.fq")
        stage1_fifo_job = Job(
            Pipeline(Command("mkfifo", str(stage1_fifo))),
            "stage1-fifo",
            1,
            task_name="hybrid-realignment",
        )
        # The haplotype job writes the fifo read by the first-stage job
        dag.add_job(stage1_fifo_job)

        stage1_driver = Driver(
            reference=ctx.reference,
            thread_count=self.cores,
            replace_rg=rg_info.replace_rg_args[0],
            input=lr_aln,
            interval=diff_bed,
            read_filter=rg_info.lr_rg_read_filter,
        )
        stage1_driver.add_algo(
            HybridStage1(
                stage1_fifo,
                model=bundle.joinpath("HybridStage1.model"),
                hap_bam="-",
                hap_bed=stage1_hap_bed,
                hap_vcf=stage1_hap_vcf,
            )
        )
        stage1_hap_job = Job(
            cmds.hybrid_stage1_hap(
                stage1_hap_bam,
                stage1_driver,
                self.cores,
            ),
            "first-stage-hap",
            # This job writes the fifo that the `first-stage` job reads, so
            # the two need to run concurrently. Jobs requesting 0 threads
            # start immediately instead of waiting for the thread budget;
            # requesting `self.cores` here could let the scheduler serialize
            # the two jobs and deadlock on the fifo.
            0,
            task_name="hybrid-realignment",
        )
        dag.add_job(stage1_hap_job, {stage1_fifo_job, cat_merge_job})

        stage1_bam = self.tmp_dir.joinpath("hybrid_stage1.bam")
        stage1_job = Job(
            cmds.hybrid_stage1(
                out_aln=stage1_bam,
                reference=ctx.reference,
                cores=self.cores,
                readgroup=f"@RG\\tID:hybrid-18893\\tSM:{self.hybrid_rg_sm}",
                ins_driver=ins_driver,
                hap_fastq_fifo=stage1_fifo,
                bwa_model=bundle.joinpath("HybridStage1_bwa.model"),
            ),
            "first-stage",
            self.cores,
            task_name="hybrid-realignment",
        )
        dag.add_job(stage1_job, {stage1_fifo_job, cat_merge_job})
        self.add_cleanup(
            dag,
            [stage1_ins_fa, stage1_ins_bed, stage1_hap_vcf],
            "rm-tmp2",
            {stage1_job, stage1_hap_job},
        )

        stage2_bed = self.tmp_dir.joinpath("hybrid_stage2.bed")
        stage2_unmap_bam = self.tmp_dir.joinpath("hybrid_stage2_unmap.bam")
        stage2_alt_bam = self.tmp_dir.joinpath("hybrid_stage2_alt.bam")
        second_stage_job = driver_job(
            ctx,
            [
                HybridStage2(
                    model=bundle.joinpath("HybridStage2.model"),
                    hap_bed=stage1_hap_bed,
                    unmap_bam=stage2_unmap_bam,
                    alt_bam=stage2_alt_bam,
                    all_bed=stage2_bed,
                )
            ],
            name="second-stage",
            task_name="hybrid-realignment",
            inputs=[stage1_bam, stage1_hap_bam],
        )
        dag.add_job(second_stage_job, {stage1_job, stage1_hap_job})
        self.add_cleanup(
            dag,
            [stage1_bam, stage1_hap_bam],
            "rm-tmp3",
            {second_stage_job},
        )

        suffix = "bam" if self.bam_format else "cram"
        stage3_aln = pathlib.Path(
            str(ctx.output_vcf).replace(".vcf.gz", f"_sr_realigned.{suffix}")
        )
        stage3_driver = Driver(
            reference=ctx.reference,
            thread_count=self.cores,
            replace_rg=rg_info.replace_rg_args[0] + rg_info.replace_rg_args[1],
            input=lr_aln + sr_aln + [stage2_unmap_bam, stage2_alt_bam],
            interval=stage2_bed,
            read_filter=rg_info.lr_rg_read_filter + rg_info.sr_rg_read_filter,
        )
        stage3_driver.add_algo(
            HybridStage3(
                "-",
                model=bundle.joinpath("HybridStage3.model"),
            )
        )
        third_stage_job = Job(
            cmds.hybrid_stage3(
                stage3_aln,
                reference=ctx.reference,
                driver=stage3_driver,
                cores=self.cores,
            ),
            "third-stage",
            self.cores,
            task_name="hybrid-realignment",
        )
        dag.add_job(third_stage_job, {second_stage_job})
        self.add_cleanup(
            dag,
            [stage2_unmap_bam, stage2_alt_bam],
            "rm-tmp4",
            {third_stage_job},
        )

        # pass 2 of variant calling
        pass2_vcf = self.tmp_dir.joinpath("hybrid_pass2.vcf.gz")
        call2 = DNAscopeStage(
            ctx=ctx,
            algos=[
                DNAscope(
                    pass2_vcf,
                    dbsnp=self.dbsnp,
                    model=hybrid_model,
                    pcr_indel_model="none",
                    emit_mode=emit_mode,
                )
            ],
            inputs=lr_aln + [stage3_aln],
            interval=stage2_bed,
            read_filter=rg_info.ultima_read_filter + rg_info.lr_rg_read_filter,
            replace_rg=rg_info.replace_rg_args[0],
            name="dnascope-2",
        ).add_to(dag, {third_stage_job})

        # Merge and normalize the VCFs
        subset_vcf = self.tmp_dir.joinpath("mix_subset.vcf.gz")
        combined_tmp_vcf = self.tmp_dir.joinpath("combined_tmp.vcf.gz")
        subset_job = Job(
            cmds.bcftools_subset(
                subset_vcf,
                combined_vcf,
                stage2_bed,
            ),
            "subset-calls",
            0,
            task_name="vcf-merge",
        )
        dag.add_job(subset_job, {second_stage_job})
        concat_job = Job(
            cmds.bcftools_concat(
                combined_tmp_vcf,
                [subset_vcf, pass2_vcf],
            ),
            "concat-calls",
            0,
            task_name="vcf-merge",
        )
        dag.add_job(concat_job, {subset_job} | call2.terminal)
        self.add_cleanup(
            dag,
            [combined_vcf, subset_vcf, pass2_vcf],
            "rm-tmp5",
            {concat_job},
        )

        # Annotate the output VCF
        hybrid_anno = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("hybrid_anno.py"))
        )
        anno_target = self.tmp_dir.joinpath("combined_tmp_anno.vcf.gz")
        if self.skip_model_apply and not self.pop_vcf:
            anno_target = ctx.output_vcf

        anno_job = Job(
            cmds.cmd_pyexec_hybrid_anno(
                anno_target,
                combined_tmp_vcf,
                stage1_hap_bed,
                hybrid_anno,
                self.cores,
            ),
            "anno-calls",
            0,
            task_name="annotation",
        )
        dag.add_job(anno_job, {concat_job})

        if not self.pop_vcf and self.skip_model_apply:
            # `anno_target` is already the output VCF
            return

        # Transfer annotations and apply the model
        transfer_spec: Optional[TransferSpec] = None
        if self.pop_vcf:
            transfer_target = self.tmp_dir.joinpath(
                "combined_tmp_transfer.vcf.gz"
            )
            if self.skip_model_apply:
                transfer_target = ctx.output_vcf
            transfer_spec = TransferSpec(
                self.transfer_config(), transfer_target
            )

        apply_spec: Optional[ApplySpec] = None
        apply_vcf = self.tmp_dir.joinpath("combined_apply.vcf.gz")
        if not self.skip_model_apply:
            apply_spec = ApplySpec(hybrid_model, apply_vcf)

        transfer_apply = TransferApplyStage(
            ctx=ctx,
            raw_vcf=anno_target,
            transfer=transfer_spec,
            apply=apply_spec,
        ).add_to(dag, {anno_job})

        if self.skip_model_apply:
            return

        # Final normalize
        norm_job = Job(
            cmds.filter_norm(
                ctx.output_vcf,
                apply_vcf,
                ctx.reference,
                exclude_homref=not self.gvcf,
            ),
            "final-norm",
            0,
            task_name="vcf-norm",
        )
        dag.add_job(norm_job, transfer_apply.terminal)
