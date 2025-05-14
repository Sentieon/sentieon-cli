"""
Functionality for the DNAscope hybrid pipeline
"""

import argparse
import copy
import multiprocessing as mp
import json
import pathlib
import shlex
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import packaging.version

from importlib_resources import files

from .archive import ar_load
from . import command_strings as cmds
from .dag import DAG
from .driver import (
    BaseDriver,
    Driver,
    DNAscope,
    DNAModelApply,
    HybridStage1,
    HybridStage2,
    HybridStage3,
)
from .dnascope import call_cnvs, DNAscopePipeline
from .dnascope_longread import DNAscopeLRPipeline
from .job import Job
from .util import (
    __version__,
    check_version,
    library_preloaded,
    parse_rg_line,
    path_arg,
    spit_alignment,
)


CALLING_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202503.01"),
    "bedtools": None,
    "bcftools": packaging.version.Version("1.10"),
    "samtools": packaging.version.Version("1.16"),
}


class DNAscopeHybridPipeline(DNAscopePipeline, DNAscopeLRPipeline):
    """The DNAscope Hybrid pipeline"""

    params: Dict[str, Dict[str, Any]] = {
        "reference": {
            "flags": ["-r", "--reference"],
            "required": True,
            "help": "fasta for reference genome.",
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
        "sr_aln": {
            "nargs": "*",
            "help": "Short-read BAM or CRAM files",
            "type": path_arg(exists=True, is_file=True),
        },
        "sr_duplicate_marking": {
            "help": "Options for duplicate marking.",
            "choices": ["markdup", "rmdup", "none"],
            "default": "markdup",
        },
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
        "dbsnp": {
            "flags": ["-d", "--dbsnp"],
            "help": (
                "dbSNP vcf file Supplying this file will annotate variants "
                "with their dbSNP refSNP ID numbers."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "bed": {
            "flags": ["-b", "--bed"],
            "help": (
                "Region BED file. Supplying this file will limit variant "
                "calling to the intervals inside the BED file."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "cores": {
            "flags": ["-t", "--cores"],
            "help": (
                "Number of threads/processes to use. Defaults to all "
                "available."
            ),
            "default": mp.cpu_count(),
        },
        "gvcf": {
            "flags": ["-g", "--gvcf"],
            "help": (
                "Generate a gVCF output file along with the VCF."
                " (default generates only the VCF)"
            ),
            "action": "store_true",
            "type": bool,
        },
        "dry_run": {
            "help": "Print the commands without running them.",
            "action": "store_true",
            "type": bool,
        },
        "skip_svs": {
            "help": "Skip SV calling",
            "action": "store_true",
            "type": bool,
        },
        "skip_metrics": {
            "help": "Skip all metrics collection and multiQC",
            "action": "store_true",
            "type": bool,
        },
        "skip_mosdepth": {
            "help": "Skip QC with mosdepth.",
            "action": "store_true",
            "type": bool,
        },
        "skip_multiqc": {
            "help": "Skip multiQC report generation.",
            "action": "store_true",
            "type": bool,
        },
        "skip_cnv": {
            "help": "Skip CNV calling.",
            "action": "store_true",
            "type": bool,
        },
        "lr_align_input": {
            "help": (
                "Align the input long-read BAM/CRAM/uBAM file to the "
                "reference genome"
            ),
            "action": "store_true",
            "type": bool,
        },
        "lr_input_ref": {
            "help": (
                "Used to decode the input long-read alignment file. Required "
                "if the input file is in the CRAM/uCRAM formats."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "bam_format": {
            "help": (
                "Use the BAM format instead of CRAM for output aligned files"
            ),
            "action": "store_true",
            "type": bool,
        },
        "rgsm": {
            "help": (
                "Overwrite the SM tag of the input readgroups for "
                "compatibility"
            )
        },
        "lr_fastq_taglist": {
            # help="A comma-separated list of tags to retain. Defaults to "
            # "'%(default)s' and the 'RG' tag is required",
            "help": argparse.SUPPRESS,
            "default": "*",
        },
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
        "minimap2_args": {
            # help="Extra arguments for sentieon minimap2",
            "help": argparse.SUPPRESS,
            "default": "-Y",
        },
        "util_sort_args": {
            # help="Extra arguments for sentieon util sort",
            "help": argparse.SUPPRESS,
            "default": "--cram_write_options version=3.0,compressor=rans",
        },
        "skip_version_check": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
            "type": bool,
        },
        "retain_tmpdir": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
            "type": bool,
        },
        "bwt_max_mem": {
            # Manually set `bwt_max_mem`
            "help": argparse.SUPPRESS,
        },
        "no_split_alignment": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
            "type": bool,
        },
        "sr_read_filter": {
            "help": argparse.SUPPRESS,
        },
        "lr_read_filter": {
            "help": argparse.SUPPRESS,
        },
    }

    positionals: Dict[str, Dict[str, Any]] = {
        "output-vcf": {
            "help": "Output VCF File. The file name must end in .vcf.gz",
            "type": path_arg(),
        },
    }

    def __init__(self) -> None:
        super().__init__()
        self.sr_r1_fastq: List[pathlib.Path] = []
        self.sr_r2_fastq: List[pathlib.Path] = []
        self.sr_readgroups: List[str] = []
        self.sr_aln: List[pathlib.Path] = []
        self.sr_duplicate_marking = "markdup"
        self.lr_aln: List[pathlib.Path] = []
        self.lr_align_input = False
        self.lr_input_ref: Optional[pathlib.Path] = None
        self.bam_format = False
        self.rgsm: Optional[str] = None
        self.lr_fastq_taglist = "*"
        self.sr_read_filter: Optional[str] = None
        self.lr_read_filter: Optional[str] = None
        self.assay = "WGS"

    def validate(self) -> None:
        self.validate_bundle()
        self.collect_readgroups()
        self.validate_readgroups()

        if not str(self.output_vcf).endswith(".vcf.gz"):
            self.logger.error("The output file should end with '.vcf.gz'")
            sys.exit(2)
        if not self.sr_aln or (self.sr_r1_fastq and self.sr_readgroups):
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

    def validate_bundle(self) -> None:
        bundle_info_bytes = ar_load(
            str(self.model_bundle) + "/bundle_info.json"
        )
        if isinstance(bundle_info_bytes, list):
            bundle_info_bytes = b"{}"

        bundle_info = json.loads(bundle_info_bytes.decode())
        self.longread_tech = bundle_info.get("longReadPlatform")
        self.shortread_tech = bundle_info.get("shortReadPlatform")
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
            self.all_readgroups[2][-1].append(
                parse_rg_line(rg.replace(r"\t", "\t"))
            )

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
            self.numa_nodes = spit_alignment(self.cores)
        n_alignment_jobs = max(1, len(self.numa_nodes))

        self.set_bwt_max_mem(0, n_alignment_jobs)

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

        assert self.output_vcf
        assert self.reference
        assert self.model_bundle

        # Short-read alignment
        sr_aln = self.sr_aln
        (
            aligned_fastq,
            align_fastq_jobs,
            fq_rm_job,
        ) = self.sr_align_fastq()
        for job in align_fastq_jobs:
            dag.add_job(job)
        sr_aln.extend(aligned_fastq)

        # Short-read dedup
        (
            deduped,
            lc_job,
            dedup_job,
            metrics_job,
            rehead_job,
        ) = self.dedup_and_metrics(
            sample_input=sr_aln,
        )
        sr_aln = deduped
        if lc_job:
            dag.add_job(lc_job, align_fastq_jobs)
            if dedup_job:
                dag.add_job(dedup_job, {lc_job})
                if metrics_job and not self.skip_metrics:
                    dag.add_job(metrics_job, {dedup_job})
                    if rehead_job:
                        dag.add_job(rehead_job, {metrics_job})
                if fq_rm_job:
                    dag.add_job(fq_rm_job, {dedup_job})

        # Long-read alignment
        realign_jobs: Set[Job] = set()
        lr_aln = self.lr_aln
        if self.lr_align_input:
            lr_aln, realign_jobs = self.lr_align_inputs()
            for job in realign_jobs:
                dag.add_job(job)

        if not self.skip_mosdepth:
            mosdpeth_jobs = self.mosdepth(sample_input=lr_aln)
            for job in mosdpeth_jobs:
                dag.add_job(job, realign_jobs)

        if not self.skip_multiqc:
            multiqc_job = self.multiqc()
            multiqc_dependencies: Set[Job] = set()
            if lc_job:
                multiqc_dependencies.add(lc_job)
            if metrics_job:
                multiqc_dependencies.add(metrics_job)
            if rehead_job:
                multiqc_dependencies.add(rehead_job)

            if multiqc_job:
                dag.add_job(multiqc_job, multiqc_dependencies)

        if not self.skip_svs:
            longreadsv_job = self.call_svs(lr_aln)
            dag.add_job(longreadsv_job, realign_jobs)

        sr_preprocessing_jobs: Set[Job] = set()
        sr_preprocessing_jobs.update(align_fastq_jobs)
        if dedup_job:
            sr_preprocessing_jobs.add(dedup_job)
        if not self.skip_cnv:
            cnvscope_job, cnvmodelapply_job = call_cnvs(
                self.tmp_dir,
                self.output_vcf,
                self.reference,
                sr_aln,
                self.model_bundle,
                self.bed,
                self.cores,
                self.skip_version_check,
            )
            dag.add_job(cnvscope_job, sr_preprocessing_jobs)
            dag.add_job(cnvmodelapply_job, {cnvscope_job})

        (
            call_job,
            select_job,
            mapq0_job,
            mapq0_slop_job,
            cat_merge_job,
            rm_job1,
            stage1_job,
            rm_job2,
            second_stage_job,
            rm_job3,
            third_stage_job,
            rm_job4,
            call2_job,
            subset_job,
            concat_job,
            rm_job5,
            anno_job,
            apply_job,
            norm_job,
        ) = self.call_variants(sr_aln, lr_aln)
        dag.add_job(call_job, realign_jobs | sr_preprocessing_jobs)
        dag.add_job(select_job, {call_job})
        dag.add_job(mapq0_job, realign_jobs | sr_preprocessing_jobs)
        dag.add_job(mapq0_slop_job, {mapq0_job})
        dag.add_job(cat_merge_job, {mapq0_slop_job, select_job})
        dag.add_job(stage1_job, {cat_merge_job})
        dag.add_job(second_stage_job, {stage1_job})
        dag.add_job(third_stage_job, {second_stage_job})
        dag.add_job(call2_job, {third_stage_job})
        dag.add_job(subset_job, {second_stage_job})
        dag.add_job(concat_job, {subset_job, call2_job})
        dag.add_job(anno_job, {concat_job})
        dag.add_job(apply_job, {anno_job})
        dag.add_job(norm_job, {apply_job})

        # Remove intermediate files during processing
        if not self.retain_tmpdir:
            dag.add_job(rm_job1, {cat_merge_job})
            dag.add_job(rm_job2, {stage1_job})
            dag.add_job(rm_job3, {second_stage_job})
            dag.add_job(rm_job4, {third_stage_job})
            dag.add_job(rm_job5, {concat_job})

        return dag

    def call_variants(
        self,
        sr_aln: List[pathlib.Path],
        lr_aln: List[pathlib.Path],
        **_kwargs: Any,
    ) -> Tuple[
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
        Job,
    ]:
        """
        Call SNVs and indels using the DNAscope hybrid pipeline
        """
        assert self.output_vcf
        assert self.reference
        assert self.model_bundle

        ref_fai = pathlib.Path(str(self.reference) + ".fai")

        if not self.skip_version_check:
            for cmd, min_version in CALLING_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        # readgroup handling for long-reads
        tr_read_filter: List[str] = []
        lr_rg_read_filter: List[str] = []
        replace_rg_args: Tuple[List[List[str]], List[List[str]]] = ([], [])
        for aln_rgs in self.lr_aln_readgroups:
            replace_rg_args[0].append([])
            for rg_line_d in aln_rgs:
                id = rg_line_d.get("ID")
                new_sm = (
                    self.hybrid_rg_sm
                    if self.hybrid_set_rg
                    else rg_line_d.get("SM")
                )
                replace_rg_args[0][-1].append(
                    f"{id}=ID:{id}\\tSM:{new_sm}\\tLR:1"
                )
                if self.lr_read_filter:
                    lr_rg_read_filter.append(
                        f"{self.lr_read_filter},rgid={id}"
                    )
                if self.longread_tech.upper() == "ONT":
                    tr_read_filter.append(
                        f"TrimRepeat,max_repeat_unit=2,min_repeat_span=6,"
                        f"rgid={id}"
                    )

        # readgroup handling for short-reads
        sr_rg_read_filter: List[str] = []
        ultima_read_filter: List[str] = []
        for aln_rgs in self.sr_aln_readgroups:
            if self.hybrid_set_rg:
                replace_rg_args[1].append([])
            for rg_line_d in aln_rgs:
                id = rg_line_d.get("ID")
                new_sm = (
                    self.hybrid_rg_sm
                    if self.hybrid_set_rg
                    else rg_line_d.get("SM")
                )
                if self.hybrid_set_rg:
                    replace_rg_args[1][-1].append(
                        f"{id}=ID:{id}\\tSM:{new_sm}"
                    )
                if self.sr_read_filter:
                    sr_rg_read_filter.append(
                        f"{self.sr_read_filter},rgid={id}"
                    )
                if self.shortread_tech.upper() == "ULTIMA":
                    ultima_read_filter.append(f"UltimaReadFilter,rgid={id}")

        # First pass - combined variant calling
        vcf_suffix = ".g.vcf.gz" if self.gvcf else ".vcf.gz"
        combined_vcf = self.tmp_dir.joinpath("initial" + vcf_suffix)
        emit_mode = "gvcf" if self.gvcf else None
        driver: BaseDriver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            replace_rg=replace_rg_args[0] + replace_rg_args[1],
            input=lr_aln + sr_aln,
            interval=self.bed,
            read_filter=tr_read_filter
            + ultima_read_filter
            + lr_rg_read_filter
            + sr_rg_read_filter,
        )
        driver.add_algo(
            DNAscope(
                combined_vcf,
                dbsnp=self.dbsnp,
                model=self.model_bundle.joinpath("hybrid.model"),
                pcr_indel_model="none",
                emit_mode=emit_mode,
            )
        )
        call_job = Job(shlex.join(driver.build_cmd()), "calling-1", self.cores)

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
        )

        mapq0_bed = self.tmp_dir.joinpath("hybrid_mapq0.bed")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            replace_rg=replace_rg_args[0] + replace_rg_args[1],
            input=lr_aln + sr_aln,
            read_filter=lr_rg_read_filter + sr_rg_read_filter,
        )
        driver.add_algo(
            HybridStage2(
                model=self.model_bundle.joinpath("HybridStage2_region.model"),
                all_bed=mapq0_bed,
            )
        )
        mapq0_job = Job(
            shlex.join(driver.build_cmd()),
            "mapq0-bed",
            self.cores,
        )

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
        )

        diff_bed = self.tmp_dir.joinpath("merged_diff.bed")
        cat_merge_job = Job(
            cmds.cmd_bedtools_cat_sort_merge(
                out_bed=diff_bed,
                in_bed=[selected_bed, mapq0_slop_bed],
                ref_fai=ref_fai,
            ),
            "concat-merge-bed",
            0,
        )
        rm_cmd = ["rm", str(selected_bed), str(mapq0_slop_bed)]
        rm_job1 = Job(shlex.join(rm_cmd), "rm-tmp1", 0, True)

        stage1_ins_fa = self.tmp_dir.joinpath("stage1_ins.fa")
        stage1_ins_bed = self.tmp_dir.joinpath("stage1_ins.bed")
        ins_driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            replace_rg=replace_rg_args[0],
            input=lr_aln,
            read_filter=lr_rg_read_filter,
        )
        ins_driver.add_algo(
            HybridStage1(
                "-",
                model=self.model_bundle.joinpath("HybridStage1_ins.model"),
                fa_file=stage1_ins_fa,
                bed_file=stage1_ins_bed,
            )
        )

        stage1_hap_bam = self.tmp_dir.joinpath("stage1_hap.bam")
        stage1_hap_bed = self.tmp_dir.joinpath("stage1_hap.bed")
        stage1_hap_vcf = self.tmp_dir.joinpath("stage1_hap.vcf")
        stage1_driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            replace_rg=replace_rg_args[0],
            input=lr_aln,
            interval=diff_bed,
            read_filter=lr_rg_read_filter,
        )
        stage1_driver.add_algo(
            HybridStage1(
                "-",
                model=self.model_bundle.joinpath("HybridStage1.model"),
                hap_bam=stage1_hap_bam,
                hap_bed=stage1_hap_bed,
                hap_vcf=stage1_hap_vcf,
            )
        )

        stage1_bam = self.tmp_dir.joinpath("hybrid_stage1.bam")
        stage1_job = Job(
            cmds.hybrid_stage1(
                out_aln=stage1_bam,
                reference=self.reference,
                cores=self.cores,
                readgroup=f"@RG\\tID:hybrid-18893\\tSM:{self.hybrid_rg_sm}",
                ins_driver=ins_driver,
                stage1_driver=stage1_driver,
                bwa_model=self.model_bundle.joinpath("HybridStage1_bwa.model"),
            ),
            "first-stage",
            self.cores,
        )
        rm_cmd = [
            "rm",
            str(stage1_ins_fa),
            str(stage1_ins_bed),
            str(stage1_hap_vcf),
        ]
        rm_job2 = Job(shlex.join(rm_cmd), "rm-tmp2", 0, True)

        stage2_bed = self.tmp_dir.joinpath("hybrid_stage2.bed")
        stage2_unmap_bam = self.tmp_dir.joinpath("hybrid_stage2_unmap.bam")
        stage2_alt_bam = self.tmp_dir.joinpath("hybrid_stage2_alt.bam")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[stage1_bam, stage1_hap_bam],
        )
        driver.add_algo(
            HybridStage2(
                model=self.model_bundle.joinpath("HybridStage2.model"),
                hap_bed=stage1_hap_bed,
                unmap_bam=stage2_unmap_bam,
                alt_bam=stage2_alt_bam,
                all_bed=stage2_bed,
            )
        )
        second_stage_job = Job(
            shlex.join(driver.build_cmd()),
            "second-stage",
            self.cores,
        )

        rm_cmd = ["rm", str(stage1_bam), str(stage1_hap_bam)]
        rm_job3 = Job(shlex.join(rm_cmd), "rm-tmp3", 0, True)

        stage3_bam = self.tmp_dir.joinpath("hybrid_stage3.bam")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            replace_rg=replace_rg_args[0] + replace_rg_args[1],
            input=lr_aln + sr_aln + [stage2_unmap_bam, stage2_alt_bam],
            interval=stage2_bed,
            read_filter=lr_rg_read_filter + sr_rg_read_filter,
        )
        driver.add_algo(
            HybridStage3(
                "-",
                model=self.model_bundle.joinpath("HybridStage3.model"),
            )
        )
        third_stage_job = Job(
            cmds.hybrid_stage3(
                stage3_bam,
                driver=driver,
                cores=self.cores,
            ),
            "third-stage",
            self.cores,
        )
        rm_cmd = ["rm", str(stage2_unmap_bam), str(stage2_alt_bam)]
        rm_job4 = Job(shlex.join(rm_cmd), "rm-tmp4", 0, True)

        # pass 2 of variant calling
        pass2_vcf = self.tmp_dir.joinpath("hybrid_pass2.vcf.gz")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            replace_rg=replace_rg_args[0],
            input=lr_aln + [stage3_bam],
            interval=stage2_bed,
            read_filter=tr_read_filter
            + ultima_read_filter
            + lr_rg_read_filter,
        )
        driver.add_algo(
            DNAscope(
                pass2_vcf,
                dbsnp=self.dbsnp,
                model=self.model_bundle.joinpath("hybrid.model"),
                pcr_indel_model="none",
                emit_mode=emit_mode,
            )
        )
        call2_job = Job(shlex.join(driver.build_cmd()), "call2", self.cores)

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
        )
        concat_job = Job(
            cmds.bcftools_concat(
                combined_tmp_vcf,
                [subset_vcf, pass2_vcf],
            ),
            "concat-calls",
            0,
        )
        rm_cmd = ["rm", str(combined_vcf), str(subset_vcf), str(pass2_vcf)]
        rm_job5 = Job(shlex.join(rm_cmd), "rm-tmp5", 0, True)

        # Annotate the output VCF
        hybrid_anno = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("hybrid_anno.py"))
        )
        combined_anno_vcf = self.tmp_dir.joinpath("combined_tmp_anno.vcf.gz")
        anno_job = Job(
            cmds.cmd_pyexec_hybrid_anno(
                combined_anno_vcf,
                combined_tmp_vcf,
                stage1_hap_bed,
                hybrid_anno,
                self.cores,
            ),
            "anno-calls",
            0,
        )

        # Model Apply
        apply_vcf = self.tmp_dir.joinpath("combined_apply.vcf.gz")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            interval=self.bed,
        )
        driver.add_algo(
            DNAModelApply(
                model=self.model_bundle.joinpath("hybrid.model"),
                vcf=combined_anno_vcf,
                output=apply_vcf,
            )
        )
        apply_job = Job(
            shlex.join(driver.build_cmd()), "model-apply", self.cores
        )

        # Final normalize
        norm_job = Job(
            cmds.filter_norm(
                self.output_vcf,
                apply_vcf,
                self.reference,
                exclude_homref=not self.gvcf,
            ),
            "final-norm",
            0,
        )
        return (
            call_job,
            select_job,
            mapq0_job,
            mapq0_slop_job,
            cat_merge_job,
            rm_job1,
            stage1_job,
            rm_job2,
            second_stage_job,
            rm_job3,
            third_stage_job,
            rm_job4,
            call2_job,
            subset_job,
            concat_job,
            rm_job5,
            anno_job,
            apply_job,
            norm_job,
        )
