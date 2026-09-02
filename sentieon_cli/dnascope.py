"""
DNAscope alignment and variant calling
"""

import argparse
import copy
from dataclasses import dataclass, field
import itertools
import os
import pathlib
import shutil
import sys
from typing import List, Optional, Set, Tuple

import packaging.version

from . import command_strings as cmds
from .dag import DAG
from .driver import (
    AlignmentStat,
    BaseAlgo,
    BaseDistributionByCycle,
    CoverageMetrics,
    DNAscope,
    GCBias,
    HsMetricAlgo,
    InsertSizeMetricAlgo,
    MeanQualityByCycle,
    QualDistribution,
    SVSolver,
    WgsMetricsAlgo,
)
from .job import Job
from .pipeline import BasePipeline
from .shell_pipeline import Command, Pipeline
from .stages.base import StageContext, driver_job, rm_job
from .stages.dedup import DedupStage
from .stages.metrics import MetricsPaths, MetricsStage
from .stages.small_variants import (
    ApplySpec,
    DNAscopeStage,
    GVCFtyperStage,
    TransferApplyStage,
)
from .util import (
    check_version,
    library_preloaded,
    parse_rg_line,
    path_arg,
    split_alignment,
    total_memory,
)

ALN_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
    "samtools": packaging.version.Version("1.16"),
}

FQ_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}

VARIANTS_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}


@dataclass
class SrPreprocessing:
    """The short-read dedup and metrics jobs added to a DAG.

    * ``deduped`` -- the alignment to call variants from; the input files
      unchanged when duplicate marking is off.
    * ``dedup_job`` -- the Dedup job, absent when duplicate marking is off.
    * ``qc_jobs`` -- every metrics-producing job, for MultiQC to wait on.
    """

    deduped: List[pathlib.Path]
    dedup_job: Optional[Job] = None
    qc_jobs: List[Job] = field(default_factory=list)


class DNAscopePipeline(BasePipeline):
    """The DNAscope pipeline"""

    params = copy.deepcopy(BasePipeline.params)
    params.update(
        {
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
            "readgroups": {
                "nargs": "*",
                "help": "Readgroup information for the fastq files.",
            },
            "sample_input": {
                "flags": ["-i", "--sample_input"],
                "nargs": "*",
                "help": "sample BAM or CRAM file.",
                "type": path_arg(exists=True, is_file=True),
            },
            "align": {
                "help": (
                    "Align the reads in the input uBAM/uCRAM file to the "
                    "reference genome. Assumes paired reads are collated in "
                    "the input file."
                ),
                "action": "store_true",
            },
            "assay": {
                "help": "The type of assay, WGS or WES.",
                "choices": ["WGS", "WES"],
                "default": "WGS",
            },
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
            "collate_align": {
                "help": (
                    "Collate and align the reads in the input BAM/CRAM file "
                    "to the reference genome. Suitable for coordinate-sorted "
                    "BAM/CRAM input."
                ),
                "action": "store_true",
            },
            "consensus": {
                "help": "Generate consensus reads during dedup",
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
            "duplicate_marking": {
                "help": "Options for duplicate marking.",
                "choices": ["markdup", "rmdup", "none"],
                "default": "markdup",
            },
            "gvcf": {
                "flags": ["-g", "--gvcf"],
                "help": (
                    "Generate a gVCF output file along with the VCF."
                    " (default generates only the VCF)"
                ),
                "action": "store_true",
            },
            "input_ref": {
                "help": (
                    "Used to decode the input alignment file. Required if the "
                    "input file is in the CRAM/uCRAM formats."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "interval_padding": {
                "help": "Amount to pad all intervals.",
                "type": int,
            },
            "pcr_free": {
                "help": "Use arguments for PCR-free data processing",
                "action": "store_true",
            },
            "skip_metrics": {
                "help": "Skip all metrics collection and multiQC",
                "action": "store_true",
            },
            "skip_multiqc": {
                "help": "Skip multiQC report generation",
                "action": "store_true",
            },
            "skip_small_variants": {
                "help": "Skip small variant (SNV/indel) calling",
                "action": "store_true",
            },
            "skip_svs": {
                "help": "Skip SV calling",
                "action": "store_true",
            },
            "bwa_args": {
                # "help": "Extra arguments for sentieon bwa",
                "help": argparse.SUPPRESS,
                "default": "",
            },
            "bwa_k": {
                # "help": "The '-K' argument in bwa",
                "help": argparse.SUPPRESS,
                "default": 100000000,
            },
            "bwt_max_mem": {
                # Manually set `bwt_max_mem`
                "help": argparse.SUPPRESS,
            },
            "no_ramdisk": {
                # Do not use /dev/shm, even on high memory machines
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "no_split_alignment": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "util_sort_args": {
                # "help": "Extra arguments for sentieon util sort",
                "help": argparse.SUPPRESS,
                "default": "--cram_write_options version=3.0,compressor=rans",
            },
        }
    )

    def __init__(self) -> None:
        super().__init__()
        self.sample_input: List[pathlib.Path] = []
        self.r1_fastq: List[pathlib.Path] = []
        self.r2_fastq: List[pathlib.Path] = []
        self.readgroups: List[str] = []
        self.model_bundle: Optional[pathlib.Path] = None
        self.dbsnp: Optional[pathlib.Path] = None
        self.bed: Optional[pathlib.Path] = None
        self.interval_padding: Optional[int] = None
        self.pcr_free = False
        self.gvcf = False
        self.duplicate_marking = "markdup"
        self.assay = "WGS"
        self.consensus = False
        self.skip_small_variants = False
        self.skip_svs = False
        self.skip_metrics = False
        self.skip_multiqc = False
        self.align = False
        self.collate_align = False
        self.input_ref: Optional[pathlib.Path] = None
        self.bam_format = False
        self.bwa_args = ""
        self.bwa_k = 100000000
        self.util_sort_args = (
            "--cram_write_options version=3.0,compressor=rans"
        )
        self.bwt_max_mem: Optional[str] = None
        self.no_ramdisk = False
        self.no_split_alignment = False

    def validate(self) -> None:
        if not self.output_vcf:
            self.logger.error("output_vcf is required")
            sys.exit(2)

        # uniquify pipeline attributes
        self.sr_r1_fastq = self.r1_fastq
        self.sr_r2_fastq = self.r2_fastq
        self.sr_readgroups = self.readgroups
        self.sr_duplicate_marking = self.duplicate_marking
        del self.r1_fastq
        del self.r2_fastq
        del self.readgroups
        del self.duplicate_marking

        # validate
        if not self.sample_input and not (
            self.sr_r1_fastq and self.sr_readgroups
        ):
            self.logger.error(
                "Please supply either the `--sample_input` or `--r1_fastq` "
                "and `--readgroups` arguments"
            )
            sys.exit(2)
        self.validate_output_vcf()
        self.skip_multiqc = True if self.skip_metrics else self.skip_multiqc

        if not library_preloaded("libjemalloc.so"):
            self.logger.warning(
                "jemalloc is recommended, but is not preloaded. See "
                "https://support.sentieon.com/appnotes/jemalloc/"
            )

        if self.bed is None:
            if self.assay == "WES":
                self.logger.warning(
                    "A BED file is recommended with WES assays to restrict "
                    "variant calling to target regions."
                )
            else:
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

        for rg in self.sr_readgroups:
            try:
                parse_rg_line(rg.replace(r"\t", "\t"))
            except ValueError as e:
                self.logger.error("Invalid --readgroups value '%s': %s", rg, e)
                sys.exit(2)

        if self.sr_r1_fastq or self.align or self.collate_align:
            self.validate_bwa_index()

    def configure(self) -> None:
        self.configure_alignment()

    def configure_alignment(self) -> None:
        self.numa_nodes: List[str] = []
        n_alignment_jobs = 1
        if not self.no_split_alignment:
            self.numa_nodes = split_alignment(self.cores)
        n_alignment_jobs = max(1, len(self.numa_nodes))

        total_input_size = 0
        shm_ok = False
        if not self.no_ramdisk:
            total_input_size = self.total_input_size()
            shm_ok = self.check_shm(total_input_size, n_alignment_jobs)
            if shm_ok:
                os.environ["SENTIEON_TMPDIR"] = "/dev/shm"

        self.set_bwt_max_mem(
            total_input_size if shm_ok else 0, n_alignment_jobs
        )

    def total_input_size(self) -> int:
        """Find the total size of all inputs"""
        total_input_size = sum([x.stat().st_size for x in self.sample_input])
        for r1, r2 in itertools.zip_longest(
            self.sr_r1_fastq, self.sr_r2_fastq
        ):
            for fq in (r1, r2):
                if isinstance(fq, pathlib.Path):
                    total_input_size += fq.stat().st_size
        return total_input_size

    def check_shm(self, total_input_size: int, n_alignment_jobs) -> bool:
        """Check if /dev/shm can be used for temporary files"""
        try:
            shm_free = shutil.disk_usage("/dev/shm").free
        except Exception:
            return False
        total_mem = total_memory()

        if (
            shm_free > total_input_size * 2.3
            and total_mem > total_input_size + 60 * 1024**3 * n_alignment_jobs
        ):
            self.logger.debug("Using /dev/shm for temporary files")
            return True
        return False

    def set_bwt_max_mem(
        self,
        total_input_size: int,
        n_alignment_jobs: int = 1,
    ):
        """Set the bwt_max_mem environment variable"""
        if self.bwt_max_mem:
            os.environ["bwt_max_mem"] = self.bwt_max_mem
            return

        total_mem = total_memory()
        total_mem_gb = total_mem / (1024.0**3)
        align_mem_gb = (
            total_mem_gb - 4 - total_input_size / (1024.0**3) * 2.3
        )  # some memory for other system processes
        bwa_mem_gb = max(
            int((align_mem_gb / n_alignment_jobs) - 6), 0
        )  # some memory for other alignment processes
        self.logger.debug("Setting bwt_max_mem to: %sG", bwa_mem_gb)
        os.environ["bwt_max_mem"] = f"{bwa_mem_gb}G"

    def build_dag(self) -> DAG:
        """Build the DAG for the pipeline"""
        self.logger.info("Building the DAG")
        dag = DAG()
        ctx = self.stage_context()

        # Alignment
        align_jobs: Set[Job] = set()
        sample_input = copy.deepcopy(self.sample_input)
        bam_rm_job = None
        if self.align or self.collate_align:
            sample_input, align_jobs, bam_rm_job = self.sr_align_inputs()
            for job in align_jobs:
                dag.add_job(job)
        aligned_fastq, align_fastq_jobs, fq_rm_job = self.sr_align_fastq()
        for job in align_fastq_jobs:
            dag.add_job(job)
        sample_input += aligned_fastq

        # Dedup and metrics
        preprocessing = self.add_sr_preprocessing(
            dag, ctx, sample_input, align_jobs.union(align_fastq_jobs)
        )
        deduped = preprocessing.deduped
        dedup_job = preprocessing.dedup_job
        if dedup_job:
            if bam_rm_job:
                dag.add_job(bam_rm_job, {dedup_job})
            if fq_rm_job:
                dag.add_job(fq_rm_job, {dedup_job})

        # Small variants
        if not self.skip_small_variants:
            call_dependencies: Set[Job] = set()
            if dedup_job:
                call_dependencies.add(dedup_job)
            else:
                call_dependencies.update(align_jobs)
                call_dependencies.update(align_fastq_jobs)
            self.add_small_variant_calling(
                dag, ctx, deduped, call_dependencies
            )

        # Multiqc
        if not self.skip_multiqc:
            multiqc_job = self.multiqc()
            if multiqc_job:
                dag.add_job(multiqc_job, set(preprocessing.qc_jobs))
        return dag

    def sr_align_inputs(self) -> Tuple[List[pathlib.Path], Set[Job], Job]:
        """Align input BAM/CRAM/uBAM/uCRAM files with bwa"""
        if not self.reference:
            self.logger.error("reference is required")
            sys.exit(2)
        if not self.model_bundle:
            self.logger.error("model_bundle is required")
            sys.exit(2)

        if not self.skip_version_check:
            for cmd, min_version in ALN_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        res: List[pathlib.Path] = []
        suffix = "bam" if self.bam_format else "cram"
        util_sort_args = self.util_sort_args
        if self.sr_duplicate_marking != "none":
            suffix = "bam"
            util_sort_args += " --bam_compression 1 "
        jobs = set()
        align_outputs: List[pathlib.Path] = []
        for i, input_aln in enumerate(self.sample_input):
            out_aln = pathlib.Path(
                str(self.output_vcf).replace(
                    ".vcf.gz", f"_bwa_sorted_{i}.{suffix}"
                )
            )
            if self.sr_duplicate_marking != "none":
                out_aln = self.tmp_dir.joinpath(f"bwa_sorted_{i}.{suffix}")
            rg_lines = cmds.get_rg_lines(input_aln, self.dry_run)
            rg_header = self.tmp_dir.joinpath(f"input_{i}.hdr")
            with open(rg_header, "w", encoding="utf-8") as rg_fh:
                for line in rg_lines:
                    print(line, file=rg_fh)
            job = Job(
                cmds.cmd_samtools_fastq_bwa(
                    out_aln,
                    input_aln,
                    self.reference,
                    self.model_bundle,
                    self.cores,
                    rg_header,
                    self.input_ref,
                    collate=self.collate_align,
                    bwa_args=self.bwa_args,
                    bwa_k=str(self.bwa_k),
                    util_sort_args=util_sort_args,
                ),
                f"bam-align-{i}",
                self.cores,
                task_name="alignment",
            )
            res.append(out_aln)
            jobs.add(job)
            align_outputs.append(out_aln)
            align_outputs.append(pathlib.Path(str(out_aln) + ".bai"))
            if suffix == "cram":
                align_outputs.append(pathlib.Path(str(out_aln) + ".crai"))

        # Create an unscheduled job to remove the aligned inputs
        rm_job = Job(
            Pipeline(
                Command("rm", *[str(x) for x in align_outputs], fail_ok=True)
            ),
            "rm-bam-aln",
            0,
            task_name="cleanup",
        )

        return (res, jobs, rm_job)

    def sr_align_fastq(
        self,
    ) -> Tuple[List[pathlib.Path], Set[Job], Optional[Job]]:
        """Align fastq files to the reference genome using bwa"""
        if not self.reference:
            self.logger.error("reference is required")
            sys.exit(2)
        if not self.model_bundle:
            self.logger.error("model_bundle is required")
            sys.exit(2)

        res: List[pathlib.Path] = []
        jobs: Set[Job] = set()
        if not self.sr_r1_fastq and not self.sr_readgroups:
            return (res, jobs, None)

        if not self.skip_version_check:
            for cmd, min_version in FQ_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        unzip = "igzip"
        if not shutil.which(unzip):
            self.logger.warning(
                "igzip is recommended for decompression, but is not "
                "available. Falling back to gzip."
            )
            unzip = "gzip"

        n_alignment_jobs = max(1, len(self.numa_nodes))
        suffix = "bam" if self.bam_format else "cram"
        util_sort_args = self.util_sort_args
        if self.sr_duplicate_marking != "none":
            suffix = "bam"
            util_sort_args += " --bam_compression 1 "
        align_outputs: List[pathlib.Path] = []
        for i, (r1, r2, rg) in enumerate(
            itertools.zip_longest(
                self.sr_r1_fastq, self.sr_r2_fastq, self.sr_readgroups
            )
        ):
            for j in range(n_alignment_jobs):
                out_aln = pathlib.Path(
                    str(self.output_vcf).replace(
                        ".vcf.gz", f"_bwa_sorted_fq_{i}_{j}.{suffix}"
                    )
                )
                if self.sr_duplicate_marking != "none":
                    out_aln = self.tmp_dir.joinpath(
                        f"bwa_sorted_fq_{i}_{j}.{suffix}"
                    )
                numa = self.numa_nodes[j] if len(self.numa_nodes) > 1 else None
                split = (
                    f"{j}/{n_alignment_jobs}"
                    if len(self.numa_nodes) > 1
                    else None
                )
                split_cores = max(1, int(self.cores / n_alignment_jobs))
                job = Job(
                    cmds.cmd_fastq_bwa(
                        out_aln,
                        r1,
                        r2,
                        rg,
                        self.reference,
                        self.model_bundle,
                        split_cores,
                        unzip,
                        self.bwa_args,
                        str(self.bwa_k),
                        util_sort_args,
                        numa,
                        split,
                    ),
                    f"bam-align-{i}-{j}",
                    split_cores,
                    resources={f"node{j}": 1},
                    task_name="alignment",
                )
                res.append(out_aln)
                jobs.add(job)
                align_outputs.append(out_aln)
                align_outputs.append(pathlib.Path(str(out_aln) + ".bai"))
                if suffix == "cram":
                    align_outputs.append(pathlib.Path(str(out_aln) + ".crai"))

        # Create an unscheduled job to remove the aligned inputs
        rm_job = Job(
            Pipeline(
                Command("rm", *[str(x) for x in align_outputs], fail_ok=True)
            ),
            "rm-fq-aln",
            0,
            task_name="cleanup",
        )

        return (res, jobs, rm_job)

    def sr_metrics_algos(self, paths: MetricsPaths) -> List[BaseAlgo]:
        """The metrics that do not need duplicate-marked reads"""
        algos: List[BaseAlgo] = [
            MeanQualityByCycle(paths.mean_qual_by_cycle),
            BaseDistributionByCycle(paths.base_distribution_by_cycle),
            QualDistribution(paths.qual_distribution),
            AlignmentStat(paths.alignment_stat),
        ]
        if self.assay == "WGS":
            algos.append(GCBias(paths.gc_bias, summary=paths.gc_bias_summary))
        return algos

    def add_sr_preprocessing(
        self,
        dag: DAG,
        ctx: StageContext,
        sample_input: List[pathlib.Path],
        upstream: Set[Job],
    ) -> SrPreprocessing:
        """Mark duplicates and collect metrics from the short reads"""
        if not self.output_vcf:
            self.logger.error("output_vcf is required")
            sys.exit(2)
        suffix = "bam" if self.bam_format else "cram"

        paths = MetricsPaths.from_output_vcf(self.output_vcf)
        paths.ensure_dir(self.dry_run)

        if self.sr_duplicate_marking == "none":
            # Without duplicate marking there is nothing to collect after
            # the fact, so every metric comes from one pass over the input
            if self.skip_metrics:
                return SrPreprocessing(
                    deduped=sample_input, dedup_job=None, qc_jobs=[]
                )
            metrics_result = MetricsStage(
                ctx=ctx,
                inputs=sample_input,
                algos=[
                    InsertSizeMetricAlgo(paths.insert_size),
                    *self.sr_metrics_algos(paths),
                ],
                name="metrics",
                task_name="metrics",
                threads=ctx.cores,
            ).add_to(dag, upstream)
            return SrPreprocessing(
                deduped=sample_input,
                dedup_job=None,
                qc_jobs=list(metrics_result.jobs),
            )

        # Metrics that ride along on the LocusCollector pass
        lc_extra_algos: List[BaseAlgo] = []
        if not self.skip_metrics:
            # Prefer to run InsertSizeMetricAlgo after duplicate marking
            if self.assay == "WES" and not self.bed:
                lc_extra_algos.append(InsertSizeMetricAlgo(paths.insert_size))
            lc_extra_algos.extend(self.sr_metrics_algos(paths))

        deduped = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", f"_deduped.{suffix}")
        )
        dedup_result = DedupStage(
            ctx=ctx,
            inputs=sample_input,
            output=deduped,
            score_file=paths.score,
            consensus=self.consensus,
            rmdup=(self.sr_duplicate_marking == "rmdup"),
            cram_write_options="version=3.0,compressor=rans",
            dedup_metrics=paths.dedup_metrics,
            lc_extra_algos=lc_extra_algos,
        ).add_to(dag, upstream)
        qc_jobs: List[Job] = [dedup_result.lc_job]

        # HsMetricAlgo and WgsMetricsAlgo run after duplicate marking to
        # account for duplicate reads
        post_algos: List[BaseAlgo] = []
        rehead_metrics: Optional[pathlib.Path] = None
        if not self.skip_metrics:
            if self.assay == "WES" and self.bed:
                post_algos = [
                    HsMetricAlgo(paths.hybrid_selection, self.bed, self.bed),
                    InsertSizeMetricAlgo(paths.insert_size),
                ]
            elif self.assay == "WGS":
                post_algos = [
                    InsertSizeMetricAlgo(paths.insert_size),
                    WgsMetricsAlgo(paths.wgs, include_unpaired="true"),
                    CoverageMetrics(paths.coverage),
                ]
                # Rehead WGS metrics so they are recognized by MultiQC
                rehead_metrics = paths.wgs
        if post_algos:
            metrics_result = MetricsStage(
                ctx=ctx,
                inputs=[deduped],
                algos=post_algos,
                interval=self.bed,
                rehead_metrics=rehead_metrics,
                threads=0,  # Run metrics in the background
            ).add_to(dag, {dedup_result.dedup_job})
            qc_jobs.extend(metrics_result.jobs)

        return SrPreprocessing(
            deduped=[deduped],
            dedup_job=dedup_result.dedup_job,
            qc_jobs=qc_jobs,
        )

    def add_small_variant_calling(
        self,
        dag: DAG,
        ctx: StageContext,
        deduped: List[pathlib.Path],
        upstream: Set[Job],
    ) -> None:
        """Call SNVs, indels, and SVs using DNAscope"""
        if not self.model_bundle:
            self.logger.error("model_bundle is required")
            sys.exit(2)
        if not self.output_vcf:
            self.logger.error("output_vcf is required")
            sys.exit(2)

        if not self.skip_version_check:
            for cmd, min_version in VARIANTS_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        out_gvcf = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", ".g.vcf.gz")
        )
        out_svs_tmp = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_svs_tmp.vcf.gz")
        )
        out_svs = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_svs.vcf.gz")
        )
        pcr_indel_model = "NONE" if self.pcr_free else "CONSERVATIVE"
        model = self.model_bundle.joinpath("dnascope.model")

        # Call variants with DNAscope
        emit_mode = "variant"
        ds_out = self.output_vcf
        tmp_vcf = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_tmp.vcf.gz")
        )
        if self.gvcf:
            emit_mode = "gvcf"
            ds_out = out_gvcf
            tmp_vcf = pathlib.Path(
                str(self.output_vcf).replace(".vcf.gz", "_tmp.g.vcf.gz")
            )

        algos: List[BaseAlgo] = [
            DNAscope(
                tmp_vcf,
                dbsnp=self.dbsnp,
                emit_mode=emit_mode,
                pcr_indel_model=pcr_indel_model,
                model=model,
            )
        ]
        if not self.skip_svs:
            algos.append(
                DNAscope(
                    out_svs_tmp,
                    dbsnp=self.dbsnp,
                    var_type="BND",
                )
            )
        call = DNAscopeStage(
            ctx=ctx,
            algos=algos,
            inputs=deduped,
            interval=self.bed,
            interval_padding=self.interval_padding,
        ).add_to(dag, upstream)

        # Genotyping and filtering with DNAModelApply
        transfer_apply = TransferApplyStage(
            ctx=ctx,
            raw_vcf=tmp_vcf,
            apply=ApplySpec(model=model, output=ds_out),
        ).add_to(dag, call.terminal)

        # Remove the tmp_vcf
        dag.add_job(
            rm_job([tmp_vcf, str(tmp_vcf) + ".tbi"], "rm-tmp-vcf"),
            transfer_apply.terminal,
        )

        # Genotype gVCFs
        if self.gvcf:
            GVCFtyperStage(
                ctx=ctx,
                gvcf=out_gvcf,
                output=self.output_vcf,
                interval=self.bed,
            ).add_to(dag, transfer_apply.terminal)

        # Call SVs
        if not self.skip_svs:
            svsolver_job = driver_job(
                ctx,
                [SVSolver(output=out_svs, vcf=out_svs_tmp)],
                interval=self.bed,
                threads=1,
                name="svsolver",
                task_name="sv-calling",
            )
            dag.add_job(svsolver_job, call.terminal)
            dag.add_job(
                rm_job([out_svs_tmp, str(out_svs_tmp) + ".tbi"], "rm-tmp-sv"),
                {svsolver_job},
            )
