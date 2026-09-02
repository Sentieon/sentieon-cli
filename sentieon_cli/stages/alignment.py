"""
Read alignment and read extraction

Alignment jobs are the roots of every pipeline DAG: they read the run's
input files and depend on nothing. Because a few pipelines need the jobs
before they insert them, these stages split the work in two -- ``build``
constructs the jobs, ``add_to`` inserts what ``build`` returned.
"""

from abc import abstractmethod
from dataclasses import dataclass, field
import itertools
import logging
import pathlib
import shutil
from typing import Iterable, List, Optional, Union

from .. import command_strings as cmds
from ..dag import DAG
from ..driver import ReadWriter
from ..job import Job
from .base import Stage, StageResult, driver_job

# The `sentieon util sort` arguments every aligner defaults to
DEFAULT_UTIL_SORT_ARGS = "--cram_write_options version=3.0,compressor=rans"


def aln_suffix(bam_format: bool) -> str:
    """The file extension for an aligned output file"""
    return "bam" if bam_format else "cram"


def find_unzip(logger: logging.Logger, level: int = logging.WARNING) -> str:
    """The decompression tool to read gzipped fastq with.

    ``igzip`` is considerably faster but is not always installed, so fall
    back to ``gzip``. Callers differ on how loudly that is worth
    reporting, hence ``level``.
    """
    unzip = "igzip"
    if not shutil.which(unzip):
        logger.log(
            level,
            "igzip is recommended for decompression, but is not "
            "available. Falling back to gzip.",
        )
        unzip = "gzip"
    return unzip


@dataclass(kw_only=True)
class AlignResult(StageResult):
    """The jobs and outputs of an `AlignmentStage`.

    * ``outputs`` -- the aligned files, in input order.
    * ``cleanup_paths`` -- files a later ``rm`` job should delete once the
      alignment has been consumed; empty when the stage writes only
      final outputs.
    """

    outputs: List[pathlib.Path]
    cleanup_paths: List[pathlib.Path]


@dataclass(kw_only=True)
class AlignmentStage(Stage):
    """A group of alignment jobs, all of them DAG roots.

    Subclasses implement :meth:`build`, which constructs the jobs and
    inserts nothing. :meth:`add_to` then adds each of them with the same
    dependencies, since no alignment job depends on another.
    """

    @abstractmethod
    def build(self) -> AlignResult:
        """Construct this stage's jobs without touching a DAG"""

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> AlignResult:
        deps = set(upstream)
        result = self.build()
        for job in result.jobs:
            dag.add_job(job, deps)
        return result


@dataclass(kw_only=True)
class BwaRealignStage(AlignmentStage):
    """Re-align BAM/CRAM/uBAM/uCRAM input to the reference with bwa.

    The @RG lines of each input are read at build time and written to a
    header file in the run's temporary directory, which the alignment
    command passes to `samtools fastq`.

    When duplicate marking follows, the alignment is an intermediate: it
    is written to the temporary directory as a lightly compressed BAM
    instead of to the output directory.
    """

    inputs: List[pathlib.Path]
    model_bundle: pathlib.Path
    bam_format: bool = False
    duplicate_marking: str = "markdup"
    input_ref: Optional[pathlib.Path] = None
    collate: bool = False
    bwa_args: str = ""
    bwa_k: str = "20000000"
    util_sort_args: str = DEFAULT_UTIL_SORT_ARGS

    def build(self) -> AlignResult:
        suffix = aln_suffix(self.bam_format)
        util_sort_args = self.util_sort_args
        dedup_follows = self.duplicate_marking != "none"
        if dedup_follows:
            suffix = "bam"
            util_sort_args += " --bam_compression 1 "

        jobs: List[Job] = []
        outputs: List[pathlib.Path] = []
        cleanup_paths: List[pathlib.Path] = []
        for i, input_aln in enumerate(self.inputs):
            out_aln = pathlib.Path(
                str(self.ctx.output_vcf).replace(
                    ".vcf.gz", f"_bwa_sorted_{i}.{suffix}"
                )
            )
            if dedup_follows:
                out_aln = self.ctx.tmp_dir.joinpath(f"bwa_sorted_{i}.{suffix}")
            rg_lines = cmds.get_rg_lines(input_aln, self.ctx.dry_run)
            rg_header = self.ctx.tmp_dir.joinpath(f"input_{i}.hdr")
            with open(rg_header, "w", encoding="utf-8") as rg_fh:
                for line in rg_lines:
                    print(line, file=rg_fh)
            jobs.append(
                Job(
                    cmds.cmd_samtools_fastq_bwa(
                        out_aln,
                        input_aln,
                        self.ctx.reference,
                        self.model_bundle,
                        self.ctx.cores,
                        rg_header,
                        self.input_ref,
                        collate=self.collate,
                        bwa_args=self.bwa_args,
                        bwa_k=self.bwa_k,
                        util_sort_args=util_sort_args,
                    ),
                    f"bam-align-{i}",
                    self.ctx.cores,
                    task_name="alignment",
                )
            )
            outputs.append(out_aln)
            cleanup_paths.append(out_aln)
            cleanup_paths.append(pathlib.Path(str(out_aln) + ".bai"))
            if suffix == "cram":
                cleanup_paths.append(pathlib.Path(str(out_aln) + ".crai"))

        return AlignResult(
            jobs=jobs,
            terminal=set(jobs),
            outputs=outputs,
            cleanup_paths=cleanup_paths,
        )


@dataclass(kw_only=True)
class BwaFastqStage(AlignmentStage):
    """Align short-read fastq pairs to the reference with bwa.

    With more than one NUMA node the reads of each pair are split across
    the nodes: one job per (fastq pair, node), each pinned with `taskset`
    and taking a share of the cores and a node token, so the scheduler
    runs exactly one alignment per node.
    """

    r1_fastq: List[pathlib.Path]
    r2_fastq: List[pathlib.Path]
    readgroups: List[str]
    model_bundle: pathlib.Path
    numa_nodes: List[str] = field(default_factory=list)
    bam_format: bool = False
    duplicate_marking: str = "markdup"
    unzip: str = "gzip"
    bwa_args: str = ""
    bwa_k: str = "20000000"
    util_sort_args: str = DEFAULT_UTIL_SORT_ARGS

    def build(self) -> AlignResult:
        n_alignment_jobs = max(1, len(self.numa_nodes))
        split_alignment = len(self.numa_nodes) > 1
        suffix = aln_suffix(self.bam_format)
        util_sort_args = self.util_sort_args
        dedup_follows = self.duplicate_marking != "none"
        if dedup_follows:
            suffix = "bam"
            util_sort_args += " --bam_compression 1 "

        jobs: List[Job] = []
        outputs: List[pathlib.Path] = []
        cleanup_paths: List[pathlib.Path] = []
        for i, (r1, r2, rg) in enumerate(
            itertools.zip_longest(
                self.r1_fastq, self.r2_fastq, self.readgroups
            )
        ):
            for j in range(n_alignment_jobs):
                out_aln = pathlib.Path(
                    str(self.ctx.output_vcf).replace(
                        ".vcf.gz", f"_bwa_sorted_fq_{i}_{j}.{suffix}"
                    )
                )
                if dedup_follows:
                    out_aln = self.ctx.tmp_dir.joinpath(
                        f"bwa_sorted_fq_{i}_{j}.{suffix}"
                    )
                numa = self.numa_nodes[j] if split_alignment else None
                split = f"{j}/{n_alignment_jobs}" if split_alignment else None
                split_cores = max(1, int(self.ctx.cores / n_alignment_jobs))
                jobs.append(
                    Job(
                        cmds.cmd_fastq_bwa(
                            out_aln,
                            r1,
                            r2,
                            rg,
                            self.ctx.reference,
                            self.model_bundle,
                            split_cores,
                            self.unzip,
                            self.bwa_args,
                            self.bwa_k,
                            util_sort_args,
                            numa,
                            split,
                        ),
                        f"bam-align-{i}-{j}",
                        split_cores,
                        resources={f"node{j}": 1},
                        task_name="alignment",
                    )
                )
                outputs.append(out_aln)
                cleanup_paths.append(out_aln)
                cleanup_paths.append(pathlib.Path(str(out_aln) + ".bai"))
                if suffix == "cram":
                    cleanup_paths.append(pathlib.Path(str(out_aln) + ".crai"))

        return AlignResult(
            jobs=jobs,
            terminal=set(jobs),
            outputs=outputs,
            cleanup_paths=cleanup_paths,
        )


@dataclass(kw_only=True)
class Minimap2RealignStage(AlignmentStage):
    """Re-align long-read BAM/CRAM/uBAM/uCRAM input with minimap2.

    ``minimap2_model`` overrides the model the aligner picks up from the
    bundle; the pangenome pipelines use it to select the long-read model.
    """

    inputs: List[pathlib.Path]
    model_bundle: pathlib.Path
    sample_name: str
    bam_format: bool = False
    input_ref: Optional[pathlib.Path] = None
    fastq_taglist: str = "*"
    minimap2_args: str = "-YL"
    util_sort_args: str = DEFAULT_UTIL_SORT_ARGS
    minimap2_model: Optional[Union[pathlib.Path, str]] = None

    def build(self) -> AlignResult:
        suffix = aln_suffix(self.bam_format)
        jobs: List[Job] = []
        outputs: List[pathlib.Path] = []
        for i, input_aln in enumerate(self.inputs):
            out_aln = pathlib.Path(
                str(self.ctx.output_vcf).replace(
                    ".vcf.gz", f"_mm2_sorted_{i}.{suffix}"
                )
            )
            rg_lines = cmds.get_rg_lines(input_aln, self.ctx.dry_run)
            jobs.append(
                Job(
                    cmds.cmd_samtools_fastq_minimap2(
                        out_aln,
                        input_aln,
                        self.ctx.reference,
                        self.model_bundle,
                        self.ctx.cores,
                        rg_lines,
                        self.sample_name,
                        self.input_ref,
                        self.fastq_taglist,
                        self.minimap2_args,
                        self.util_sort_args,
                        self.minimap2_model,
                    ),
                    f"bam-realign-{i}",
                    self.ctx.cores,
                    task_name="alignment",
                )
            )
            outputs.append(out_aln)

        return AlignResult(
            jobs=jobs,
            terminal=set(jobs),
            outputs=outputs,
            cleanup_paths=[],
        )


@dataclass(kw_only=True)
class Minimap2FastqStage(AlignmentStage):
    """Align long-read fastq to the reference with minimap2"""

    fastq: List[pathlib.Path]
    readgroups: List[str]
    model_bundle: pathlib.Path
    bam_format: bool = False
    unzip: str = "gzip"
    minimap2_args: str = "-YL"
    util_sort_args: str = DEFAULT_UTIL_SORT_ARGS

    def build(self) -> AlignResult:
        suffix = aln_suffix(self.bam_format)
        jobs: List[Job] = []
        outputs: List[pathlib.Path] = []
        for i, (fq, rg) in enumerate(zip(self.fastq, self.readgroups)):
            out_aln = pathlib.Path(
                str(self.ctx.output_vcf).replace(
                    ".vcf.gz", f"_mm2_sorted_fq_{i}.{suffix}"
                )
            )
            jobs.append(
                Job(
                    cmds.cmd_fastq_minimap2(
                        out_aln,
                        fq,
                        rg,
                        self.ctx.reference,
                        self.model_bundle,
                        self.ctx.cores,
                        self.unzip,
                        self.minimap2_args,
                        self.util_sort_args,
                    ),
                    f"align-{i}",
                    self.ctx.cores,
                    task_name="alignment",
                )
            )
            outputs.append(out_aln)

        return AlignResult(
            jobs=jobs,
            terminal=set(jobs),
            outputs=outputs,
            cleanup_paths=[],
        )


@dataclass(kw_only=True)
class BwaExtractStage(AlignmentStage):
    """Align short-read fastq with bwa and extract the reads a pangenome
    pipeline needs, in a single pass.

    ``readgroup`` is the rendered `@RG` string; the pangenome pipelines
    derive it differently from the input readgroup.
    """

    output_bam: pathlib.Path
    output_fastq: pathlib.Path
    r1_fastq: List[pathlib.Path]
    r2_fastq: List[pathlib.Path]
    readgroup: str
    extract_model: pathlib.Path
    bwa_model: pathlib.Path
    unzip: str = "gzip"
    name: str = "bwa-extract"
    task_name: str = "alignment"

    def build(self) -> AlignResult:
        job = Job(
            cmds.cmd_bwa_extract(
                self.output_bam,
                self.output_fastq,
                self.ctx.reference,
                self.r1_fastq,
                self.r2_fastq,
                self.readgroup,
                self.extract_model,
                self.bwa_model,
                self.ctx.cores,
                unzip=self.unzip,
            ),
            self.name,
            self.ctx.cores,
            task_name=self.task_name,
        )
        return AlignResult(
            jobs=[job],
            terminal={job},
            outputs=[self.output_bam, self.output_fastq],
            cleanup_paths=[],
        )


@dataclass(kw_only=True)
class ReadWriterStage(AlignmentStage):
    """Rewrite alignments with `sentieon driver --algo ReadWriter`.

    With several inputs the driver merges them; with an ``interval`` it
    extracts the reads overlapping that region.
    """

    inputs: List[pathlib.Path]
    output: pathlib.Path
    name: str
    task_name: str
    threads: Optional[int] = None
    interval: Optional[Union[pathlib.Path, str]] = None

    def build(self) -> AlignResult:
        job = driver_job(
            self.ctx,
            [ReadWriter(self.output)],
            inputs=self.inputs,
            interval=self.interval,
            name=self.name,
            task_name=self.task_name,
            threads=self.threads,
        )
        return AlignResult(
            jobs=[job],
            terminal={job},
            outputs=[self.output],
            cleanup_paths=[],
        )
