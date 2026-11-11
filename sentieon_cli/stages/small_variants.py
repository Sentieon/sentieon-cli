"""
Small-variant calling, annotation transfer, and genotyping
"""

import argparse
from dataclasses import dataclass
import pathlib
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Set, Union

from ..dag import DAG
from ..driver import BaseAlgo, DNAModelApply, GVCFtyper, LuaScript
from ..job import Job
from .base import Stage, StageResult, driver_job
from .transfer import TransferConfig, TransferResult, TransferStage

DEFAULT_GQ_BANDS = "10-60/10"

# The `--gq_bands` argument of every pipeline that can emit a gVCF
GQ_BANDS_PARAM: Dict[str, Any] = {
    "help": (
        "The GQ bands of the reference blocks in the gVCF output, passed to "
        "DNAscope's --gq_bands. Only used with --gvcf. "
        "(default: %(default)s)"
    ),
    "default": DEFAULT_GQ_BANDS,
}

# The hidden `--pcr_indel_model` argument, which overrides `--pcr_free`
PCR_INDEL_MODEL_PARAM: Dict[str, Any] = {
    "help": argparse.SUPPRESS,
}


def dnascope_pcr_indel_model(
    pcr_free: bool, override: Optional[str] = None
) -> str:
    """The DNAscope `--pcr_indel_model` value for a run"""
    if override is not None:
        return override
    return "NONE" if pcr_free else "CONSERVATIVE"


@dataclass(kw_only=True)
class DriverResult(StageResult):
    """The job of a single-job `sentieon driver` stage"""

    job: Job


@dataclass(kw_only=True)
class DriverStage(Stage):
    """A single `sentieon driver` invocation, wired into the DAG.

    The driver arguments are passed through untouched; a stage that needs
    more structure than "run these algos" subclasses this one and fixes
    the defaults it cares about.
    """

    algos: List[BaseAlgo]
    name: str
    task_name: str
    threads: Optional[int] = None
    inputs: Optional[List[pathlib.Path]] = None
    interval: Optional[Union[pathlib.Path, str]] = None
    interval_padding: Optional[int] = None
    read_filter: Optional[List[str]] = None
    replace_rg: Optional[List[List[str]]] = None
    env: Optional[Dict[str, str]] = None

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> DriverResult:
        job = driver_job(
            self.ctx,
            self.algos,
            name=self.name,
            task_name=self.task_name,
            threads=self.threads,
            inputs=self.inputs,
            interval=self.interval,
            interval_padding=self.interval_padding,
            read_filter=self.read_filter,
            replace_rg=self.replace_rg,
            env=self.env,
        )
        dag.add_job(job, set(upstream))
        return DriverResult(jobs=[job], terminal={job}, job=job)


@dataclass(kw_only=True)
class DNAscopeStage(DriverStage):
    """The DNAscope-family calling pass.

    The caller supplies the algos in the order the driver should see
    them -- SNV/indel calling first, then a BND or PangenomeSV pass in
    the same pileup.
    """

    name: str = "dnascope"
    task_name: str = "variant-calling"


@dataclass(kw_only=True)
class ModelApplyResult(StageResult):
    """The job and output of `ModelApplyStage`"""

    job: Job
    output: pathlib.Path


@dataclass(kw_only=True)
class ModelApplyStage(Stage):
    """Genotype and filter a raw VCF with DNAModelApply"""

    model: pathlib.Path
    vcf: pathlib.Path
    output: pathlib.Path
    name: str = "model-apply"
    task_name: str = "model-apply"

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> ModelApplyResult:
        job = driver_job(
            self.ctx,
            [DNAModelApply(self.model, self.vcf, self.output)],
            name=self.name,
            task_name=self.task_name,
        )
        dag.add_job(job, set(upstream))
        return ModelApplyResult(
            jobs=[job],
            terminal={job},
            job=job,
            output=self.output,
        )


@dataclass(kw_only=True)
class GVCFtyperResult(StageResult):
    """The job and output of `GVCFtyperStage`"""

    job: Job
    output: pathlib.Path


@dataclass(kw_only=True)
class GVCFtyperStage(Stage):
    """Genotype a gVCF to produce a VCF"""

    gvcf: pathlib.Path
    output: pathlib.Path
    interval: Optional[pathlib.Path] = None
    dbsnp: Optional[pathlib.Path] = None
    name: str = "gvcftyper"
    task_name: str = "gvcftyper"

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> GVCFtyperResult:
        job = driver_job(
            self.ctx,
            [GVCFtyper(output=self.output, vcf=self.gvcf, dbsnp=self.dbsnp)],
            interval=self.interval,
            name=self.name,
            task_name=self.task_name,
        )
        dag.add_job(job, set(upstream))
        return GVCFtyperResult(
            jobs=[job],
            terminal={job},
            job=job,
            output=self.output,
        )


@dataclass(frozen=True)
class TransferSpec:
    """Ask `TransferApplyStage` for an annotation transfer"""

    config: TransferConfig
    out_vcf: pathlib.Path
    tag: Optional[str] = None


@dataclass(frozen=True)
class ApplySpec:
    """Ask `TransferApplyStage` for a DNAModelApply pass"""

    model: pathlib.Path
    output: pathlib.Path
    name: str = "model-apply"


@dataclass(frozen=True)
class PostSpec:
    """A post-processing script run before DNAModelApply"""

    script: pathlib.Path
    repeat_model: pathlib.Path
    out_vcf: pathlib.Path
    hap_table: Optional[pathlib.Path] = None
    interval: Optional[Union[pathlib.Path, str]] = None
    interval_padding: Optional[int] = None
    after: FrozenSet[Job] = frozenset()


@dataclass(kw_only=True)
class TransferApplyResult(StageResult):
    """The jobs and output of `TransferApplyStage`"""

    transfer: Optional[TransferResult]
    apply_job: Optional[Job]
    output: pathlib.Path
    post_job: Optional[Job] = None


@dataclass(kw_only=True)
class TransferApplyStage(Stage):
    """Post-process a raw VCF: transfer annotations, then apply the model.

    Either half may be omitted. When both run, DNAModelApply reads the
    transfer output and waits on the transfer alone -- the transfer's
    concat job already depends on the raw-VCF producer. Without a
    transfer, DNAModelApply waits on ``upstream``. ``post`` runs before
    DNAModelApply.
    """

    raw_vcf: pathlib.Path
    transfer: Optional[TransferSpec] = None
    post: Optional[PostSpec] = None
    apply: Optional[ApplySpec] = None

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> TransferApplyResult:
        deps = set(upstream)
        if self.transfer is None and self.apply is None and self.post is None:
            raise ValueError(
                "TransferApplyStage needs a transfer, a post pass or an apply"
            )

        jobs: List[Job] = []
        transfer_result: Optional[TransferResult] = None
        apply_in = self.raw_vcf
        apply_deps = deps
        if self.transfer is not None:
            transfer_result = TransferStage(
                ctx=self.ctx,
                config=self.transfer.config,
                raw_vcf=self.raw_vcf,
                out_vcf=self.transfer.out_vcf,
                tag=self.transfer.tag,
            ).add_to(dag, deps)
            jobs.extend(transfer_result.jobs)
            apply_in = self.transfer.out_vcf
            apply_deps = transfer_result.terminal

        post_job: Optional[Job] = None
        if self.post is not None:
            post_job = self.build_post_job(apply_in)
            dag.add_job(post_job, set(apply_deps) | self.post.after)
            jobs.append(post_job)
            apply_in = self.post.out_vcf
            apply_deps = {post_job}

        apply_job: Optional[Job] = None
        if self.apply is not None:
            apply_result = ModelApplyStage(
                ctx=self.ctx,
                model=self.apply.model,
                vcf=apply_in,
                output=self.apply.output,
                name=self.apply.name,
            ).add_to(dag, apply_deps)
            apply_job = apply_result.job
            jobs.append(apply_job)
            terminal: Set[Job] = {apply_job}
            output = self.apply.output
        elif post_job is not None:
            terminal = {post_job}
            output = apply_in
        else:
            assert transfer_result is not None
            terminal = transfer_result.terminal
            output = transfer_result.out_vcf

        return TransferApplyResult(
            jobs=jobs,
            terminal=terminal,
            transfer=transfer_result,
            apply_job=apply_job,
            output=output,
            post_job=post_job,
        )

    def build_post_job(self, in_vcf: pathlib.Path) -> Job:
        """Run the post-processing script on ``in_vcf``"""
        assert self.post is not None
        args: List[Union[pathlib.Path, str]] = [
            "--repeat-model",
            self.post.repeat_model,
        ]
        if self.post.hap_table is not None:
            args.extend(["--hap-table", self.post.hap_table])
        args.extend(["--lambda", "2", "--nofilt", in_vcf, self.post.out_vcf])
        return driver_job(
            self.ctx,
            [LuaScript(self.post.script, args)],
            name="post-lua",
            task_name="post-lua",
            interval=self.post.interval,
            interval_padding=self.post.interval_padding,
        )
