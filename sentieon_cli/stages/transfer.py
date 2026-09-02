"""
Annotation transfer from a population VCF
"""

from dataclasses import dataclass
import pathlib
from typing import Dict, Iterable, List, Optional

from ..dag import DAG
from ..job import Job
from ..shard import Shard
from ..transfer import build_transfer_jobs
from .base import Stage, StageResult


@dataclass(frozen=True)
class TransferConfig:
    """The run-wide inputs the annotation transfer needs.

    Every pipeline that transfers annotations builds the same fan-out from
    the same four values, so they travel together instead of being threaded
    through each call site.
    """

    pop_vcf: pathlib.Path
    shards: List[Shard]
    pop_vcf_contigs: Dict[str, Optional[int]]
    fai_data: Dict[str, Dict[str, int]]


@dataclass(kw_only=True)
class TransferResult(StageResult):
    """The jobs and output of `TransferStage`"""

    shard_jobs: List[Job]
    concat_job: Job
    out_vcf: pathlib.Path


@dataclass(kw_only=True)
class TransferStage(Stage):
    """Transfer annotations from the population VCF onto a raw VCF.

    One job per genomic shard merges the shard's annotations, then a
    single concat job assembles ``out_vcf``. ``tag`` disambiguates the
    job names when a pipeline transfers more than once in a run.
    """

    config: TransferConfig
    raw_vcf: pathlib.Path
    out_vcf: pathlib.Path
    tag: Optional[str] = None

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> TransferResult:
        deps = set(upstream)

        shard_jobs, concat_job = build_transfer_jobs(
            self.out_vcf,
            self.config.pop_vcf,
            self.raw_vcf,
            self.ctx.tmp_dir,
            self.config.shards,
            self.config.pop_vcf_contigs,
            self.config.fai_data,
            self.ctx.dry_run,
            self.ctx.cores,
            tag=self.tag,
        )

        for job in shard_jobs:
            dag.add_job(job, deps)
        dag.add_job(concat_job, set(shard_jobs))

        return TransferResult(
            jobs=[*shard_jobs, concat_job],
            terminal={concat_job},
            shard_jobs=shard_jobs,
            concat_job=concat_job,
            out_vcf=self.out_vcf,
        )
