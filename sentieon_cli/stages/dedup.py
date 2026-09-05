"""
Duplicate marking with LocusCollector and Dedup
"""

from dataclasses import dataclass, field
import pathlib
from typing import Iterable, List, Optional

from ..dag import DAG
from ..driver import BaseAlgo, Dedup, LocusCollector
from ..job import Job
from .base import Stage, StageResult, driver_job


@dataclass(kw_only=True)
class DedupResult(StageResult):
    """The jobs and output of `DedupStage`"""

    lc_job: Job
    dedup_job: Job
    output: pathlib.Path


@dataclass(kw_only=True)
class DedupStage(Stage):
    """Mark duplicates: LocusCollector collects the scores, Dedup applies
    them.

    Both passes read the same input files with the same read filters, so
    the two jobs only differ in their algo. Metrics that are cheap to
    collect before duplicate marking can ride along on the LocusCollector
    pass through ``lc_extra_algos``.
    """

    inputs: List[pathlib.Path]
    output: pathlib.Path
    score_file: pathlib.Path
    tag: Optional[str] = None
    consensus: bool = False
    rmdup: bool = False
    cram_write_options: Optional[str] = None
    dedup_metrics: Optional[pathlib.Path] = None
    read_filters: List[str] = field(default_factory=list)
    lc_extra_algos: List[BaseAlgo] = field(default_factory=list)
    threads: Optional[int] = None
    task_name: str = "dedup"

    def _name(self, base: str) -> str:
        return base if self.tag is None else f"{base}-{self.tag}"

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> DedupResult:
        deps = set(upstream)
        read_filter = self.read_filters or None

        lc_job = driver_job(
            self.ctx,
            [
                LocusCollector(self.score_file, consensus=self.consensus),
                *self.lc_extra_algos,
            ],
            inputs=self.inputs,
            read_filter=read_filter,
            name=self._name("locuscollector"),
            task_name=self.task_name,
            threads=self.threads,
        )
        dedup_job = driver_job(
            self.ctx,
            [
                Dedup(
                    self.output,
                    self.score_file,
                    cram_write_options=self.cram_write_options,
                    metrics=self.dedup_metrics,
                    rmdup=self.rmdup,
                )
            ],
            inputs=self.inputs,
            read_filter=read_filter,
            name=self._name("dedup"),
            task_name=self.task_name,
            threads=self.threads,
        )

        dag.add_job(lc_job, deps)
        dag.add_job(dedup_job, {lc_job})

        return DedupResult(
            jobs=[lc_job, dedup_job],
            terminal={dedup_job},
            lc_job=lc_job,
            dedup_job=dedup_job,
            output=self.output,
        )
