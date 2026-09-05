"""
Building blocks for pipeline stages

A *stage* is a reusable piece of a pipeline: it knows how to build the jobs
for one operation (dedup, metrics, small-variant calling, ...) and how to
wire them into a :class:`~sentieon_cli.dag.DAG`. Stages are configured with
plain dataclass fields plus a :class:`StageContext` of run-wide settings.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
import pathlib
from typing import Iterable, List, Optional, Sequence, Set, Union

from ..dag import DAG
from ..driver import BaseAlgo, Driver
from ..job import Job
from ..shell_pipeline import Command, Pipeline


@dataclass(frozen=True)
class StageContext:
    """Run-wide settings every stage needs."""

    reference: pathlib.Path
    output_vcf: pathlib.Path
    tmp_dir: pathlib.Path
    cores: int
    dry_run: bool
    skip_version_check: bool


@dataclass(kw_only=True)
class StageResult:
    """The jobs a stage added to the DAG.

    * ``jobs`` -- every job the stage inserted, in insertion order. It
      exists for the occasional edge that has to attach to a job in the
      middle of a stage rather than to its output.
    * ``terminal`` -- the stage's sink job(s), the ones that carry its
      declared outputs. Downstream stages depend on these.

    Subclasses add named fields for the files the stage produces.
    """

    jobs: List[Job]
    terminal: Set[Job]


@dataclass(kw_only=True)
class Stage(ABC):
    """A reusable group of jobs that inserts itself into a DAG.

    Subclasses are dataclasses (``@dataclass(kw_only=True)``) that add
    their own configuration fields and return a :class:`StageResult`
    subclass naming the outputs they produce.

    The :meth:`add_to` contract:

    * ``upstream`` is any iterable of jobs already in the DAG; it is
      consumed once into a set.
    * Every *entry* job -- an in-stage job with no in-stage dependency --
      is inserted with that set as its dependencies. Edges between the
      stage's own jobs are added on top of that.
    * Insertion is topological (a job is added only after the jobs it
      depends on) and deterministic: the same stage and DAG produce the
      same jobs, in the same order, every time.
    * The jobs are created inside ``add_to``; constructing a stage has no
      side effects and assigns no job ids.
    * ``add_to`` may be called only once per stage instance. Job identity
      is the shell pipeline, so a second call raises in
      :meth:`~sentieon_cli.dag.DAG.add_job`.
    """

    ctx: StageContext

    @abstractmethod
    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> StageResult:
        """Add this stage's jobs to ``dag``, downstream of ``upstream``"""


def driver_job(
    ctx: StageContext,
    algos: Sequence[BaseAlgo],
    *,
    name: str,
    task_name: str,
    threads: Optional[int] = None,
    inputs: Optional[List[pathlib.Path]] = None,
    interval: Optional[Union[pathlib.Path, str]] = None,
    interval_padding: Optional[int] = None,
    read_filter: Optional[List[str]] = None,
    replace_rg: Optional[List[List[str]]] = None,
) -> Job:
    """A job running `sentieon driver` with ``algos``.

    ``threads`` defaults to the run's core count; the driver arguments are
    passed through untouched so the command matches what the pipelines
    build by hand.
    """
    driver = Driver(
        reference=ctx.reference,
        thread_count=ctx.cores,
        interval=interval,
        interval_padding=interval_padding,
        read_filter=read_filter,
        replace_rg=replace_rg,
        input=inputs,
    )
    for algo in algos:
        driver.add_algo(algo)
    return Job(
        Pipeline(Command(*driver.build_cmd())),
        name,
        ctx.cores if threads is None else threads,
        task_name=task_name,
    )


def rm_job(paths: Iterable[Union[pathlib.Path, str]], name: str) -> Job:
    """A job removing intermediate files, tolerating missing ones"""
    return Job(
        Pipeline(Command("rm", *[str(p) for p in paths], fail_ok=True)),
        name,
        0,
        task_name="cleanup",
    )
