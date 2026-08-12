"""
A directed acyclic graph of jobs to execute
"""

from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Optional, Set

from .exceptions import DagExecutionError
from .job import Job
from .logging import get_logger

logger = get_logger(__name__)


class DAG:
    """A directed acyclic graph of jobs"""

    def __init__(self) -> None:
        # waiting_jobs = {job: {dependencies}}
        self.waiting_jobs: Dict[Job, Set[Job]] = {}
        # map dependencies to waiting jobs
        # planned_jobs = {dependency: [downstream_jobs]}
        self.planned_jobs: Dict[Job, List[Job]] = {}
        self.ready_jobs: Dict[Job, None] = {}
        self.finished_jobs: List[Job] = []

    def add_job(
        self, job: Job, dependencies: Optional[Iterable[Job]] = None
    ) -> None:
        """Add a job to the DAG.

        ``dependencies`` may be any iterable of jobs already in the DAG;
        it is consumed once and stored as a private set.
        """
        deps: Set[Job] = (
            set(dependencies) if dependencies is not None else set()
        )
        if job in self.waiting_jobs or job in self.ready_jobs:
            raise ValueError(
                f"{job} is already in the DAG; jobs are identified by their "
                f"shell pipeline, so the same pipeline cannot be added twice."
            )
        for dependency in deps:
            if (
                dependency not in self.waiting_jobs
                and dependency not in self.ready_jobs
            ):
                raise ValueError(
                    f"Dependency '{dependency}' is not in the DAG"
                )

        if deps:
            self.waiting_jobs[job] = deps
            for dependency in deps:
                downstream = self.planned_jobs.setdefault(dependency, [])
                downstream.append(job)
        else:
            self.ready_jobs[job] = None

    def mark_finished(self, job: Job) -> Dict[Job, None]:
        """Record a finished job and return the jobs it unblocked.

        The returned mapping is used as an ordered set (values are ``None``).
        Raises :class:`~sentieon_cli.exceptions.DagExecutionError` if ``job``
        was not currently ready.
        """
        if job not in self.ready_jobs:
            raise DagExecutionError(
                f"Finished job '{job}' was not ready for execution"
            )

        logger.debug("Marking finished: %s", job)
        del self.ready_jobs[job]
        self.finished_jobs.append(job)

        new_ready_jobs: Dict[Job, None] = {}
        for dependency in self.planned_jobs.get(job, []):
            upstream = self.waiting_jobs[dependency]
            upstream.remove(job)
            if len(upstream) < 1:
                self.ready_jobs[dependency] = None
                new_ready_jobs[dependency] = None
                del self.waiting_jobs[dependency]
        return new_ready_jobs

    def skip_satisfied(self, is_satisfied: Callable[[Job], bool]) -> List[Job]:
        """Mark already-satisfied jobs finished before a run (resume).

        Walks the ready frontier: a job is skipped iff ``is_satisfied(job)``
        returns True and every one of its dependencies was itself skipped.
        A job behind an unsatisfied ancestor is never offered to the
        predicate -- it only becomes ready once all its dependencies are
        finished, and during this walk "finished" means "skipped".

        Skipped jobs land in :attr:`finished_jobs`, so an executor never
        sees them; unsatisfied jobs are left in place untouched. Call this
        before handing the DAG to a scheduler.

        Returns the skipped jobs in topological order.
        """
        skipped: List[Job] = []
        frontier: List[Job] = list(self.ready_jobs)
        i = 0
        while i < len(frontier):
            job = frontier[i]
            i += 1
            if is_satisfied(job):
                skipped.append(job)
                frontier.extend(self.mark_finished(job))
        return skipped


def has_all_outputs(job: Job) -> bool:
    """The default resume predicate for :meth:`DAG.skip_satisfied`.

    True iff ``job`` declares at least one output and every declared
    output exists on the local filesystem. Jobs with no declared outputs
    always return False -- they cannot be verified, so they must rerun.
    Zero-byte files count as present; a dangling symlink does not
    (``Path.exists()`` follows links) -- matching
    :class:`~sentieon_cli.storage.LocalStorageProvider` semantics.
    """
    return bool(job.outputs) and all(path.exists() for path in job.outputs)
