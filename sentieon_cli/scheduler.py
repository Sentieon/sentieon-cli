"""Schedule jobs"""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Set

from .dag import DAG
from .exceptions import DagExecutionError
from .job import Job
from .logging import get_logger

logger = get_logger(__name__)


class BaseScheduler(ABC):
    """Selects which ready jobs an executor may run, subject to a policy.

    An executor drives a scheduler with two calls:

    1. :meth:`start` once, to get the first batch of jobs to run.
    2. :meth:`job_finished` for each job as it completes, to get the jobs
       that its completion unblocked (possibly an empty set).

    Contract / invariants:

    * :meth:`start` is called exactly once, before any :meth:`job_finished`.
    * Every job passed to :meth:`job_finished` was previously returned by
      :meth:`start` or :meth:`job_finished`, and is reported at most once.
    * A scheduler instance drives a single run; do not reuse it.
    * :meth:`start` may raise
      :class:`~sentieon_cli.exceptions.DagExecutionError` if the DAG
      contains a job that can never be scheduled.
    * The executor stops when it has no running jobs and the scheduler
      returns no new ones.
    """

    @abstractmethod
    def start(self) -> Set[Job]:
        """Return the first batch of jobs to run."""

    @abstractmethod
    def job_finished(self, job: Job) -> Set[Job]:
        """Record a finished job and return the jobs it unblocked."""


class ThreadScheduler(BaseScheduler):
    """Schedule jobs as threads are available"""

    def __init__(
        self,
        dag: DAG,
        threads: int = 1,
        resources: Optional[Dict[str, int]] = None,
    ) -> None:
        self.dag = dag
        self.threads = threads
        self.available_threads = threads
        self.resources = {} if resources is None else resources
        # Jobs that are ready to run but have not yet been scheduled.
        self._ready: Dict[Job, None] = {}

    def _check_feasible(self) -> None:
        """Ensure every job can eventually be scheduled.

        A job that requests more threads than the total budget, or more of a
        managed resource than exists, can never be released and would
        otherwise stall the run silently.
        """
        jobs = list(self.dag.waiting_jobs) + list(self.dag.ready_jobs)
        for job in jobs:
            if job.threads > self.threads:
                raise DagExecutionError(
                    f"{job} requests {job.threads} threads but the scheduler "
                    f"budget is only {self.threads}"
                )
            for resource, needed in job.resources.items():
                capacity = self.resources.get(resource)
                if capacity is not None and needed > capacity:
                    raise DagExecutionError(
                        f"{job} requests {needed} of resource '{resource}' "
                        f"but the scheduler capacity is only {capacity}"
                    )

    def _fit(self) -> Set[Job]:
        """Reserve budget for as many ready jobs as fit; return that set."""
        logger.debug("Ready pool: %s", self._ready)

        scheduled_jobs: Set[Job] = set()
        for ready_job in self._ready:
            if self.available_threads - ready_job.threads >= 0 and all(
                available - ready_job.resources.get(resource, 0) >= 0
                for resource, available in self.resources.items()
            ):
                self.available_threads -= ready_job.threads
                for resource, used in ready_job.resources.items():
                    if resource in self.resources:
                        self.resources[resource] -= used
                scheduled_jobs.add(ready_job)

        for job in scheduled_jobs:
            del self._ready[job]
        return scheduled_jobs

    def start(self) -> Set[Job]:
        self._check_feasible()
        self._ready = dict(self.dag.ready_jobs)
        return self._fit()

    def job_finished(self, job: Job) -> Set[Job]:
        self.available_threads += job.threads
        for resource, used in job.resources.items():
            if resource in self.resources:
                self.resources[resource] += used
        self._ready.update(self.dag.mark_finished(job))
        return self._fit()
