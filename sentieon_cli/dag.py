"""
A directed acyclic graph of jobs to execute
"""

from __future__ import annotations

from typing import Dict, Generator, List, Optional, Set

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
        self.ready_jobs: Set[Job] = set()
        self.finished_jobs: List[Job] = []

    def add_job(self, job: Job, dependencies: Optional[Set[Job]] = None):
        """Add a job to the DAG"""
        if dependencies:
            for dependency in dependencies:
                assert (
                    dependency in self.waiting_jobs
                    or dependency in self.ready_jobs
                )

        if isinstance(dependencies, set) and len(dependencies) > 0:
            self.waiting_jobs[job] = dependencies
            for dependency in dependencies:
                downstream = self.planned_jobs.setdefault(dependency, [])
                downstream.append(job)
        else:
            self.ready_jobs.add(job)

    def update_dag(
        self,
    ) -> Generator[Set[Job], Optional[Job], None]:
        """Update the DAG with a finished job"""
        finished_job = yield self.ready_jobs.copy()

        while True:
            logger.debug("Waiting jobs: %s", self.waiting_jobs)
            logger.debug("Planned jobs: %s", self.planned_jobs)
            logger.debug("Ready jobs: %s", self.ready_jobs)
            logger.debug("Finished jobs: %s", self.finished_jobs)
            logger.debug("Newly finished: %s", finished_job)

            new_ready_jobs: set[Job] = set()
            if isinstance(finished_job, Job):
                if finished_job not in self.ready_jobs:
                    raise DagExecutionError(
                        f"Finished job '{finished_job}' was not ready for "
                        "execution"
                    )

                self.ready_jobs.remove(finished_job)
                self.finished_jobs.append(finished_job)
                for dependency in self.planned_jobs.get(finished_job, []):
                    upstream = self.waiting_jobs[dependency]
                    upstream.remove(finished_job)
                    if len(upstream) < 1:
                        self.ready_jobs.add(dependency)
                        new_ready_jobs.add(dependency)
                        del self.waiting_jobs[dependency]
            finished_job = yield new_ready_jobs
