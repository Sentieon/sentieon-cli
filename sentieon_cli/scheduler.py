"""Schedule jobs"""

from typing import Dict, Generator, Optional, Set

from .dag import DAG
from .job import Job
from .logging import get_logger

logger = get_logger(__name__)


class ThreadScheduler:
    """Schedule jobs as threads are available"""

    def __init__(
        self,
        dag: DAG,
        threads: int = 1,
        resources: Optional[Dict[str, int]] = None,
    ):
        self.dag = dag
        self.threads = threads
        self.available_threads = threads
        self.resources = {} if resources is None else resources

    def schedule(self) -> Generator[Set[Job], Optional[Job], None]:
        """Schedule a job for execution"""
        dag_gen = self.dag.update_dag()
        ready_jobs = dag_gen.send(None)

        while True:
            logger.debug("Ready jobs: %s", ready_jobs)

            scheduled_jobs: Set[Job] = set()
            for ready_job in ready_jobs:
                if self.available_threads - ready_job.threads >= 0 and all(
                    [
                        available - ready_job.resources.get(resource, 0) >= 0
                        for resource, available in self.resources.items()
                    ]
                ):
                    self.available_threads -= ready_job.threads
                    for resource, used in ready_job.resources.items():
                        if resource in self.resources:
                            self.resources[resource] -= used
                    scheduled_jobs.add(ready_job)

            ready_jobs -= scheduled_jobs

            finished_job = yield scheduled_jobs
            if isinstance(finished_job, Job):
                self.available_threads += finished_job.threads
                for resource, used in finished_job.resources.items():
                    if resource in self.resources:
                        self.resources[resource] += used
            ready_jobs.update(dag_gen.send(finished_job))
