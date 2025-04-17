"""Execute jobs"""

import asyncio
import asyncio.subprocess
import sys
from typing import Any, List, Tuple

from .job import Job
from .logging import get_logger
from .scheduler import ThreadScheduler

logger = get_logger(__name__)


class BaseExecutor:
    """Execute jobs"""

    def __init__(self, scheduler: ThreadScheduler):
        self.scheduler = scheduler
        self.jobs_with_errors: List[Job] = []

    def execute(self) -> None:
        """Execute jobs from the DAG"""
        raise NotImplementedError


class DryRunExecutor(BaseExecutor):
    """Dry-run execution"""

    def run_job(self, job: Job) -> None:
        """Dry-run a job"""
        print(job.shell)

    def execute(self) -> None:
        scheduler_gen = self.scheduler.schedule()
        ready_jobs = scheduler_gen.send(None)
        for job in ready_jobs:
            self.run_job(job)

        while ready_jobs:
            finished_jobs = ready_jobs.copy()
            ready_jobs = {
                new_job
                for completed_job in finished_jobs
                for new_job in scheduler_gen.send(completed_job)
            }
            for job in ready_jobs:
                self.run_job(job)


class LocalExecutor(BaseExecutor):
    """Run jobs locally"""

    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self.running: List[
            Tuple[
                Job,
                asyncio.subprocess.Process,
                asyncio.Task[int],
            ]
        ] = []

    async def run_job(self, job: Job) -> None:
        """Run a job"""
        cmd = job.shell
        logger.info("Running: %s", cmd)
        proc = await asyncio.create_subprocess_shell(
            cmd,
            stdout=sys.stdout,
            stderr=sys.stderr,
            executable="/bin/bash",
        )
        self.running.append(
            (
                job,
                proc,
                asyncio.create_task(proc.wait()),
            )
        )

    def execute(self) -> None:
        """Execute jobs from the DAG"""
        asyncio.run(self._execute())

    async def _execute(self) -> None:
        """Execute jobs from the DAG"""
        self.jobs_with_errors: List[Job] = []
        scheduler_gen = self.scheduler.schedule()
        ready_jobs = scheduler_gen.send(None)
        for job in ready_jobs:
            await self.run_job(job)

        while self.running:
            done, _running = await asyncio.wait(
                [job[2] for job in self.running],
                return_when=asyncio.FIRST_COMPLETED,
            )

            finished_jobs = [
                self.running.pop(i)
                for i in reversed(range(len(self.running)))
                if self.running[i][2] in done
            ]

            # Check job execution
            for job, proc, _task in finished_jobs:
                if proc.returncode != 0 and not job.fail_ok:
                    logger.error("Error running command, '%s'", job.shell)
                    self.jobs_with_errors.append(job)

            if self.jobs_with_errors:
                # Don't start new jobs
                continue

            ready_jobs = {
                new_job
                for completed_job in finished_jobs
                for new_job in scheduler_gen.send(completed_job[0])
            }

            # Run the ready jobs
            for job in ready_jobs:
                await self.run_job(job)
