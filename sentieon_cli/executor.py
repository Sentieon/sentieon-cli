"""Execute jobs"""

import asyncio
import asyncio.subprocess
import signal
import sys
import time
from typing import Any, List, Tuple

from .job import Job
from .logging import get_logger
from .scheduler import ThreadScheduler
from .shell_pipeline import Context

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
        self.start_new_jobs = True
        self.running: List[
            Tuple[
                Job,
                Context,
                asyncio.Task[int],
                int,
            ]
        ] = []
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, signum, frame) -> None:
        logger.error("Termination signal detected. Terminating jobs.")
        self.start_new_jobs = False
        for running_job in self.running:
            for subcommand in running_job[1].commands:
                assert subcommand.proc
                subcommand.proc.send_signal(signal=signal.SIGTERM)
        time.sleep(10)
        for running_job in self.running:
            context = running_job[1]
            for subcommand in context.commands:
                assert subcommand.proc
                if subcommand.proc.returncode is None:
                    subcommand.proc.kill()
            asyncio.run(context.cleanup())

    async def run_job(self, job: Job) -> None:
        """Run a job"""
        cmd = job.shell
        logger.info("Running: %s", cmd)
        context = Context()
        start_time = time.monotonic_ns()
        try:
            proc = await cmd.run(
                context,
                stderr=sys.stderr,
            )
            self.running.append(
                (
                    job,
                    context,
                    asyncio.create_task(proc.wait()),
                    start_time,
                )
            )
        except Exception as e:
            logger.error("failed to start command: %s", job.shell)
            logger.error("Error: %s", str(e))
            await context.cleanup()

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
            for job, context, _task, start_time in finished_jobs:
                end_time = time.monotonic_ns()
                total_seconds = (end_time - start_time) / 1e9

                # Check if the command failed
                cmd_failed = False
                for subcommand in context.commands:
                    assert subcommand.proc
                    ret = (
                        await subcommand.proc.wait()
                    )  # Wait on all sub-commands
                    if ret != 0 and not subcommand.fail_ok:
                        logger.error(
                            f"Return code {ret} for sub-command: {subcommand}"
                        )
                        cmd_failed = True
                        if ret == -9:
                            logger.error(
                                "Sub-command received SIGKILL. "
                                "Possible out-of-memory?"
                            )

                if cmd_failed:
                    logger.error("Command failure: %s", job.shell)
                    self.jobs_with_errors.append(job)
                    self.start_new_jobs = False
                else:
                    logger.info(
                        f"Finished command in "
                        f"{total_seconds:.2f}: {job.shell}"
                    )
                await context.cleanup()

            if not self.start_new_jobs:
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
