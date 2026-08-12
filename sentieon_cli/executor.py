"""Execute jobs"""

import asyncio
import contextlib
import signal
import sys
import threading
import time
from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, List, Optional, Tuple

from .job import Job
from .logging import get_logger
from .scheduler import BaseScheduler
from .shell_pipeline import Context
from .storage import LocalStorageProvider, StorageProvider

logger = get_logger(__name__)


def _signal_procs(context: Context, sig: int) -> None:
    """Send a signal to every live sub-process in a context."""
    for subcommand in context.commands:
        proc = subcommand.proc
        if proc is not None and proc.returncode is None:
            try:
                proc.send_signal(sig)
            except ProcessLookupError:
                pass


def _kill_survivors(context: Context) -> None:
    """SIGKILL any sub-process that is still running."""
    for subcommand in context.commands:
        proc = subcommand.proc
        if proc is not None and proc.returncode is None:
            try:
                proc.kill()
            except ProcessLookupError:
                pass


async def _cleanup_quietly(context: Context) -> None:
    """Release a context, logging -- not raising -- a surfaced launch failure.

    A proc-sub inner command that fails to spawn surfaces its OSError from
    ``context.cleanup()``'s gather. On the abort and shutdown paths the owning
    job is already recorded (or the whole run is being torn down), so a stored
    OSError must not escape ``execute()``; log it and release the rest. A
    non-OSError stays a bug and propagates, matching ``_drive``'s per-job
    handling.
    """
    try:
        await context.cleanup()
    except OSError as exc:
        logger.error("error releasing resources during teardown: %s", exc)


class BaseExecutor(ABC):
    """Execute jobs"""

    def __init__(self, scheduler: BaseScheduler):
        self.scheduler = scheduler
        self.jobs_with_errors: List[Job] = []

    @abstractmethod
    def execute(self) -> None:
        """Execute jobs from the DAG"""


class DryRunExecutor(BaseExecutor):
    """Dry-run execution"""

    def run_job(self, job: Job) -> None:
        """Dry-run a job"""
        print(job.shell)

    def execute(self) -> None:
        ready_jobs = self.scheduler.start()
        for job in ready_jobs:
            self.run_job(job)

        while ready_jobs:
            finished_jobs = ready_jobs.copy()
            ready_jobs = {
                new_job
                for completed_job in finished_jobs
                for new_job in self.scheduler.job_finished(completed_job)
            }
            for job in ready_jobs:
                self.run_job(job)


class AsyncExecutor(BaseExecutor, ABC):
    """Base for executors that run an asyncio event loop.

    Provides the ``execute()`` entry point and opt-in SIGINT/SIGTERM handling
    with graceful shutdown; subclasses implement :meth:`_drive` (the run loop)
    and :meth:`_shutdown` (teardown after an interrupt).

    Set ``install_signal_handlers=True`` to install handlers (on the main
    thread only) for the duration of a run; the previous handlers are
    restored afterwards, so a caller's signal state is left untouched.
    """

    def __init__(
        self,
        scheduler: BaseScheduler,
        *,
        install_signal_handlers: bool = False,
    ) -> None:
        super().__init__(scheduler)
        self.install_signal_handlers = install_signal_handlers
        self.start_new_jobs = True

    def _install_signal_handlers(
        self,
        loop: asyncio.AbstractEventLoop,
        stop_event: asyncio.Event,
    ) -> Callable[[], None]:
        """Install SIGINT/SIGTERM handlers on the running loop.

        Returns a callable that removes the handlers and restores the
        previous ones. A no-op is returned when handlers are disabled, or
        when we are not on the main thread (Python can only manage signals
        there, so we skip rather than crash).
        """
        if not self.install_signal_handlers:
            return lambda: None
        if threading.current_thread() is not threading.main_thread():
            logger.debug("Signal handlers require the main thread; skipping.")
            return lambda: None

        def on_signal(signum: int) -> None:
            logger.error(
                "Received %s; terminating jobs.",
                signal.Signals(signum).name,
            )
            self.start_new_jobs = False
            stop_event.set()

        signums = (signal.SIGINT, signal.SIGTERM)
        previous: Dict[int, Any] = {
            signum: signal.getsignal(signum) for signum in signums
        }
        installed: List[int] = []
        for signum in signums:
            try:
                loop.add_signal_handler(signum, on_signal, signum)
                installed.append(signum)
            except (NotImplementedError, RuntimeError, ValueError) as exc:
                logger.debug(
                    "Could not install handler for %s: %s", signum, exc
                )

        def restore() -> None:
            for signum in installed:
                try:
                    loop.remove_signal_handler(signum)
                except (RuntimeError, ValueError):
                    pass
                handler = previous[signum]
                if handler is not None:
                    try:
                        signal.signal(signum, handler)
                    except (OSError, ValueError, TypeError):
                        pass

        return restore

    def execute(self) -> None:
        """Execute jobs from the DAG"""
        asyncio.run(self._execute())

    async def _execute(self) -> None:
        """Execute jobs from the DAG"""
        self.jobs_with_errors = []
        loop = asyncio.get_running_loop()
        stop_event = asyncio.Event()
        restore = self._install_signal_handlers(loop, stop_event)
        try:
            await self._drive(stop_event)
            if stop_event.is_set():
                await self._shutdown()
        finally:
            restore()

    @abstractmethod
    async def _drive(self, stop_event: asyncio.Event) -> None:
        """Submit/run jobs until the DAG drains or a stop is requested."""

    @abstractmethod
    async def _shutdown(self) -> None:
        """Tear down whatever is still running after a stop request."""


class LocalExecutor(AsyncExecutor):
    """Run jobs locally as async subprocesses."""

    def __init__(
        self,
        scheduler: BaseScheduler,
        *,
        install_signal_handlers: bool = False,
        shutdown_grace_period: float = 10.0,
        storage: Optional[StorageProvider] = None,
    ) -> None:
        super().__init__(
            scheduler,
            install_signal_handlers=install_signal_handlers,
        )
        self.shutdown_grace_period = shutdown_grace_period
        self.storage: StorageProvider = (
            storage if storage is not None else LocalStorageProvider()
        )
        self.running: List[
            Tuple[
                Job,
                Context,
                asyncio.Task[int],
                int,
            ]
        ] = []

    async def run_job(self, job: Job) -> None:
        """Run a job"""
        cmd = job.shell
        logger.info("Running: %s", cmd)
        context = Context()
        start_time = time.monotonic_ns()
        try:
            for path in job.inputs:
                await self.storage.stage_in(path)
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
        except OSError as e:
            # OSError covers the expected launch failures (command not found,
            # permission denied, pipe/FIFO errors). Anything else is a bug in
            # the caller's pipeline and is left to propagate.
            logger.error("failed to start command: %s", job.shell)
            logger.error("Error: %s", str(e))
            self.jobs_with_errors.append(job)
            self.start_new_jobs = False
            # A later pipeline stage can fail to spawn after earlier stages
            # are already running; terminate them so they do not outlive the
            # workflow, then release the remaining resources. A proc-sub inner
            # launch failure can also surface from cleanup() here (narrow
            # race): the job is already recorded, so swallow it.
            await self._terminate_context(context)
            await _cleanup_quietly(context)
            self._remove_outputs(job)

    async def _terminate_context(self, context: Context) -> None:
        """Tear down processes already spawned in an aborted context.

        Used when a job aborts mid-launch (a later pipeline stage failed to
        spawn): earlier stages are already running and would otherwise be
        leaked. SIGTERM them, allow a grace period, then SIGKILL any
        survivors, waiting for each to be reaped.
        """
        live = [
            subcommand.proc
            for subcommand in context.commands
            if subcommand.proc is not None
            and subcommand.proc.returncode is None
        ]
        if not live:
            return
        _signal_procs(context, signal.SIGTERM)
        waits = [asyncio.create_task(proc.wait()) for proc in live]
        await asyncio.wait(waits, timeout=self.shutdown_grace_period)
        _kill_survivors(context)
        await asyncio.wait(waits)

    def _remove_outputs(self, job: Job) -> None:
        """Delete a failed job's declared outputs.

        A partial or unverified output must not survive a failure: resume
        (:meth:`~sentieon_cli.dag.DAG.skip_satisfied`) treats an existing
        declared output as proof that its producer succeeded. Directories are
        deliberately not removed (never recurse into user paths); a warning
        is logged instead.
        """
        for path in job.outputs:
            try:
                path.unlink()
            except FileNotFoundError:
                continue
            except OSError as exc:
                logger.warning(
                    "could not remove output %s of failed job %s: %s",
                    path,
                    job,
                    exc,
                )
            else:
                logger.info("Removed output of failed job %s: %s", job, path)

    async def _drive(self, stop_event: asyncio.Event) -> None:
        """Schedule and run jobs until the DAG drains or we are stopped."""
        ready_jobs = self.scheduler.start()
        for job in ready_jobs:
            if not self.start_new_jobs:
                # A job in this batch failed to launch; drain, don't start
                # the rest -- matching the policy the main loop applies.
                break
            await self.run_job(job)

        stop_task = asyncio.ensure_future(stop_event.wait())
        try:
            while self.running and not stop_event.is_set():
                wait_on: List[Any] = [job[2] for job in self.running]
                wait_on.append(stop_task)
                done, _pending = await asyncio.wait(
                    wait_on,
                    return_when=asyncio.FIRST_COMPLETED,
                )

                if stop_event.is_set():
                    break

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
                        if not subcommand.proc:
                            logger.error(
                                "Sub-command has no process: %s",
                                subcommand,
                            )
                            cmd_failed = True
                            continue
                        ret = (
                            await subcommand.proc.wait()
                        )  # Wait on all sub-commands
                        # A subcommand killed by SIGPIPE is not a failure:
                        # it was writing to a pipe whose reader exited early
                        # (e.g. `... | head`), which is normal in default
                        # (non-pipefail) bash. Real crashes surface as other
                        # non-zero codes and still fail the job.
                        if ret == -signal.SIGPIPE:
                            logger.debug(
                                "Sub-command exited on SIGPIPE (ok): %s",
                                subcommand,
                            )
                        elif ret != 0 and not subcommand.fail_ok:
                            logger.error(
                                f"Return code {ret} for sub-command: "
                                f"{subcommand}"
                            )
                            cmd_failed = True
                            if ret == -9:
                                logger.error(
                                    "Sub-command received SIGKILL. "
                                    "Possible out-of-memory?"
                                )

                    if not cmd_failed:
                        for path in job.outputs:
                            try:
                                await self.storage.stage_out(path)
                            except OSError as exc:
                                logger.error(
                                    "missing output %s for %s: %s",
                                    path,
                                    job.shell,
                                    exc,
                                )
                                cmd_failed = True

                    # Release resources before recording the verdict. A
                    # proc-sub inner command that failed to launch surfaces
                    # its OSError from the background task here; treat it as
                    # a launch failure for this job rather than letting it
                    # abort the run while siblings are still mid-flight.
                    try:
                        await context.cleanup()
                    except OSError as exc:
                        logger.error("failed to start command: %s", job.shell)
                        logger.error("Error: %s", str(exc))
                        cmd_failed = True

                    if cmd_failed:
                        logger.error("Command failure: %s", job.shell)
                        self.jobs_with_errors.append(job)
                        self.start_new_jobs = False
                        self._remove_outputs(job)
                    else:
                        logger.info(
                            f"Finished command in "
                            f"{total_seconds:.2f}: {job.shell}"
                        )

                if not self.start_new_jobs:
                    # Don't start new jobs
                    continue

                ready_jobs = {
                    new_job
                    for completed_job in finished_jobs
                    for new_job in self.scheduler.job_finished(
                        completed_job[0]
                    )
                }

                # Run the ready jobs
                for job in ready_jobs:
                    if not self.start_new_jobs:
                        # A job in this batch failed to launch; drain the
                        # rest of the run without starting more.
                        break
                    await self.run_job(job)
        finally:
            stop_task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await stop_task

    async def _shutdown(self) -> None:
        """Terminate running jobs, escalate to SIGKILL, then clean up."""
        running = self.running
        self.running = []

        for _job, context, _task, _start in running:
            _signal_procs(context, signal.SIGTERM)

        tasks = [tup[2] for tup in running]
        if tasks:
            await asyncio.wait(tasks, timeout=self.shutdown_grace_period)

        for _job, context, _task, _start in running:
            _kill_survivors(context)
        if tasks:
            await asyncio.wait(tasks)

        for _job, context, _task, _start in running:
            # A proc-sub inner launch failure can surface from cleanup() as we
            # tear down after an interrupt; log it rather than let it escape,
            # and keep releasing the remaining contexts.
            await _cleanup_quietly(context)

        # Interrupted jobs never had their outputs verified; a partial file
        # must not survive to fool a later resume. Jobs that completed in the
        # same wake-up as the stop signal are conservatively included -- they
        # were never verified either.
        for job, _context, _task, _start in running:
            self._remove_outputs(job)
