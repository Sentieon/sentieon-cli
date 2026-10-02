"""Execute jobs"""

import asyncio
import concurrent.futures
import contextlib
import os
import pathlib
import resource
import signal
import sys
import threading
import time
from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, List, Optional, Tuple

from .job import Job
from .logging import get_logger
from .run_logs import JobLogSink, RunLogs
from .scheduler import BaseScheduler
from .shell_pipeline import Command, Context

logger = get_logger(__name__)

# Lines of a failing process's log to quote in the failure report
TAIL_LINES = 20
# Only the tail of a log is read; tool logs can be arbitrarily large
TAIL_BYTES = 64 * 1024


def _log_tail(path: pathlib.Path, lines: int = TAIL_LINES) -> List[str]:
    """Return the last ``lines`` lines of a log file.

    Reads at most ``TAIL_BYTES`` from the end, so quoting the tail of a
    multi-gigabyte log is cheap. A line split by the seek is dropped with the
    rest of the truncated prefix.
    """
    try:
        with open(path, "rb") as handle:
            size = handle.seek(0, os.SEEK_END)
            handle.seek(max(0, size - TAIL_BYTES))
            data = handle.read()
    except OSError as exc:
        logger.debug("could not read %s: %s", path, exc)
        return []
    return data.decode("utf-8", errors="replace").splitlines()[-lines:]


def _maxrss_bytes(ru: resource.struct_rusage) -> int:
    """Return a child's peak resident set size in bytes.

    ``ru_maxrss`` is already bytes on macOS; every other platform we run on
    reports kibibytes.
    """
    if sys.platform == "darwin":
        return int(ru.ru_maxrss)
    return int(ru.ru_maxrss) * 1024


def _mib(nbytes: int) -> str:
    """Render a byte count as MiB for the human-readable log lines."""
    return f"{nbytes / (1024 * 1024):.1f} MiB"


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

    ``thread_pool_size`` sizes the loop's default thread pool for the run;
    ``None`` leaves the loop's own default alone, which is what an executor
    that spawns no processes wants.
    """

    def __init__(
        self,
        scheduler: BaseScheduler,
        *,
        install_signal_handlers: bool = False,
        thread_pool_size: Optional[int] = None,
    ) -> None:
        super().__init__(scheduler)
        self.install_signal_handlers = install_signal_handlers
        self.thread_pool_size = thread_pool_size
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
        if self.thread_pool_size is not None:
            # Threads are created lazily, so a generous cap costs nothing
            # unless it is used, and asyncio.run() shuts the executor down on
            # exit.
            loop.set_default_executor(
                concurrent.futures.ThreadPoolExecutor(
                    max_workers=self.thread_pool_size,
                    thread_name_prefix="sentieon-wait",
                )
            )
        stop_event = asyncio.Event()
        restore = self._install_signal_handlers(loop, stop_event)
        try:
            try:
                await self._drive(stop_event)
            except BaseException:
                # Tear the children down before the exception unwinds into
                # asyncio.run's Runner.close, which joins the wait threads --
                # and those threads are blocked in wait4 on live children.
                # BaseException so cancellation and KeyboardInterrupt are
                # covered too.
                self.start_new_jobs = False
                try:
                    await self._shutdown()
                except Exception as exc:
                    # A failed teardown must not mask the original error.
                    logger.error("teardown after failure also failed: %s", exc)
                raise
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
    """Run jobs locally as sub-processes, driven by an asyncio loop.

    ``thread_pool_size`` defaults high because the loop's default pool serves
    two kinds of blocking work at once: one thread per in-flight ``wait()``
    and the proc-sub FIFO opens. A FIFO open queued behind saturated wait
    threads would deadlock the run, so the pool must comfortably exceed the
    number of processes a run keeps in flight.
    """

    def __init__(
        self,
        scheduler: BaseScheduler,
        *,
        install_signal_handlers: bool = False,
        shutdown_grace_period: float = 10.0,
        run_logs: Optional[RunLogs] = None,
        thread_pool_size: int = 512,
    ) -> None:
        super().__init__(
            scheduler,
            install_signal_handlers=install_signal_handlers,
            thread_pool_size=thread_pool_size,
        )
        self.shutdown_grace_period = shutdown_grace_period
        self.run_logs = run_logs
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
        sink = self.run_logs.job_sink(job) if self.run_logs else None
        context = Context() if sink is None else Context(log_sink=sink)
        start_time = time.monotonic_ns()
        try:
            if sink is None:
                # No log directory: the processes inherit our stderr
                proc = await cmd.run(
                    context,
                    stderr=sys.stderr,
                )
            else:
                # Unset streams are resolved against the sink
                proc = await cmd.run(context)
            self.running.append(
                (
                    job,
                    context,
                    asyncio.create_task(proc.async_wait()),
                    start_time,
                )
            )
        except OSError as e:
            # OSError covers the expected launch failures (command not found,
            # permission denied, pipe/FIFO errors). Anything else is a bug in
            # the caller's pipeline and is left to propagate.
            logger.error("failed to start command: %s", job.shell)
            logger.error("Error: %s", str(e))
            # Stages that did spawn have logs; keep them and point at them.
            self._report_logs(job, sink)
            self.jobs_with_errors.append(job)
            self.start_new_jobs = False
            # A later pipeline stage can fail to spawn after earlier stages
            # are already running; terminate them so they do not outlive the
            # workflow, then release the remaining resources. A proc-sub inner
            # launch failure can also surface from cleanup() here (narrow
            # race): the job is already recorded, so swallow it.
            await self._terminate_context(context)
            await _cleanup_quietly(context)

    def _report_logs(self, job: Job, sink: Optional[JobLogSink]) -> None:
        """Point at the logs of a job that was aborted rather than reaped."""
        if sink is None:
            return
        paths = sink.log_paths()
        if paths:
            logger.error(
                "Logs from %s: %s",
                job,
                ", ".join(str(path) for path in paths),
            )

    def _report_failure(
        self,
        job: Job,
        subcommand: Command,
        ret: int,
        sink: Optional[JobLogSink],
    ) -> None:
        """Report a failed sub-command with its pid, log path and log tail.

        Without a log -- no log directory, or a process that never spawned --
        the caller's messages are the whole report.
        """
        log_path = sink.stderr_log_for(subcommand) if sink else None
        if log_path is None:
            return
        pid = subcommand.proc.pid if subcommand.proc else None
        report = [
            f"Failed sub-command of {job}: {subcommand}",
            f"  exit code: {ret}, pid: {pid}",
            f"  log: {log_path}",
        ]
        tail = _log_tail(log_path)
        if tail:
            report.append(f"  last {len(tail)} line(s) of the log:")
            report.extend(f"    {line}" for line in tail)
        else:
            report.append("  the log is empty")
        logger.error("\n".join(report))

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
        waits = [asyncio.create_task(proc.async_wait()) for proc in live]
        await asyncio.wait(waits, timeout=self.shutdown_grace_period)
        _kill_survivors(context)
        await asyncio.wait(waits)

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
                    failures: List[Tuple[Command, int]] = []
                    # Kernel resource usage, aggregated over the job's
                    # processes. RSS is the per-process peak rather than a
                    # sum: concurrent pipeline stages peak at different
                    # times, so a sum would overstate the job.
                    have_rusage = False
                    job_utime = 0.0
                    job_stime = 0.0
                    job_maxrss = 0
                    n_procs = 0
                    sink = context.log_sink
                    for subcommand in context.commands:
                        if not subcommand.proc:
                            logger.error(
                                "Sub-command has no process: %s",
                                subcommand,
                            )
                            cmd_failed = True
                            continue
                        n_procs += 1
                        ret = (
                            await subcommand.proc.async_wait()
                        )  # Wait on all sub-commands
                        # Reported for every reaped process, whatever the
                        # job's outcome; signal-killed children have rusage
                        # too.
                        rusage = subcommand.proc.rusage
                        rss = (
                            _maxrss_bytes(rusage)
                            if rusage is not None
                            else None
                        )
                        if rusage is not None and rss is not None:
                            have_rusage = True
                            job_utime += rusage.ru_utime
                            job_stime += rusage.ru_stime
                            job_maxrss = max(job_maxrss, rss)
                            logger.debug(
                                "rusage for %s [pid %d, %s]: utime=%.2fs "
                                "stime=%.2fs maxrss=%.1fMiB minflt=%d "
                                "majflt=%d inblock=%d oublock=%d "
                                "nvcsw=%d nivcsw=%d",
                                job,
                                subcommand.proc.pid,
                                subcommand,
                                rusage.ru_utime,
                                rusage.ru_stime,
                                rss / (1024 * 1024),
                                rusage.ru_minflt,
                                rusage.ru_majflt,
                                rusage.ru_inblock,
                                rusage.ru_oublock,
                                rusage.ru_nvcsw,
                                rusage.ru_nivcsw,
                            )
                        if self.run_logs is not None:
                            self.run_logs.record_process_metrics(
                                job=job,
                                command=subcommand,
                                pid=subcommand.proc.pid,
                                exit_code=ret,
                                stage=(
                                    sink.stage_index_for(subcommand)
                                    if sink is not None
                                    else None
                                ),
                                rusage=rusage,
                                maxrss_bytes=rss,
                            )
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
                            failures.append((subcommand, ret))
                            if ret == -9:
                                logger.error(
                                    "Sub-command received SIGKILL. "
                                    "Possible out-of-memory?"
                                )

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
                        # The sink's files are closed by cleanup(), so their
                        # tails can be read now.
                        for failed, code in failures:
                            self._report_failure(job, failed, code, sink)
                        if self.run_logs is not None:
                            logger.error(
                                "Task logs for debugging: %s",
                                self.run_logs.task_logs,
                            )
                        self.jobs_with_errors.append(job)
                        self.start_new_jobs = False
                    elif have_rusage:
                        logger.info(
                            "Finished command in %.2fs (user %.1fs, "
                            "sys %.1fs, max proc RSS %s): %s",
                            total_seconds,
                            job_utime,
                            job_stime,
                            _mib(job_maxrss),
                            job.shell,
                        )
                    else:
                        logger.info(
                            f"Finished command in "
                            f"{total_seconds:.2f}: {job.shell}"
                        )
                    if self.run_logs is not None:
                        self.run_logs.record_job_metrics(
                            job=job,
                            failed=cmd_failed,
                            wall_s=total_seconds,
                            user_s=job_utime if have_rusage else None,
                            sys_s=job_stime if have_rusage else None,
                            max_proc_rss_bytes=(
                                job_maxrss if have_rusage else None
                            ),
                            processes=n_procs,
                        )
                    if sink is not None:
                        sink.finalize(success=not cmd_failed)

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

        for job, context, _task, _start in running:
            # A proc-sub inner launch failure can surface from cleanup() as we
            # tear down after an interrupt; log it rather than let it escape,
            # and keep releasing the remaining contexts.
            await _cleanup_quietly(context)
            # An interrupted job keeps all of its logs; they explain how far
            # it got.
            self._report_logs(job, context.log_sink)
