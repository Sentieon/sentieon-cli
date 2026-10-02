"""
The log directory for a single run
"""

from __future__ import annotations

import dataclasses
import datetime
import logging
import pathlib
import shlex
import shutil
import sys
from typing import IO, Iterable, List, Optional, TYPE_CHECKING

from .logging import add_file_handler, get_logger, remove_file_handler
from .util import __version__, sanitize

if TYPE_CHECKING:
    import resource

    from .job import Job
    from .shell_pipeline import Command

logger = get_logger(__name__)

_STDERR = "stderr"
_STDOUT = "stdout"

# The rusage-derived columns of a process row, empty as a group when a
# process was reaped without rusage.
_RUSAGE_COLUMNS = [
    "user_s",
    "sys_s",
    "maxrss_mib",
    "minflt",
    "majflt",
    "inblock",
    "oublock",
    "nvcsw",
    "nivcsw",
]

PROCESS_METRICS_COLUMNS = [
    "job_id",
    "task_name",
    "pid",
    "stage",
    "exit_code",
    *_RUSAGE_COLUMNS,
    "command",
]

JOB_METRICS_COLUMNS = [
    "job_id",
    "task_name",
    "threads",
    "status",
    "wall_s",
    "user_s",
    "sys_s",
    "max_proc_rss_mib",
    "processes",
]

_MIB = 1024 * 1024


def _tsv_field(value: object) -> str:
    """Render a value as a TSV field that cannot break the row layout"""
    text = "" if value is None else str(value)
    for char in ("\t", "\n", "\r"):
        text = text.replace(char, " ")
    return text


def _write_row(path: pathlib.Path, fields: Iterable[object]) -> None:
    """Append one row to a metrics file.

    A metrics file is a by-product of the run: a write failure (full disk,
    read-only log directory) is reported once and otherwise ignored.
    """
    row = "\t".join(_tsv_field(field) for field in fields)
    try:
        with open(path, "a") as handle:
            handle.write(row + "\n")
    except OSError as exc:
        logger.warning("could not write %s: %s", path, exc)


@dataclasses.dataclass
class _ProcessLog:
    """A log file opened for one spawned process"""

    command: Command
    stage_index: int
    path: pathlib.Path
    kind: str
    header_bytes: int


class JobLogSink:
    """The per-process log files of a single job.

    Attached to the job's ``Context``, so every process spawned for the job --
    including process substitutions -- writes its stderr (and its stdout, when
    that would otherwise be inherited) to its own file.
    """

    def __init__(self, task_logs: pathlib.Path, job: Job) -> None:
        self.job_id = job.job_id
        self.name = job.name
        self.task_name = job.task_name
        self.log_dir = task_logs / sanitize(job.task_name)
        self._prefix = sanitize(job.job_id)
        self._stage_index = 0  # Next unallocated stage index
        self._logs: List[_ProcessLog] = []

    def open_stderr(self, command: Command) -> IO[bytes]:
        """Open the stderr log of the next process to spawn"""
        return self._open(command, _STDERR)

    def open_stdout(self, command: Command) -> IO[bytes]:
        """Open the stdout log of the next process to spawn"""
        return self._open(command, _STDOUT)

    def _stage_for(self, command: Command) -> int:
        """The stage index of a command, allocated on first use.

        Both streams of one process share its index, so an already-seen
        command reuses it. Matched by identity: equal-but-distinct commands
        are distinct processes.
        """
        for log in self._logs:
            if log.command is command:
                return log.stage_index
        stage_index = self._stage_index
        self._stage_index += 1
        return stage_index

    def stage_index_for(self, command: Command) -> Optional[int]:
        """The stage index already allocated to a command, if any.

        Unlike ``_stage_for`` this never allocates: a command whose streams
        were all explicit (so no log was opened) has no stage index.
        """
        for log in self._logs:
            if log.command is command:
                return log.stage_index
        return None

    def _open(self, command: Command, kind: str) -> IO[bytes]:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        stage_index = self._stage_for(command)
        suffix = ".log" if kind == _STDERR else ".stdout.log"
        path = self.log_dir / f"{self._prefix}.{stage_index}{suffix}"
        handle = open(path, "wb")
        header_bytes = 0
        if kind == _STDERR:
            # The child inherits the file offset, so the header must be on
            # disk before the FD is handed over or it would interleave with
            # (or land after) the process's own output. Stdout logs get no
            # header: they hold exactly what the tool would have printed.
            header = self._header(command, stage_index)
            handle.write(header)
            handle.flush()
            header_bytes = len(header)
        self._logs.append(
            _ProcessLog(command, stage_index, path, kind, header_bytes)
        )
        return handle

    def _header(self, command: Command, stage_index: int) -> bytes:
        timestamp = datetime.datetime.now().astimezone()
        return (
            f"# timestamp: {timestamp.isoformat(timespec='seconds')}\n"
            f"# task: {self.task_name}\n"
            f"# job: {self.name} ({self.job_id})\n"
            f"# stage: {stage_index}\n"
            f"# command: {command}\n"
        ).encode()

    def stderr_log_for(self, command: Command) -> Optional[pathlib.Path]:
        """The stderr log of a command, matched by identity.

        Equal-but-distinct ``Command`` objects can share a context, so the
        lookup must not fall back to ``__eq__``.
        """
        for log in self._logs:
            if log.kind == _STDERR and log.command is command:
                return log.path
        return None

    def log_paths(self) -> List[pathlib.Path]:
        """Every log file opened for this job"""
        return [log.path for log in self._logs]

    def finalize(self, success: bool) -> None:
        """Drop the logs of a successful job that hold nothing but a header.

        Empty logs are common (cleanup jobs, unused stdout) and only add
        clutter. After a failure everything is kept -- an empty log is itself
        informative.
        """
        if not success:
            return
        for log in self._logs:
            try:
                if log.path.stat().st_size <= log.header_bytes:
                    log.path.unlink()
            except OSError:
                pass


class RunLogs:
    """The log directory of a single sentieon-cli invocation"""

    def __init__(self, log_dir: pathlib.Path) -> None:
        self.log_dir = log_dir
        self.run_log = log_dir / "run.log"
        self.command_txt = log_dir / "command.txt"
        self.task_logs = log_dir / "task_logs"
        self.process_metrics = log_dir / "process_metrics.txt"
        self.job_metrics = log_dir / "job_metrics.txt"
        self.file_handler: Optional[logging.FileHandler] = None

    def setup(self) -> None:
        """Prepare the log directory and start writing `run.log`"""
        self.create_dirs()
        self.write_command()
        self.file_handler = add_file_handler(self.run_log)

    def create_dirs(self) -> None:
        """Create the log directory, clearing logs from any previous run.

        The metrics files are (re-)written with their header row here, so a
        rerun into the same directory replaces the previous run's rows and a
        run with no finished jobs still leaves a readable, empty table.
        """
        # The log directory itself is never removed - the user may point
        # `--log_dir` at an existing directory holding unrelated files.
        self.log_dir.mkdir(parents=True, exist_ok=True)
        if self.task_logs.exists():
            shutil.rmtree(self.task_logs)
        self.task_logs.mkdir(parents=True)
        self.process_metrics.write_text(
            "\t".join(PROCESS_METRICS_COLUMNS) + "\n"
        )
        self.job_metrics.write_text("\t".join(JOB_METRICS_COLUMNS) + "\n")

    def job_sink(self, job: Job) -> JobLogSink:
        """Create the log sink for a job"""
        return JobLogSink(self.task_logs, job)

    def record_process_metrics(
        self,
        *,
        job: Job,
        command: Command,
        pid: int,
        exit_code: int,
        stage: Optional[int],
        rusage: Optional["resource.struct_rusage"],
        maxrss_bytes: Optional[int],
    ) -> None:
        """Append one process to `process_metrics.txt`.

        Written for every reaped process, whatever the job's outcome; the
        rusage columns are empty for a process reaped without one.
        """
        if rusage is None:
            usage: List[object] = [None] * len(_RUSAGE_COLUMNS)
        else:
            rss = "" if maxrss_bytes is None else f"{maxrss_bytes / _MIB:.1f}"
            usage = [
                f"{rusage.ru_utime:.2f}",
                f"{rusage.ru_stime:.2f}",
                rss,
                rusage.ru_minflt,
                rusage.ru_majflt,
                rusage.ru_inblock,
                rusage.ru_oublock,
                rusage.ru_nvcsw,
                rusage.ru_nivcsw,
            ]
        _write_row(
            self.process_metrics,
            [job.job_id, job.task_name, pid, stage, exit_code]
            + usage
            + [command],
        )

    def record_job_metrics(
        self,
        *,
        job: Job,
        failed: bool,
        wall_s: float,
        user_s: Optional[float],
        sys_s: Optional[float],
        max_proc_rss_bytes: Optional[int],
        processes: int,
    ) -> None:
        """Append one job to `job_metrics.txt`.

        Only jobs that reach accounting get a row: a job that failed to
        launch, or one cut short by an interrupt, is absent.
        """
        rss = (
            ""
            if max_proc_rss_bytes is None
            else f"{max_proc_rss_bytes / _MIB:.1f}"
        )
        _write_row(
            self.job_metrics,
            [
                job.job_id,
                job.task_name,
                job.threads,
                "failed" if failed else "ok",
                f"{wall_s:.2f}",
                "" if user_s is None else f"{user_s:.2f}",
                "" if sys_s is None else f"{sys_s:.2f}",
                rss,
                processes,
            ],
        )

    def write_command(self) -> None:
        """Record the invocation so the run can be reproduced"""
        timestamp = datetime.datetime.now().astimezone()
        self.command_txt.write_text(
            f"command: {shlex.join(sys.argv)}\n"
            f"version: {__version__}\n"
            f"directory: {pathlib.Path.cwd()}\n"
            f"timestamp: {timestamp.isoformat(timespec='seconds')}\n"
        )

    def close(self) -> None:
        """Stop writing `run.log`"""
        if self.file_handler is not None:
            remove_file_handler(self.file_handler)
            self.file_handler = None
