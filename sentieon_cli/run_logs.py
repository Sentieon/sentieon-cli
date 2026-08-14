"""
The log directory for a single run
"""

from __future__ import annotations

import dataclasses
import datetime
import logging
import pathlib
import re
import shlex
import shutil
import sys
from typing import IO, List, Optional, TYPE_CHECKING

from .logging import add_file_handler, remove_file_handler
from .util import __version__

if TYPE_CHECKING:
    from .job import Job
    from .shell_pipeline import Command

_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")

_STDERR = "stderr"
_STDOUT = "stdout"


def sanitize(component: str) -> str:
    """Restrict a path component to filesystem-safe characters"""
    return _UNSAFE.sub("-", component)


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
        self.file_handler: Optional[logging.FileHandler] = None

    def setup(self) -> None:
        """Prepare the log directory and start writing `run.log`"""
        self.create_dirs()
        self.write_command()
        self.file_handler = add_file_handler(self.run_log)

    def create_dirs(self) -> None:
        """Create the log directory, clearing logs from any previous run"""
        # The log directory itself is never removed - the user may point
        # `--log_dir` at an existing directory holding unrelated files.
        self.log_dir.mkdir(parents=True, exist_ok=True)
        if self.task_logs.exists():
            shutil.rmtree(self.task_logs)
        self.task_logs.mkdir(parents=True)

    def job_sink(self, job: Job) -> JobLogSink:
        """Create the log sink for a job"""
        return JobLogSink(self.task_logs, job)

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
