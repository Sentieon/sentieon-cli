"""Complex shell pipelines with process substitution.

POSIX only: pipelines are wired together with OS pipes and FIFOs, and process
substitution relies on named pipes, so this module does not work on Windows.
"""

from __future__ import annotations
import asyncio
import dataclasses
import fcntl
import os
import pathlib
import tempfile
import shlex
from abc import ABC, abstractmethod
from typing import Any, Dict, IO, List, Optional, TYPE_CHECKING, Union

from .logging import get_logger

if TYPE_CHECKING:
    from .run_logs import JobLogSink

logger = get_logger(__name__)


def _freeze(value: Any) -> Any:
    """Return an order-independent, hashable view of ``value``."""
    if isinstance(value, dict):
        return frozenset((k, _freeze(v)) for k, v in value.items())
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(v) for v in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(v) for v in value)
    return value


def _set_pipe_size(fd: int, size: int) -> None:
    """Best-effort enlarge the OS pipe buffer backing ``fd``.

    Lets a throughput-bound pipeline hold more in-flight data between
    stages.
    """
    set_pipe_sz = getattr(fcntl, "F_SETPIPE_SZ", None)
    if set_pipe_sz is None:
        return
    try:
        fcntl.fcntl(fd, set_pipe_sz, size)
    except OSError as exc:
        logger.warning(
            "could not set pipe size to %d (%s); using the default. "
            "Raise fs.pipe-max-size to allow larger buffers.",
            size,
            exc,
        )


@dataclasses.dataclass
class _ProcSubWait:
    """A proc-sub FIFO whose blocking open may need unblocking."""

    fifo_path: str
    opened: bool = False  # run_inner's open() returned (set pre-await)
    unblocked: bool = False  # cleanup already opened both ends


class Context:
    """Holds shared resources."""

    def __init__(self, log_sink: Optional[JobLogSink] = None) -> None:
        # This directory is already unique per pipeline run
        self.temp_dir = tempfile.TemporaryDirectory()
        # Background tasks to wait for (e.g., proc sub processes)
        self.tasks: List[asyncio.Task[Any]] = []
        # Hold the commands run in this context
        self.commands: List[Command] = []
        self._counter = 0  # Monotonic counter
        self.file_handles: List[IO[Any]] = []
        # Per-process log files; without a sink, stdout/stderr are inherited
        self.log_sink = log_sink
        # Set by cleanup(); a proc sub whose FIFO open was unblocked by
        # cleanup must not spawn its inner command.
        self.closing: bool = False
        self._procsub_waits: List[_ProcSubWait] = []

    def get_new_fifo(self) -> str:
        """Creates a new unique FIFO and returns its path."""
        self._counter += 1
        # Simple, readable names like fifo_1, fifo_2, etc.
        path = os.path.join(self.temp_dir.name, f"fifo_{self._counter}")
        os.mkfifo(path)
        return path

    async def cleanup(self) -> None:
        """Wait for background tasks, then release all resources.

        A proc sub's blocking FIFO open() never returns if the outer
        command exited without opening the other end, so cleanup first
        unblocks those opens, then waits for the tasks. Task exceptions
        are re-raised only after every resource has been released.
        """
        self.closing = True
        unblock_fds: List[int] = []
        try:
            while True:
                for wait in self._procsub_waits:
                    if wait.opened or wait.unblocked:
                        continue
                    wait.unblocked = True
                    # Hold both ends open: this unblocks a blocked open()
                    # of either end, and read-then-write cannot raise
                    # ENXIO.
                    unblock_fds.append(
                        os.open(wait.fifo_path, os.O_RDONLY | os.O_NONBLOCK)
                    )
                    unblock_fds.append(
                        os.open(wait.fifo_path, os.O_WRONLY | os.O_NONBLOCK)
                    )
                # Tasks can register new proc subs (nested pipelines), so
                # re-snapshot until everything is done.
                pending = [t for t in self.tasks if not t.done()]
                if not pending:
                    break
                await asyncio.gather(*pending, return_exceptions=True)
        finally:
            for fd in unblock_fds:
                os.close(fd)
        for fh in self.file_handles:
            fh.close()
        self.temp_dir.cleanup()
        for task in self.tasks:
            if task.cancelled():
                continue
            exc = task.exception()
            if exc is not None and not isinstance(exc, asyncio.CancelledError):
                raise exc


class ShellNode(ABC):
    """Abstract base class for any node in the shell syntax tree."""


class Command(ShellNode):
    """Represents a single command (e.g., 'grep', 'cat')."""

    def __init__(
        self,
        executable: str,
        *args: str | ProcSub,
        fail_ok: bool = False,
        exec_kwargs: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.executable = executable
        self.args = list(args)
        self.fail_ok = fail_ok
        self.proc: Union[asyncio.subprocess.Process, None] = None
        self.exec_kwargs: Dict[str, Any] = exec_kwargs if exec_kwargs else {}

    async def run(
        self,
        context: Context,
        stdin: Union[IO[Any], int, None] = None,
        stdout: Union[IO[Any], int, None] = None,
        stderr: Union[IO[Any], int, None] = None,
    ) -> asyncio.subprocess.Process:
        # 1. Resolve Arguments (Handle Process Substitutions)
        final_args = [self.executable]

        for arg in self.args:
            if isinstance(arg, ProcSub):
                # Recursively setup the process substitution
                fifo_path = await arg.setup(context)
                final_args.append(fifo_path)
            else:
                final_args.append(str(arg))

        # 2. Capture the streams the caller left unset
        # An explicit stream (a pipe between stages, a file_output, a proc-sub
        # FIFO) always wins; only an inherited stream goes to the sink.
        if context.log_sink is not None:
            if stderr is None:
                stderr_log = context.log_sink.open_stderr(self)
                context.file_handles.append(stderr_log)
                stderr = stderr_log
            if stdout is None:
                stdout_log = context.log_sink.open_stdout(self)
                context.file_handles.append(stdout_log)
                stdout = stdout_log

        # 3. Run the process
        # We allow the caller to wait for this process
        self.proc = await asyncio.create_subprocess_exec(
            *final_args,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            **self.exec_kwargs,
        )
        context.commands.append(self)
        return self.proc

    def __str__(self) -> str:
        rendered = shlex.join([self.executable] + [str(x) for x in self.args])
        # Render env/cwd so the command is self-contained bash (backends
        # submit str(job.shell)). Other exec_kwargs (close_fds, ...) have no
        # shell rendering and are omitted. env is a conventional K=V prefix
        # that merges into the inherited environment.
        env = self.exec_kwargs.get("env")
        if env:
            assignments = " ".join(
                f"{key}={shlex.quote(str(val))}" for key, val in env.items()
            )
            rendered = f"{assignments} {rendered}"
        cwd = self.exec_kwargs.get("cwd")
        if cwd:
            rendered = f"(cd {shlex.quote(str(cwd))} && {rendered})"
        return rendered

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({self.executable}, "
            + ", ".join([repr(x) for x in self.args])
            + f", fail_ok={self.fail_ok}, "
            + f"exec_kwargs={self.exec_kwargs!r})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Command):
            return NotImplemented
        return (
            self.executable == other.executable
            and self.args == other.args
            and self.fail_ok == other.fail_ok
            and self.exec_kwargs == other.exec_kwargs
        )

    def __hash__(self) -> int:
        return hash(
            tuple([self.executable, self.fail_ok] + self.args)
            + (_freeze(self.exec_kwargs),)
        )


class Pipeline(ShellNode):
    """Represents a sequence of commands connected by pipes (|)."""

    def __init__(
        self,
        *nodes: Command,
        pipe_size: Optional[int] = None,
        file_input: Optional[pathlib.Path] = None,
        file_output: Optional[pathlib.Path] = None,
    ):
        self.nodes = list(nodes)
        if not self.nodes:
            raise ValueError("Pipeline nodes cannot be empty")
        # Best-effort OS pipe-buffer size (bytes) for every internal pipe;
        # None keeps the kernel default. A local-execution tuning hint, so
        # it is not part of pipeline identity (like Job.threads).
        self.pipe_size = pipe_size
        self.file_input = file_input
        self.file_output = file_output

    async def run(
        self,
        context: Context,
        stdin: Union[IO[Any], int, None] = None,
        stdout: Union[IO[Any], int, None] = None,
        stderr: Union[IO[Any], int, None] = None,
    ) -> asyncio.subprocess.Process:
        # We cannot have handles from both files and to run().
        if self.file_input is not None and stdin is not None:
            raise ValueError(
                f"Pipeline {self} was given both a file_input and a stdin "
                f"stream ({stdin!r}); provide only one."
            )
        if self.file_output is not None and stdout is not None:
            raise ValueError(
                f"Pipeline {self} was given both a file_output and a stdout "
                f"stream ({stdout!r}); provide only one."
            )

        processes = []
        # We need to close these FDs in the parent after passing them to
        # children
        fds_to_close = []

        # pass the file handle as input to the pipeline
        current_stdin = stdin
        if self.file_input is not None:
            current_stdin = open(self.file_input, "rb")
            context.file_handles.append(current_stdin)

        try:
            # Chain all commands except the last one
            for i in range(len(self.nodes) - 1):
                # Create the OS pipe linking this command to the next.
                read_fd, write_fd = os.pipe()
                fds_to_close.extend([read_fd, write_fd])
                if self.pipe_size is not None:
                    _set_pipe_size(write_fd, self.pipe_size)

                # Run the node (Command or sub-Pipeline)
                # stdout goes to the write-end of the pipe
                proc = await self.nodes[i].run(
                    context,
                    stdin=current_stdin,
                    stdout=write_fd,
                    stderr=stderr,
                )
                processes.append(proc)

                # Close write_fd in parent (child has it now)
                os.close(write_fd)
                fds_to_close.remove(write_fd)

                # Next command reads from the read-end
                current_stdin = read_fd

            # Send the pipeline output to a file
            if self.file_output is not None:
                stdout = open(self.file_output, "wb")
                context.file_handles.append(stdout)

            # Run the final command
            last_proc = await self.nodes[-1].run(
                context,
                stdin=current_stdin,
                stdout=stdout,
                stderr=stderr,
            )
            processes.append(last_proc)

        finally:
            # Cleanup FDs
            for fd in fds_to_close:
                try:
                    os.close(fd)
                except OSError:
                    pass

        return processes[-1]  # Return the last process handle

    def __str__(self) -> str:
        res = []
        if self.file_input:
            res.append(f"<'{self.file_input}' ")
        res.append(" | ".join(str(node) for node in self.nodes))
        if self.file_output:
            res.append(f" >'{self.file_output}'")
        return "".join(res)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            + ", ".join(repr(x) for x in self.nodes)
            + f", pipe_size={self.pipe_size}, "
            + f"file_input={self.file_input}, "
            + f"file_output={self.file_output})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Pipeline):
            return NotImplemented
        # pipe_size is a tuning hint, not part of identity (see __init__).
        return (
            self.nodes == other.nodes
            and self.file_input == other.file_input
            and self.file_output == other.file_output
        )

    def __hash__(self) -> int:
        attrs: List[Any] = self.nodes + [self.file_input, self.file_output]
        return hash(tuple(attrs))


class ProcSub(ABC):
    """Abstract base class for Process Substitutions (<(...) and >(...))."""

    def __init__(self, node: Pipeline) -> None:
        self.node = node
        self.fifo_path: str | None = None

    @abstractmethod
    async def setup(self, context: Context) -> str:
        """Prepares the FIFO and starts the background process."""
        pass

    def __repr__(self) -> str:
        return self.__class__.__name__ + "(" + repr(self.node) + ")"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ProcSub):
            return NotImplemented
        return type(self) is type(other) and self.node == other.node

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.node))


class InputProcSub(ProcSub):
    """
    Represents <(cmd).
    The inner command WRITES to a FIFO.
    The outer command receives the filename.
    """

    async def setup(self, context: Context) -> str:
        fifo_path = context.get_new_fifo()
        self.fifo_path = fifo_path
        wait = _ProcSubWait(fifo_path)

        # We must open the FIFO for writing.
        # WARNING: opening a FIFO blocks until the other end is opened.
        # We run the open() in a thread to prevent blocking the async event
        # loop.
        def open_fifo_write() -> IO[Any]:
            return open(fifo_path, "w")

        # Start the task that waits for the open, then runs the process
        async def run_inner() -> None:
            # This awaits until the outer command opens the file for
            # reading (or until cleanup() unblocks the open)!
            write_handle = await asyncio.to_thread(open_fifo_write)
            wait.opened = True
            try:
                if context.closing:
                    return  # cleanup unblocked us; do not spawn
                proc = await self.node.run(
                    context,
                    stdout=write_handle,
                )
                await proc.wait()
            finally:
                write_handle.close()

        # Schedule this task in the background
        context.tasks.append(asyncio.create_task(run_inner()))
        context._procsub_waits.append(wait)

        return fifo_path

    def __str__(self) -> str:
        return "<(" + str(self.node) + ")"


class OutputProcSub(ProcSub):
    """
    Represents >(cmd).
    The inner command READS from a FIFO.
    The outer command receives the filename.
    """

    async def setup(self, context: Context) -> str:
        fifo_path = context.get_new_fifo()
        self.fifo_path = fifo_path
        wait = _ProcSubWait(fifo_path)

        def open_fifo_read() -> IO[Any]:
            return open(fifo_path, "r")

        async def run_inner() -> None:
            # This awaits until the outer command opens the file for
            # writing (or until cleanup() unblocks the open)!
            read_handle = await asyncio.to_thread(open_fifo_read)
            wait.opened = True
            try:
                if context.closing:
                    return  # cleanup unblocked us; do not spawn
                proc = await self.node.run(context, stdin=read_handle)
                await proc.wait()
            finally:
                read_handle.close()

        context.tasks.append(asyncio.create_task(run_inner()))
        context._procsub_waits.append(wait)

        return fifo_path

    def __str__(self) -> str:
        return ">(" + str(self.node) + ")"
