"""Complex shell pipelines with process substitution"""

from __future__ import annotations
import asyncio
import io
import os
import pathlib
import tempfile
import shlex
import sys
from abc import ABC, abstractmethod
from typing import Any, Dict, IO, Iterable, List, Optional, Union

from .logging import get_logger

logger = get_logger(__name__)


class Context:
    """Holds shared resources."""

    def __init__(self) -> None:
        # This directory is already unique per pipeline run
        self.temp_dir = tempfile.TemporaryDirectory()
        # Background tasks to wait for (e.g., proc sub processes)
        self.tasks: List[asyncio.Task] = []
        # Hold the commands run in this context
        self.commands: List[Command] = []
        self._counter = 0  # Monotonic counter
        self.file_handles: List[io.IOBase] = []

    def get_new_fifo(self) -> str:
        """Creates a new unique FIFO and returns its path."""
        self._counter += 1
        # Simple, readable names like fifo_1, fifo_2, etc.
        path = os.path.join(self.temp_dir.name, f"fifo_{self._counter}")
        os.mkfifo(path)
        return path

    async def cleanup(self) -> None:
        if self.tasks:
            await asyncio.gather(*self.tasks)
        for fh in self.file_handles:
            fh.close()
        self.temp_dir.cleanup()


class ShellNode(ABC):
    """Abstract base class for any node in the shell syntax tree."""


class Command(ShellNode):
    """Represents a single command (e.g., 'grep', 'cat')."""

    def __init__(
        self,
        executable: str,
        *args: str | ProcSub,
        fail_ok=False,
        exec_kwargs: Optional[Dict] = None,
    ) -> None:
        self.executable = executable
        self.args = list(args)
        self.fail_ok = fail_ok
        self.proc: Union[asyncio.subprocess.Process, None] = None
        self.exec_kwargs: Dict[str, Any] = exec_kwargs if exec_kwargs else {}

    async def run(
        self,
        context: Context,
        stdin: Union[IO, int, None] = None,
        stdout: Union[IO, int, None] = None,
        stderr: Union[IO, int, None] = None,
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

        # 2. Run the process
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
        return shlex.join([self.executable] + [str(x) for x in self.args])

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({self.executable}, "
            + ", ".join([repr(x) for x in self.args])
            + f", fail_ok={self.fail_ok})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Command):
            return NotImplemented
        return (
            self.executable == other.executable
            and self.args == other.args
            and self.fail_ok == other.fail_ok
            and self.proc == other.proc
        )

    def __hash__(self) -> int:
        return hash(tuple([self.executable, self.fail_ok] + self.args))


class Pipeline(ShellNode):
    """Represents a sequence of commands connected by pipes (|)."""

    def __init__(
        self,
        *nodes: Command,
        skip_pipe: Optional[Iterable[int]] = None,
        file_input: Optional[pathlib.Path] = None,
        file_output: Optional[pathlib.Path] = None,
    ):
        self.nodes = list(nodes)
        assert self.nodes  # Nodes cannot be empty
        self.skip_pipe = set(skip_pipe) if skip_pipe else set()
        self.file_input = file_input
        self.file_output = file_output

    async def run(
        self,
        context: Context,
        stdin: Union[IO, int, None] = None,
        stdout: Union[IO, int, None] = None,
        stderr: Union[IO, int, None] = None,
    ) -> asyncio.subprocess.Process:
        # We cannot have handles from both files and to run()
        if self.file_input and stdin:
            logger.error(
                "Pipeline %s with both file and stdin input. %s",
                self,
                str(stdin),
            )
            sys.exit(1)
        if self.file_output and stdout:
            logger.error(
                "Pipeline %s with both file and stdout output. %s",
                self,
                str(stdout),
            )
            sys.exit(1)

        processes = []
        # We need to close these FDs in the parent after passing them to
        # children
        fds_to_close = []

        # pass the file handle as input to the pipeline
        current_stdin = stdin
        if self.file_input:
            current_stdin = open(self.file_input, "rb")
            context.file_handles.append(current_stdin)

        try:
            # Chain all commands except the last one
            for i in range(len(self.nodes) - 1):
                read_fd = None
                write_fd = None
                if i not in self.skip_pipe:
                    # Create OS pipe
                    read_fd, write_fd = os.pipe()
                    fds_to_close.extend([read_fd, write_fd])

                # Run the node (Command or sub-Pipeline)
                # stdout goes to the write-end of the pipe
                proc = await self.nodes[i].run(
                    context,
                    stdin=current_stdin,
                    stdout=write_fd,
                    stderr=stderr,
                )
                processes.append(proc)

                if write_fd:
                    # Close write_fd in parent (child has it now)
                    os.close(write_fd)
                    fds_to_close.remove(write_fd)

                # Next command reads from the read-end
                current_stdin = read_fd

            # Send the pipeline output to a file
            if self.file_output:
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
        res.append(str(self.nodes[0]))
        for i, node in enumerate(self.nodes[1:]):
            if i in self.skip_pipe:
                res.append("; ")
            else:
                res.append(" | ")
            res.append(str(node))
        if self.file_output:
            res.append(f" >'{self.file_output}'")
        return "".join(res)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            + ", ".join(repr(x) for x in self.nodes)
            + f"skip_pipe={self.skip_pipe},file_input={self.file_input},"
            + f"file_output={self.file_output})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Pipeline):
            return NotImplemented
        return (
            self.nodes == other.nodes
            and self.skip_pipe == other.skip_pipe
            and self.file_input == other.file_input
            and self.file_output == other.file_output
        )

    def __hash__(self) -> int:
        attrs: List[Any] = self.nodes + [self.file_input, self.file_output]
        if self.skip_pipe is None:
            attrs.append(None)
        else:
            attrs.append(tuple(self.skip_pipe))
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
        return self.node == other.node

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.node))


class InputProcSub(ProcSub):
    """
    Represents <(cmd).
    The inner command WRITES to a FIFO.
    The outer command receives the filename.
    """

    async def setup(self, context: Context) -> str:
        self.fifo_path = context.get_new_fifo()

        # We must open the FIFO for writing.
        # WARNING: opening a FIFO blocks until the other end is opened.
        # We run the open() in a thread to prevent blocking the async event
        # loop.
        def open_fifo_write():
            assert self.fifo_path
            return open(self.fifo_path, "w")

        # Start the task that waits for the open, then runs the process
        async def run_inner():
            # This awaits until the outer command opens the file for reading!
            write_handle = await asyncio.to_thread(open_fifo_write)
            try:
                proc = await self.node.run(
                    context,
                    stdout=write_handle,
                )
                await proc.wait()
            finally:
                write_handle.close()

        # Schedule this task in the background
        task = asyncio.create_task(run_inner())
        context.tasks.append(task)

        return self.fifo_path

    def __str__(self) -> str:
        return "<(" + str(self.node) + ")"


class OutputProcSub(ProcSub):
    """
    Represents >(cmd).
    The inner command READS from a FIFO.
    The outer command receives the filename.
    """

    async def setup(self, context: Context) -> str:
        self.fifo_path = context.get_new_fifo()

        def open_fifo_read():
            assert self.fifo_path
            return open(self.fifo_path, "r")

        async def run_inner():
            # This awaits until the outer command opens the file for writing!
            read_handle = await asyncio.to_thread(open_fifo_read)
            try:
                proc = await self.node.run(context, stdin=read_handle)
                await proc.wait()
            finally:
                read_handle.close()

        task = asyncio.create_task(run_inner())
        context.tasks.append(task)

        return self.fifo_path

    def __str__(self) -> str:
        return ">(" + str(self.node) + ")"
