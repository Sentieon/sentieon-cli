"""
The log directory for a single run
"""

import datetime
import logging
import pathlib
import shlex
import shutil
import sys
from typing import Optional

from .logging import add_file_handler, remove_file_handler
from .util import __version__


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
