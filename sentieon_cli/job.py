"""
Job objects
"""

import subprocess as sp
import sys
import time

from .logging import get_logger

logger = get_logger(__name__)


class Job:
    """A job for execution"""

    def __init__(
        self, shell: str, name: str, threads: int = 1, fail_ok: bool = False
    ):
        self.shell = shell
        self.name = name
        self.threads = threads
        self.fail_ok = fail_ok

    def __hash__(self):
        return hash(self.shell)

    def __eq__(self, other: object):
        if isinstance(other, Job):
            return self.shell == other.shell
        return False

    def __ne__(self, other: object):
        return not self == other

    def __repr__(self):
        return f"Job({self.name})"

    def __str__(self):
        return f"Job({self.name})"

    def run(self, dry_run: bool = False):
        """Run a command"""
        if dry_run:
            print(self.shell)
            return

        logger.info("running command: %s", self.shell)
        t0 = time.time()
        sp.run(
            self.shell,
            shell=True,
            check=False,
            stdout=sys.stdout,
            stderr=sys.stderr,
            executable="/bin/bash",
        )
        logger.info("finished in: %s seconds", f"{time.time() - t0:.1f}")
