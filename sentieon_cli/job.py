"""
Job objects
"""

import asyncio
import sys
import time
from typing import Dict, Optional

from .logging import get_logger
from .shell_pipeline import Context, Pipeline

logger = get_logger(__name__)


class Job:
    """A job for execution"""

    def __init__(
        self,
        pipeline: Pipeline,
        name: str,
        threads: int = 1,
        resources: Optional[Dict[str, int]] = None,
    ):
        self.shell = pipeline
        self.name = name
        self.threads = threads
        self.resources = {} if resources is None else resources

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
        context = Context()
        try:
            asyncio.run(
                self.shell.run(context, stdout=sys.stdout, stderr=sys.stderr)
            )
            for subcommand in context.commands:
                assert subcommand.proc
                ret = asyncio.run(subcommand.proc.wait())
                if ret != 0 and not subcommand.fail_ok:
                    logger.error(
                        f"subcommand failed with code {ret}: {subcommand}"
                    )
        except Exception as e:
            logger.error(f"Failed to run command: {e}")
        finally:
            asyncio.run(context.cleanup())
        logger.info("finished in: %s seconds", f"{time.time() - t0:.1f}")
