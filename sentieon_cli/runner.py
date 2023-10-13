import subprocess as sp
import sys
import time

from .logging import get_logger

logger = get_logger(__name__)


def run(cmd: str):
    """Run a command."""
    logger.debug("running: %s", cmd)
    t0 = time.time()
    sp.run(cmd, shell=True, check=True, stdout=sys.stdout, stderr=sys.stderr)
    logger.debug("finished in: %s", tdiff = f"{time.time() - t0: .1f}")
