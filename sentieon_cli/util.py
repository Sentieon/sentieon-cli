"""
Utility functions
"""

import argparse
import os
import pathlib
import re
import shutil
import subprocess as sp
import tempfile
from typing import Callable, List, Optional

import packaging.version

from .logging import get_logger

__version__ = "1.1.0"

logger = get_logger(__name__)

PRELOAD_SEP = r":| "
PRELOAD_SEP_PAT = re.compile(PRELOAD_SEP)

NUMA_NODE_PAT = re.compile(r"^NUMA node. CPU\(s\):   (?P<cpus>.*)$")


def tmp():
    """Create a temporary directory for the current process."""
    tmp_base = os.getenv("SENTIEON_TMPDIR")
    tmp_dir = tempfile.mkdtemp(dir=tmp_base)
    return tmp_dir


def check_version(
    cmd: str,
    version: Optional[packaging.version.Version],
) -> bool:
    """Check the version of an executable"""
    cmd_list: List[str] = cmd.split()
    exec_file = shutil.which(cmd_list[0])
    if not exec_file:
        logger.error("Error: no '%s' found in the PATH", cmd)
        return False

    if version is None:
        return True

    cmd_list.append("--version")
    cmd_version_str = (
        sp.check_output(cmd_list).decode("utf-8", "ignore").strip()
    )
    if cmd_list[0] == "sentieon":
        cmd_version_str = cmd_version_str.split("-")[-1]
    else:
        # handle, e.g. bcftools which outputs multiple lines.
        cmd_version_str = (
            cmd_version_str.split("\n")[0].split()[-1].split("-")[0]
        )
    cmd_version = packaging.version.Version(cmd_version_str)
    if cmd_version < version:
        logger.error(
            "Error: the pipeline requires %s version '%s' or later "
            "but %s '%s' was found in the PATH",
            cmd,
            version,
            cmd,
            cmd_version,
        )
        return False
    return True


def path_arg(
    exists: Optional[bool] = None,
    is_dir: Optional[bool] = None,
    is_file: Optional[bool] = None,
    is_fifo: Optional[bool] = None,
) -> Callable[[str], pathlib.Path]:
    """pathlib checked types for argparse"""

    def _path_arg(arg: str) -> pathlib.Path:
        p = pathlib.Path(arg)

        attrs = [exists, is_dir, is_file, is_fifo]
        attr_names = ["exists", "is_dir", "is_file", "is_fifo"]

        for attr_val, attr_name in zip(attrs, attr_names):
            if attr_val is None:  # Skip attributes that are not defined
                continue

            m = getattr(p, attr_name)
            if m() != attr_val:
                raise argparse.ArgumentTypeError(
                    "The supplied path argument needs the attribute"
                    f" {attr_name}={attr_val}, but {attr_name}={m()}"
                )
        return p

    return _path_arg


def library_preloaded(library_name: str) -> bool:
    """Check if a shared library is preloaded through LD_PRELOAD"""
    ld_preload = os.getenv("LD_PRELOAD", "")
    for lib in PRELOAD_SEP_PAT.split(ld_preload):
        lib_base = os.path.basename(lib)
        if library_name in lib_base:
            return True
    return False


def total_memory() -> int:
    """The total memory accounting for cgroup limits"""
    total_mem = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")

    # Attempt to find cgroup limits
    cgroup_mem_limit = 10000 * 1024**3
    cgroup_files = [
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    ]
    for cgroup_file_s in cgroup_files:
        cgroup_file = pathlib.Path(cgroup_file_s)
        if cgroup_file.is_file():
            try:
                cgroup_mem_limit = int(open(cgroup_file).read().rstrip())
            except Exception:
                pass
    total_mem = min(total_mem, cgroup_mem_limit)
    return total_mem


def find_numa_nodes() -> List[str]:
    """Find NUMA nodes on the system, if available"""
    res = sp.run("lscpu", capture_output=True, text=True)
    numa_nodes = []
    for line in res.stdout.split("\n"):
        m = NUMA_NODE_PAT.match(line)
        if m:
            cpus = m.groupdict()["cpus"]
            numa_nodes.append(cpus)
    logger.debug("Identified NUMA nods: %s", numa_nodes)
    return numa_nodes
