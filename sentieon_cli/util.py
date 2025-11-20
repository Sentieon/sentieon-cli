"""
Utility functions
"""

import argparse
import multiprocessing as mp
import os
import pathlib
import re
import shlex
import shutil
import subprocess as sp
import tempfile
from typing import Callable, Dict, List, Optional

import packaging.version

from .logging import get_logger

__version__ = "1.4.0"

logger = get_logger(__name__)

PRELOAD_SEP = r":| "
PRELOAD_SEP_PAT = re.compile(PRELOAD_SEP)

NUMA_NODE_PAT = re.compile(r"^NUMA node. CPU\(s\):\s+(?P<cpus>.*)$")
READ_LENGTH_PAT = re.compile(r"SN\taverage length:\t(?P<length>\d*)$")


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
    elif cmd_list[0] == "pbsv":
        cmd_version_str = cmd_version_str.split(" ")[1]
    elif cmd_list[0] == "hificnv":
        cmd_version_str = cmd_version_str.split(" ")[1].split("-")[0]
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
    numa_nodes = []
    try:
        res = sp.run("lscpu", capture_output=True, text=True)
        for line in res.stdout.split("\n"):
            m = NUMA_NODE_PAT.match(line)
            if m:
                cpus = m.groupdict()["cpus"]
                numa_nodes.append(cpus)
        logger.debug("Identified NUMA nodes: %s", numa_nodes)
    except FileNotFoundError:
        numa_nodes = ["0-" + str(mp.cpu_count() - 1)]
    return numa_nodes


def split_numa_nodes(numa_nodes: List[str]) -> List[str]:
    """Split numa nodes in half"""
    new_numa_nodes = []
    for numa_node in numa_nodes:
        if "," in numa_node:
            ranges = numa_node.split(",")
            mid = len(ranges) // 2
            new_numa_nodes.append(",".join(ranges[:mid]))
            new_numa_nodes.append(",".join(ranges[mid:]))
        else:
            start, end = map(int, numa_node.split("-"))
            mid = (start + end) // 2
            new_numa_nodes.append(f"{start}-{mid}")
            new_numa_nodes.append(f"{mid + 1}-{end}")
    return new_numa_nodes


def split_alignment(cores: int) -> List[str]:
    """split large alignment tasks into smaller parts on large machines"""
    numa_nodes = find_numa_nodes()
    while len(numa_nodes) > 0 and cores / len(numa_nodes) > 64:
        numa_nodes = split_numa_nodes(numa_nodes)
    if cores > 32 and numa_nodes:
        return numa_nodes
    else:
        return []


def parse_rg_line(rg_line: str) -> Dict[str, str]:
    """Parse an @RG line"""
    tags = rg_line.split("\t")[1:]
    return {tag.split(":", 1)[0]: tag.split(":", 1)[1] for tag in tags}


def get_read_length_aln(
    aln: pathlib.Path,
    reference: pathlib.Path,
    n_reads: int = 100000,
) -> int:
    """Get the average read length for an alignment file"""
    cmds = []
    cmds.append(
        [
            "samtools",
            "view",
            "-h",
            "--reference",
            shlex.quote(str(reference)),
            shlex.quote(str(aln)),
        ]
    )
    cmds.append(
        [
            "head",
            "-n",
            str(n_reads),
        ]
    )
    cmds.append(["samtools", "stats", "-"])
    all_cmds = [shlex.join(x) for x in cmds]
    cmd = " | ".join(all_cmds)
    res = sp.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        executable="/bin/bash",
    )

    for line in res.stdout.split("\n"):
        m = READ_LENGTH_PAT.match(line)
        if m:
            return int(m.groupdict()["length"])
    return 151
