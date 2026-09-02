"""
Utility functions
"""

import argparse
from enum import Enum
from importlib.metadata import PackageNotFoundError, version
import multiprocessing as mp
import os
import pathlib
import re
import shutil
import subprocess as sp
import sys
import tempfile
from typing import Callable, Dict, List, Mapping, Optional, Tuple

import packaging.version

from .logging import get_logger

try:
    __version__ = version("sentieon_cli")
except PackageNotFoundError:
    __version__ = "0.0.0"

logger = get_logger(__name__)

PRELOAD_SEP = r":| "
PRELOAD_SEP_PAT = re.compile(PRELOAD_SEP)

NUMA_NODE_PAT = re.compile(r"^NUMA node. CPU\(s\):\s+(?P<cpus>.*)$")

_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")


class SampleSex(Enum):
    FEMALE = 1
    MALE = 2
    UNKNOWN = 3


def sample_sex_arg(value: str) -> SampleSex:
    """Parse the `--sample_sex` argument"""
    sex = value.strip().lower()
    if sex == "male":
        return SampleSex.MALE
    if sex == "female":
        return SampleSex.FEMALE
    raise argparse.ArgumentTypeError(
        f"invalid sample sex '{value}'. Please supply 'male' or 'female'"
    )


def cnvscope_sex_args(
    sample_sex: Optional[SampleSex],
    par_bed: Optional[pathlib.Path],
) -> Tuple[Optional[str], Optional[pathlib.Path]]:
    """The CNVscope `--sex` and `--par` arguments for a sample.

    Male samples are called with the pseudo-autosomal regions (PAR) BED
    file. Female samples do not need a PAR BED file. When the sample sex
    is not known, both arguments are omitted and CNVscope assumes a
    diploid genome, matching the behavior of previous releases.
    """
    if sample_sex is SampleSex.MALE:
        return ("M", par_bed)
    if sample_sex is SampleSex.FEMALE:
        return ("F", None)
    logger.warning(
        "The sample sex is not known. CNVscope will assume a diploid "
        "genome, matching the behavior of previous releases. Supply "
        "`--sample_sex` for sex-aware CNV calling."
    )
    return (None, None)


def sanitize(component: str) -> str:
    """Restrict a path component to filesystem-safe characters"""
    return _UNSAFE.sub("-", component)


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
    try:
        cmd_version_str = (
            sp.check_output(cmd_list).decode("utf-8", "ignore").strip()
        )
    except (sp.CalledProcessError, OSError) as e:
        logger.error(
            "Error: could not determine the version of '%s': %s", cmd, e
        )
        return False
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


def require_versions(
    min_versions: Mapping[str, Optional[packaging.version.Version]],
    *,
    skip: bool = False,
) -> None:
    """Exit unless every executable meets its minimum version.

    `check_version` has already logged the reason, so the exit is silent.
    Pass `skip=True` (the pipelines' `--skip_version_check`) to do nothing.

    A `Mapping` rather than a `Dict`, so the pipelines' module-level
    `*_MIN_VERSIONS` constants -- some of which mypy infers as
    `Dict[str, Version]`, with no `None` entry -- are accepted.
    """
    if skip:
        return
    for cmd, min_version in min_versions.items():
        if not check_version(cmd, min_version):
            sys.exit(2)


def versions_available(
    min_versions: Mapping[str, Optional[packaging.version.Version]],
    *,
    skip: bool = False,
) -> bool:
    """Whether every executable meets its minimum version.

    For optional tools, where a missing or outdated executable skips a step
    rather than ending the run. `skip=True` reports them as available.
    """
    if skip:
        return True
    return all(
        check_version(cmd, min_version)
        for cmd, min_version in min_versions.items()
    )


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
    """Parse an @RG line.

    Each tab-separated field after ``@RG`` must be of the form
    ``KEY:VALUE``. Raises ``ValueError`` with the offending token if a
    field lacks the ``KEY:`` prefix (e.g. a bare ``HG002`` instead of
    ``SM:HG002``).
    """
    tags = rg_line.split("\t")[1:]
    parsed: Dict[str, str] = {}
    for tag in tags:
        key, sep, value = tag.partition(":")
        if not sep or not key:
            raise ValueError(
                f"malformed @RG field '{tag}' — expected 'KEY:VALUE' "
                f"(e.g. 'SM:sample')"
            )
        parsed[key] = value
    return parsed


def vcf_id(in_vcf: pathlib.Path) -> Optional[str]:
    """Collect the SentieonVcfID header"""
    cmd = ["bcftools", "view", "-h", str(in_vcf)]
    p = sp.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        logger.error(
            "`%s` failed with return code %d: %s",
            " ".join(cmd),
            p.returncode,
            p.stderr.strip(),
        )
        return None
    for line in p.stdout.split("\n"):
        if line.startswith("##SentieonVcfID="):
            i = line.index("=")
            return line[i + 1 :]  # noqa: E203
    return None


def check_kmc_patch(kmc_cmd: str = "kmc") -> bool:
    """Check if the KMC version supports the required patch"""
    # Create a temporary directory for the test
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = pathlib.Path(temp_dir)
        output_prefix = temp_path / "kmc_test"

        # Test input sequence
        test_input = (
            ">206B4ABXX100825:7:1:1360:6029/1\n"
            "TGATTTTNNNNNNNNNNNTGAAGAACGCACCCATGTTAAAGAGCATGACAAANNNANNACAAGGCTAAGNGGCGNG\n"  # noqa: E501
            ">206B4ABXX100825:7:1:1362:4449/1\n"
            "ATTCCCCNNNNNNNNNNNCCACAGCCGGAGGAGCTGACCAACATCCTGGAGATNTGNAATGTGGTCTTANCCAGNA"  # noqa: E501
        )

        cmd = [
            kmc_cmd,
            "-k29",
            "-m4",
            "-okff",
            "-t1",
            "-fa",
            "/dev/stdin",
            str(output_prefix),
            str(temp_path),
        ]

        try:
            sp.run(
                cmd,
                input=test_input,
                text=True,
                check=True,
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL,
            )
            return True
        except (sp.CalledProcessError, FileNotFoundError):
            return False
