"""
Unit tests for the KMC thread cap
"""

import pathlib
import re

from sentieon_cli.command_strings import (
    KMC_MAX_THREADS,
    cmd_extract_kmc,
    cmd_hybrid_kmc,
    cmd_kmc,
)

TMP = pathlib.Path("/tmp/kmc")
PREFIX = pathlib.Path("/tmp/sample")


def kmc_t(cmd_str: str) -> str:
    """The -t value of the kmc command in a pipeline"""
    match = re.search(r"\bkmc -k\S+ -m\S+ -okff -t(\d+)", cmd_str)
    assert match, cmd_str
    return match.group(1)


def test_kmc_threads_capped():
    """More cores than KMC accepts are capped at its limit"""
    cmd = cmd_kmc(PREFIX, pathlib.Path("/tmp/fq.list"), TMP, threads=192)
    assert kmc_t(str(cmd)) == str(KMC_MAX_THREADS)


def test_kmc_threads_below_cap_unchanged():
    """A thread count KMC accepts is passed through"""
    cmd = cmd_kmc(PREFIX, pathlib.Path("/tmp/fq.list"), TMP, threads=64)
    assert kmc_t(str(cmd)) == "64"


def test_extract_kmc_caps_only_kmc():
    """The driver and pgutil keep every core; only KMC is capped"""
    cmd_str = str(
        cmd_extract_kmc(
            PREFIX,
            pathlib.Path("/tmp/ext.fq.gz"),
            [pathlib.Path("/data/sample.bam")],
            pathlib.Path("/ref/hg38.fa"),
            pathlib.Path("/bundle/extract.model"),
            TMP,
            pathlib.Path("/tmp/rw.bam"),
            threads=192,
        )
    )
    assert kmc_t(cmd_str) == str(KMC_MAX_THREADS)
    assert "-t 192" in cmd_str


def test_hybrid_kmc_threads_capped():
    """The stdin KMC used by the hybrid pipelines is capped too"""
    cmd = cmd_hybrid_kmc(
        PREFIX, [pathlib.Path("/data/r1.fq.gz")], [], TMP, threads=192
    )
    assert kmc_t(str(cmd)) == str(KMC_MAX_THREADS)
