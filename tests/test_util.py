import os
import subprocess as sp
import sys
from unittest.mock import patch

import packaging.version

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

import sentieon_cli.util  # NOQA


def test_split_numa_nodes():
    """Test the NUMA split func"""
    res = sentieon_cli.util.split_numa_nodes(["0-97"])
    assert res == ["0-48", "49-97"]
    res = sentieon_cli.util.split_numa_nodes(res)
    assert res == ["0-24", "25-48", "49-73", "74-97"]

    res = sentieon_cli.util.split_numa_nodes(["0-15,32-47", "16-31,48-63"])
    assert res == ["0-15", "32-47", "16-31", "48-63"]
    res = sentieon_cli.util.split_numa_nodes(res)
    assert res == [
        "0-7",
        "8-15",
        "32-39",
        "40-47",
        "16-23",
        "24-31",
        "48-55",
        "56-63",
    ]


def test_check_version_parses_version():
    """A well-behaved `--version` is compared against the minimum"""
    min_version = packaging.version.Version("0.9.0")
    with (
        patch.object(
            sentieon_cli.util.shutil,
            "which",
            return_value="/bin/segdup-caller",
        ),
        patch.object(
            sentieon_cli.util.sp,
            "check_output",
            return_value=b"segdup-caller 0.9.0\n",
        ),
    ):
        assert sentieon_cli.util.check_version("segdup-caller", min_version)

    with (
        patch.object(
            sentieon_cli.util.shutil,
            "which",
            return_value="/bin/segdup-caller",
        ),
        patch.object(
            sentieon_cli.util.sp,
            "check_output",
            return_value=b"segdup-caller 0.8.0\n",
        ),
    ):
        assert not sentieon_cli.util.check_version(
            "segdup-caller", min_version
        )


def test_check_version_handles_failed_command():
    """A non-zero `--version` reports False rather than raising.

    segdup-caller >=0.6.0 exits non-zero on `--version` when
    SENTIEON_LICENSE is unset.
    """
    min_version = packaging.version.Version("0.9.0")
    for err in (
        sp.CalledProcessError(1, ["segdup-caller", "--version"]),
        OSError("Exec format error"),
    ):
        with (
            patch.object(
                sentieon_cli.util.shutil,
                "which",
                return_value="/bin/segdup-caller",
            ),
            patch.object(
                sentieon_cli.util.sp, "check_output", side_effect=err
            ),
        ):
            assert not sentieon_cli.util.check_version(
                "segdup-caller", min_version
            )
