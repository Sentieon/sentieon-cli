"""
Unit tests for the executable version helpers in `sentieon_cli.util`
"""

import subprocess as sp
from unittest.mock import MagicMock

import packaging.version

from sentieon_cli import util


class TestExecutableVersion:
    """Test `util.executable_version`"""

    def test_sentieon_driver_version(self, monkeypatch):
        """The sentieon release is taken from the package name"""
        check_output = MagicMock(return_value=b"sentieon-genomics-202503.04\n")
        monkeypatch.setattr(util.sp, "check_output", check_output)

        assert util.executable_version(
            "sentieon driver"
        ) == packaging.version.Version("202503.04")
        assert check_output.call_args[0][0] == [
            "sentieon",
            "driver",
            "--version",
        ]

    def test_multiline_version(self, monkeypatch):
        """Only the first line of a multi-line `--version` is parsed"""
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(
                return_value=(
                    b"bcftools 1.22\n"
                    b"Using htslib 1.22\n"
                    b"Copyright (C) 2024 Genome Research Ltd.\n"
                )
            ),
        )

        assert util.executable_version(
            "bcftools"
        ) == packaging.version.Version("1.22")

    def test_called_process_error(self, monkeypatch):
        """A failing `--version` call reports no version"""
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(side_effect=sp.CalledProcessError(1, "bcftools")),
        )

        assert util.executable_version("bcftools") is None

    def test_os_error(self, monkeypatch):
        """An un-runnable command reports no version"""
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(side_effect=OSError("no such file")),
        )

        assert util.executable_version("bcftools") is None

    def test_unparseable_version(self, monkeypatch):
        """Output that is not a version reports no version"""
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(return_value=b"bcftools unknown-build\n"),
        )

        assert util.executable_version("bcftools") is None


class TestCheckVersion:
    """Test `util.check_version`"""

    def test_missing_executable(self, monkeypatch):
        """A command that is not in the PATH fails the check"""
        monkeypatch.setattr(util.shutil, "which", MagicMock(return_value=None))
        check_output = MagicMock()
        monkeypatch.setattr(util.sp, "check_output", check_output)

        assert (
            util.check_version("bcftools", packaging.version.Version("1.22"))
            is False
        )
        check_output.assert_not_called()

    def test_no_minimum_version_skips_probe(self, monkeypatch):
        """Tools listed without a minimum version are not probed"""
        monkeypatch.setattr(
            util.shutil, "which", MagicMock(return_value="/usr/bin/vg")
        )
        check_output = MagicMock()
        monkeypatch.setattr(util.sp, "check_output", check_output)

        assert util.check_version("vg", None) is True
        check_output.assert_not_called()

    def test_old_version(self, monkeypatch):
        """An outdated executable fails the check"""
        monkeypatch.setattr(
            util.shutil, "which", MagicMock(return_value="/usr/bin/bcftools")
        )
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(return_value=b"bcftools 1.16\n"),
        )

        assert (
            util.check_version("bcftools", packaging.version.Version("1.22"))
            is False
        )

    def test_current_version(self, monkeypatch):
        """An up-to-date executable passes the check"""
        monkeypatch.setattr(
            util.shutil, "which", MagicMock(return_value="/usr/bin/bcftools")
        )
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(return_value=b"bcftools 1.22\n"),
        )

        assert (
            util.check_version("bcftools", packaging.version.Version("1.22"))
            is True
        )

    def test_unparseable_version_fails(self, monkeypatch):
        """Unparseable `--version` output fails the check"""
        monkeypatch.setattr(
            util.shutil, "which", MagicMock(return_value="/usr/bin/bcftools")
        )
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(return_value=b"bcftools unknown-build\n"),
        )

        assert (
            util.check_version("bcftools", packaging.version.Version("1.22"))
            is False
        )
