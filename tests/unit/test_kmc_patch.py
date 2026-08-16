"""
Tests for KMC patch check
"""

import subprocess as sp
from unittest.mock import patch

from sentieon_cli.util import check_kmc_patch


def kmc_stdout(n_kmers: int) -> str:
    """The stats block printed by KMC after a successful run"""
    return (
        "1st stage: 0.068111s\n"
        "2nd stage: 0.01257s\n"
        "Total    : 0.080681s\n"
        "Tmp size : 0MB\n"
        "\n"
        "Stats:\n"
        f"   No. of k-mers below min. threshold : {n_kmers:>12}\n"
        "   No. of k-mers above max. threshold :            0\n"
        f"   No. of unique k-mers               : {n_kmers:>12}\n"
        "   No. of unique counted k-mers       :            0\n"
        f"   Total no. of k-mers                : {n_kmers:>12}\n"
        "   Total no. of reads                 :            2\n"
        "   Total no. of super-k-mers          :            3\n"
    )


def test_check_kmc_patch_success():
    """Test check_kmc_patch returns True when k-mers are counted"""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = kmc_stdout(13)
        assert check_kmc_patch() is True
        mock_run.assert_called_once()
        args, kwargs = mock_run.call_args
        assert "kmc" == args[0][0]
        assert kwargs["input"] is not None
        # The k-mer count is parsed from the captured stdout
        assert kwargs["stdout"] == sp.PIPE


def test_check_kmc_patch_zero_kmers():
    """An unpatched KMC reads nothing from stdin but still exits 0"""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = kmc_stdout(0)
        assert check_kmc_patch() is False


def test_check_kmc_patch_missing_stats():
    """Test check_kmc_patch returns False without the k-mer count"""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = "Stage 1: 100%\nStage 2: 100%\n"
        assert check_kmc_patch() is False


def test_check_kmc_patch_no_stdout():
    """Test check_kmc_patch returns False without any output"""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = None
        assert check_kmc_patch() is False


def test_check_kmc_patch_custom_command():
    """Test check_kmc_patch runs the supplied kmc executable"""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = kmc_stdout(13)
        assert check_kmc_patch("/opt/kmc/kmc") is True
        args, _kwargs = mock_run.call_args
        assert args[0][0] == "/opt/kmc/kmc"
        # KMC reads the test input from stdin
        assert "/dev/stdin" in args[0]


def test_check_kmc_patch_failure():
    """Test check_kmc_patch returns False on failure"""
    with patch("subprocess.run") as mock_run:
        mock_run.side_effect = sp.CalledProcessError(1, "kmc")
        assert check_kmc_patch() is False


def test_check_kmc_patch_not_found():
    """Test check_kmc_patch returns False if kmc is missing"""
    with patch("subprocess.run") as mock_run:
        mock_run.side_effect = FileNotFoundError
        assert check_kmc_patch() is False
