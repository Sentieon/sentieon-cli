"""
Tests for KMC patch check
"""

import subprocess as sp
from unittest.mock import MagicMock, patch

from sentieon_cli.util import check_kmc_patch


def test_check_kmc_patch_success():
    """Test check_kmc_patch returns True on success"""
    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        assert check_kmc_patch() is True
        mock_run.assert_called_once()
        args, kwargs = mock_run.call_args
        assert "kmc" == args[0][0]
        assert kwargs["input"] is not None


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
