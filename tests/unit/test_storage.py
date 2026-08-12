"""
Unit tests for the storage providers.
"""

import os
import sys

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.storage import LocalStorageProvider  # noqa: E402


@pytest.mark.asyncio
async def test_stage_in_present_file_ok(tmp_path):
    f = tmp_path / "in.txt"
    f.write_text("data")
    await LocalStorageProvider().stage_in(f)  # no raise


@pytest.mark.asyncio
async def test_stage_in_empty_file_ok(tmp_path):
    f = tmp_path / "in.txt"
    f.touch()  # a zero-byte file still counts as present
    await LocalStorageProvider().stage_in(f)


@pytest.mark.asyncio
async def test_stage_in_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        await LocalStorageProvider().stage_in(tmp_path / "nope.txt")


@pytest.mark.asyncio
async def test_stage_out_present_file_ok(tmp_path):
    f = tmp_path / "out.txt"
    f.touch()
    await LocalStorageProvider().stage_out(f)


@pytest.mark.asyncio
async def test_stage_out_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        await LocalStorageProvider().stage_out(tmp_path / "nope.txt")
