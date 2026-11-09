"""
Unit tests for the check for Sentieon's ExpansionHunter
"""

from sentieon_cli.util import check_expansion_hunter_read_mode


def stub_help(tmp_path, text, exit_code=0):
    """A stub executable printing `text` to stderr, as ExpansionHunter does"""
    path = tmp_path / "ExpansionHunter"
    path.write_text(f"#!/bin/sh\necho '{text}' >&2\nexit {exit_code}\n")
    path.chmod(0o755)
    return str(path)


def test_read_mode_found_in_stderr_help(tmp_path):
    cmd = stub_help(tmp_path, "  --read-mode arg (=auto)  Read model")
    assert check_expansion_hunter_read_mode(cmd)


def test_read_mode_missing_from_help(tmp_path):
    cmd = stub_help(tmp_path, "  --variant-catalog arg  JSON file")
    assert not check_expansion_hunter_read_mode(cmd)


def test_read_mode_ignores_the_exit_code(tmp_path):
    cmd = stub_help(tmp_path, "  --read-mode arg", exit_code=1)
    assert check_expansion_hunter_read_mode(cmd)


def test_read_mode_missing_executable(tmp_path):
    assert not check_expansion_hunter_read_mode(str(tmp_path / "missing"))
