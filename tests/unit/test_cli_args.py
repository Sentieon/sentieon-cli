"""
Unit tests for the top-level CLI arguments
"""

import argparse
import os
import sys

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli import add_logging_args, resolve_loglevel  # noqa: E402


def _loglevel(argv):
    parser = argparse.ArgumentParser()
    add_logging_args(parser)
    return resolve_loglevel(parser.parse_args(argv))


@pytest.mark.parametrize(
    "argv,expected",
    [
        ([], "INFO"),
        (["-v"], "INFO"),
        (["-q"], "WARNING"),
        (["-d"], "DEBUG"),
        # `-d` wins wherever it appears; it used to be undone by a later
        # `-v` because both flags wrote to one destination.
        (["-d", "-v"], "DEBUG"),
        (["-q", "-v"], "INFO"),
    ],
)
def test_the_verbosity_flags_resolve_to_a_console_level(argv, expected):
    assert _loglevel(argv) == expected
