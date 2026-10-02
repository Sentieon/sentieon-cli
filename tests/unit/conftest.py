"""
Shared unit-test fixtures
"""

import logging

import pytest

from sentieon_cli import logging as cli_logging
from sentieon_cli.job import Job


@pytest.fixture(autouse=True)
def reset_job_ids():
    """Restart job id numbering for every test.

    ``job_id`` comes from a class-level per-name counter that the CLI resets
    once per invocation, so without this a test's ids would depend on which
    tests ran before it.
    """
    Job.reset_ids()


@pytest.fixture(autouse=True)
def no_handler_leaks():
    """Drop any handler a test leaves on the package logger.

    `setup_logging` attaches a `run.log` FileHandler to the package logger,
    so any test that reaches it with an output VCF (or log_dir) set would
    otherwise leak that handler into later tests.
    """
    package_logger = logging.getLogger("sentieon_cli")
    before = list(package_logger.handlers)
    console_level = cli_logging._console_handler.level
    yield
    for handler in list(package_logger.handlers):
        if handler in before:
            continue
        if isinstance(handler, logging.FileHandler):
            cli_logging.remove_file_handler(handler)
        else:
            package_logger.removeHandler(handler)
    cli_logging.set_console_level(console_level)
