"""
Shared integration-test fixtures
"""

import pytest

from sentieon_cli.job import Job

# `pytest --doctest-modules` collects .py files as modules, and two
# conftest.py files outside a package collide on the module name. This one
# holds no doctests, so skip collecting it.
collect_ignore = ["conftest.py"]


@pytest.fixture(autouse=True)
def reset_job_ids():
    """Restart job id numbering for every test.

    ``job_id`` comes from a class-level per-name counter that the CLI resets
    once per invocation, so without this a test's ids would depend on which
    tests ran before it.
    """
    Job.reset_ids()
