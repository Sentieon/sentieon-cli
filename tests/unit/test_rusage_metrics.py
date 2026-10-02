"""
Unit tests for the per-process resource-usage metrics
"""

import logging
import os
import pathlib
import resource
import sys
from typing import List

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.executor import (  # noqa: E402
    LocalExecutor,
    _maxrss_bytes,
)
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.run_logs import RunLogs  # noqa: E402
from sentieon_cli.scheduler import ThreadScheduler  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402


class _RecordingHandler(logging.Handler):
    """Collects formatted messages from the package logger."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


@pytest.fixture
def messages():
    """Capture the package logger, which does not propagate to the root."""
    handler = _RecordingHandler()
    package_logger = logging.getLogger("sentieon_cli")
    package_logger.addHandler(handler)
    yield handler.messages
    package_logger.removeHandler(handler)


def _run_one_job(log_dir: pathlib.Path) -> None:
    """Run a single trivial job with its output captured under ``log_dir``"""
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("true")), "metrics", task_name="metrics"))
    run_logs = RunLogs(log_dir)
    run_logs.create_dirs()
    LocalExecutor(ThreadScheduler(dag, 1), run_logs=run_logs).execute()


def test_each_process_of_a_job_reports_its_rusage(tmp_path, messages):
    _run_one_job(tmp_path / "logs")

    rusage_lines = [m for m in messages if m.startswith("rusage for ")]
    assert len(rusage_lines) == 1
    line = rusage_lines[0]
    assert "Job(metrics-1)" in line
    for field in ("utime=", "stime=", "maxrss=", "minflt=", "nivcsw="):
        assert field in line


def test_a_finished_job_reports_its_aggregated_metrics(tmp_path, messages):
    _run_one_job(tmp_path / "logs")

    finished = [m for m in messages if m.startswith("Finished command in")]
    assert len(finished) == 1
    assert "user " in finished[0]
    assert "sys " in finished[0]
    assert "max proc RSS" in finished[0]


def test_maxrss_is_reported_in_bytes():
    ru = resource.getrusage(resource.RUSAGE_SELF)

    maxrss = _maxrss_bytes(ru)

    assert isinstance(maxrss, int)
    assert maxrss > 0


def test_maxrss_units_follow_the_platform(monkeypatch):
    """macOS reports bytes; everywhere else reports kibibytes."""
    ru = resource.getrusage(resource.RUSAGE_SELF)

    monkeypatch.setattr(sys, "platform", "darwin")
    assert _maxrss_bytes(ru) == ru.ru_maxrss

    monkeypatch.setattr(sys, "platform", "linux")
    assert _maxrss_bytes(ru) == ru.ru_maxrss * 1024
