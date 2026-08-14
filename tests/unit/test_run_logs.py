"""
Unit tests for the run log directory and the logging plumbing
"""

import argparse
import logging
import os
import pathlib
import sys
from typing import List

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli import logging as cli_logging  # noqa: E402
from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.exceptions import DagExecutionError  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.pipeline import BasePipeline  # noqa: E402
from sentieon_cli.run_logs import RunLogs  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402
from sentieon_cli.util import __version__  # noqa: E402

PACKAGE_LOGGER = "sentieon_cli"


class _DummyPipeline(BasePipeline):
    """A minimal concrete pipeline for exercising the logging plumbing."""

    def validate(self) -> None:
        pass

    def configure(self) -> None:
        pass

    def build_dag(self) -> DAG:
        return DAG()


class _FailingPipeline(_DummyPipeline):
    """A pipeline whose DAG construction raises."""

    def build_dag(self) -> DAG:
        raise RuntimeError("boom")


class _ErroredExecutor:
    """A stand-in executor that finished with failed jobs."""

    def __init__(self, jobs: List[Job]) -> None:
        self.jobs_with_errors = jobs


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
    package_logger = logging.getLogger(PACKAGE_LOGGER)
    package_logger.addHandler(handler)
    yield handler.messages
    package_logger.removeHandler(handler)


def _args(loglevel: str = "INFO") -> argparse.Namespace:
    return argparse.Namespace(loglevel=loglevel)


def test_module_loggers_delegate_to_the_package_logger():
    module_logger = cli_logging.get_logger("sentieon_cli.example")

    assert module_logger.handlers == []
    assert module_logger.propagate
    assert not logging.getLogger(PACKAGE_LOGGER).propagate


def test_default_log_dir_is_derived_from_the_output_vcf(tmp_path):
    pipeline = _DummyPipeline()
    pipeline.output_vcf = tmp_path / "sample.vcf.gz"
    pipeline.setup_logging(_args())

    assert pipeline.run_logs is not None
    assert pipeline.run_logs.log_dir == tmp_path / "sample_logs"


def test_default_log_dir_only_strips_the_trailing_suffix(tmp_path):
    # `str.replace` would mangle a name with '.vcf.gz' in the middle.
    pipeline = _DummyPipeline()
    pipeline.output_vcf = tmp_path / "a.vcf.gz.rerun.vcf.gz"
    pipeline.setup_logging(_args())

    assert pipeline.run_logs.log_dir == tmp_path / "a.vcf.gz.rerun_logs"


def test_explicit_log_dir_skips_the_output_suffix_check(tmp_path):
    pipeline = _DummyPipeline()
    pipeline.output_vcf = tmp_path / "sample.bcf"
    pipeline.log_dir = tmp_path / "elsewhere"
    pipeline.setup_logging(_args())

    assert pipeline.run_logs.log_dir == tmp_path / "elsewhere"
    assert (tmp_path / "elsewhere" / "run.log").is_file()


def test_invalid_output_suffix_exits_before_creating_a_log_dir(tmp_path):
    pipeline = _DummyPipeline()
    pipeline.output_vcf = tmp_path / "sample.bcf"

    with pytest.raises(SystemExit) as excinfo:
        pipeline.setup_logging(_args())

    assert excinfo.value.code == 2
    assert list(tmp_path.iterdir()) == []


def test_setup_wipes_stale_task_logs_but_keeps_the_log_dir(tmp_path):
    log_dir = tmp_path / "logs"
    stale = log_dir / "task_logs" / "alignment"
    stale.mkdir(parents=True)
    (stale / "bwa-1.0.log").write_text("from the previous run")
    (log_dir / "unrelated.txt").write_text("keep me")

    RunLogs(log_dir).create_dirs()

    assert (log_dir / "unrelated.txt").read_text() == "keep me"
    assert (log_dir / "task_logs").is_dir()
    assert list((log_dir / "task_logs").iterdir()) == []


def test_command_txt_records_the_invocation(tmp_path, monkeypatch):
    monkeypatch.setattr(
        sys, "argv", ["sentieon-cli", "dnascope", "a b.vcf.gz"]
    )
    run_logs = RunLogs(tmp_path / "logs")
    run_logs.create_dirs()
    run_logs.write_command()

    contents = run_logs.command_txt.read_text()
    assert "command: sentieon-cli dnascope 'a b.vcf.gz'" in contents
    assert f"version: {__version__}" in contents
    assert f"directory: {pathlib.Path.cwd()}" in contents
    assert "timestamp: " in contents


def test_run_log_captures_debug_while_the_console_stays_at_info(tmp_path):
    pipeline = _DummyPipeline()
    pipeline.output_vcf = tmp_path / "sample.vcf.gz"
    pipeline.setup_logging(_args("INFO"))
    logging.getLogger("sentieon_cli.example").debug("a debug record")
    run_log = pipeline.run_logs.run_log
    pipeline.run_logs.close()

    contents = run_log.read_text()
    assert "a debug record" in contents
    assert "Starting sentieon-cli version" in contents
    assert f"Writing logs to: {pipeline.run_logs.log_dir}" in contents
    assert cli_logging._console_handler.level == logging.INFO


def test_repeated_setup_does_not_accumulate_handlers(tmp_path):
    package_logger = logging.getLogger(PACKAGE_LOGGER)
    before = len(package_logger.handlers)

    pipeline = None
    for i in range(3):
        pipeline = _DummyPipeline()
        pipeline.output_vcf = tmp_path / f"sample{i}.vcf.gz"
        pipeline.setup_logging(_args())
        assert len(package_logger.handlers) == before + 1

    pipeline.run_logs.close()
    assert len(package_logger.handlers) == before


def test_rerun_truncates_the_previous_run_log(tmp_path):
    log_dir = tmp_path / "logs"
    pipeline = None
    for _ in range(2):
        pipeline = _DummyPipeline()
        pipeline.log_dir = log_dir
        pipeline.setup_logging(_args())
    run_log = pipeline.run_logs.run_log
    pipeline.run_logs.close()

    banner = "Starting sentieon-cli version"
    assert run_log.read_text().count(banner) == 1


def test_dry_run_skips_file_logging(tmp_path):
    pipeline = _DummyPipeline()
    pipeline.dry_run = True
    pipeline.output_vcf = tmp_path / "sample.vcf.gz"
    pipeline.setup_logging(_args())

    assert pipeline.run_logs is None
    assert list(tmp_path.iterdir()) == []


def test_bare_pipeline_setup_logging_creates_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pipeline = _DummyPipeline()
    pipeline.setup_logging(_args())

    assert pipeline.run_logs is None
    assert list(tmp_path.iterdir()) == []


def test_check_execution_names_the_task_log_dir(tmp_path):
    pipeline = _DummyPipeline()
    pipeline.log_dir = tmp_path / "logs"
    pipeline.setup_logging(_args())
    job = Job(Pipeline(Command("false")), "boom", task_name="failing")

    with pytest.raises(DagExecutionError) as excinfo:
        pipeline.check_execution(DAG(), _ErroredExecutor([job]))
    pipeline.run_logs.close()

    message = str(excinfo.value)
    assert "Job(boom-1)" in message
    assert str(tmp_path / "logs" / "task_logs") in message


def test_end_of_run_message_on_success(tmp_path, monkeypatch, messages):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _DummyPipeline()
    pipeline.dry_run = True

    pipeline.main(_args())

    assert any("status: succeeded" in msg for msg in messages)


def test_end_of_run_message_on_failure(tmp_path, monkeypatch, messages):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _FailingPipeline()
    pipeline.log_dir = tmp_path / "logs"

    with pytest.raises(RuntimeError, match="boom"):
        pipeline.main(_args())

    run_log = pipeline.run_logs.run_log
    pipeline.run_logs.close()

    assert any("status: failed" in msg for msg in messages)
    assert any(str(tmp_path / "logs") in msg for msg in messages)
    assert "status: failed" in run_log.read_text()
