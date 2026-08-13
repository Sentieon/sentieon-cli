"""
Unit tests for the BasePipeline lifecycle (argument registration, logging
scope, temp-dir cleanup, and post-run checks).
"""

import argparse
import logging
import os
import sys

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.exceptions import DagExecutionError  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.pipeline import BasePipeline  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402


class _DummyPipeline(BasePipeline):
    """A minimal concrete pipeline for exercising BasePipeline plumbing."""

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


def _option_strings(parser):
    return [opt for action in parser._actions for opt in action.option_strings]


def test_add_arguments_does_not_mutate_the_class_specs():
    before_reference = dict(BasePipeline.params["reference"])
    before_cores = dict(BasePipeline.params["cores"])

    parser = argparse.ArgumentParser()
    BasePipeline.add_arguments(parser)

    assert BasePipeline.params["reference"] == before_reference
    assert BasePipeline.params["cores"] == before_cores
    # In particular the specs keep their flags and gain no inferred type.
    assert "flags" in BasePipeline.params["reference"]
    assert "type" not in BasePipeline.params["cores"]


def test_add_arguments_keeps_short_flags_across_parsers():
    # Registering the same pipeline class on a second parser must produce
    # the same options; the short flags used to be lost because the first
    # call popped them from the shared class-level spec.
    p1 = argparse.ArgumentParser()
    p2 = argparse.ArgumentParser()
    BasePipeline.add_arguments(p1)
    BasePipeline.add_arguments(p2)

    for parser in (p1, p2):
        opts = _option_strings(parser)
        assert "-r" in opts
        assert "--reference" in opts
        assert "-t" in opts
        assert "--cores" in opts


def test_setup_logging_does_not_touch_the_root_logger():
    root_level_before = logging.getLogger().level
    package_logger = logging.getLogger("sentieon_cli")
    package_level_before = package_logger.level
    try:
        pipeline = _DummyPipeline()
        pipeline.setup_logging(argparse.Namespace(loglevel=logging.DEBUG))

        assert logging.getLogger().level == root_level_before
        assert package_logger.level == logging.DEBUG
    finally:
        package_logger.setLevel(package_level_before)


def test_main_cleans_up_tmpdir_when_build_dag_raises(tmp_path, monkeypatch):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _FailingPipeline()

    with pytest.raises(RuntimeError, match="boom"):
        pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert list(tmp_path.iterdir()) == []  # no leaked temp dir


def test_main_retains_tmpdir_on_failure_when_requested(tmp_path, monkeypatch):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _FailingPipeline()
    pipeline.retain_tmpdir = True

    with pytest.raises(RuntimeError, match="boom"):
        pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert len(list(tmp_path.iterdir())) == 1  # kept for inspection


def test_main_cleans_up_tmpdir_on_success(tmp_path, monkeypatch):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _DummyPipeline()
    pipeline.dry_run = True

    pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert list(tmp_path.iterdir()) == []


class _StubExecutor:
    def __init__(self, jobs_with_errors=None):
        self.jobs_with_errors = jobs_with_errors or []


def test_check_execution_names_the_failed_jobs():
    pipeline = _DummyPipeline()
    pipeline.setup_logging(argparse.Namespace(loglevel="WARNING"))
    job = Job(Pipeline(Command("false")), "broken-step")

    with pytest.raises(DagExecutionError, match="broken-step"):
        pipeline.check_execution(DAG(), _StubExecutor([job]))


def test_check_execution_flags_unexecuted_jobs():
    pipeline = _DummyPipeline()
    pipeline.setup_logging(argparse.Namespace(loglevel="WARNING"))
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("echo", "x")), "unexecuted"))

    with pytest.raises(DagExecutionError, match="unexecuted"):
        pipeline.check_execution(dag, _StubExecutor())
