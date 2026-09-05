"""
Unit tests for the BasePipeline lifecycle (argument registration, logging
scope, temp-dir cleanup, and post-run checks).
"""

import argparse
import logging
import os
import sys
from typing import Optional

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.dnascope import DNAscopePipeline  # noqa: E402
from sentieon_cli.dnascope_hybrid import (  # noqa: E402
    DNAscopeHybridPipeline,
)
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline  # noqa: E402
from sentieon_cli.exceptions import DagExecutionError  # noqa: E402
from sentieon_cli.hybrid_pangenome import HybridPangenome  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.pipeline import BasePipeline  # noqa: E402
from sentieon_cli.sentieon_pangenome import SentieonPangenome  # noqa: E402
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


class _TwoDagPipeline(_DummyPipeline):
    """A pipeline that counts the DAGs run through the `main` hook.

    `second_dag` is the DAG returned by the hook; `fail_on` is the 1-based
    index of the `check_execution` call that raises.
    """

    def __init__(self, second_dag=None, fail_on=None):
        super().__init__()
        self.second_dag = second_dag
        self.fail_on = fail_on
        self.dags_run = 0
        self.check_calls = 0

    def build_second_dag(self):
        return self.second_dag

    def run(self, dag):
        self.dags_run += 1
        return _StubExecutor()

    def check_execution(self, dag, executor):
        self.check_calls += 1
        if self.fail_on == self.check_calls:
            raise DagExecutionError("a job failed")


def test_main_runs_one_dag_when_the_second_dag_hook_returns_none(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _TwoDagPipeline()

    pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert pipeline.dags_run == 1
    assert pipeline.check_calls == 1
    assert list(tmp_path.iterdir()) == []


def test_main_runs_the_dag_returned_by_the_second_dag_hook(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _TwoDagPipeline(second_dag=DAG())

    pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert pipeline.dags_run == 2
    assert pipeline.check_calls == 2
    assert list(tmp_path.iterdir()) == []


def test_main_cleans_up_tmpdir_when_the_second_dag_fails(
    tmp_path, monkeypatch
):
    # The second DAG runs inside the same try/finally, so the temp
    # directory is removed on its failure path too
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    pipeline = _TwoDagPipeline(second_dag=DAG(), fail_on=2)

    with pytest.raises(DagExecutionError):
        pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert pipeline.dags_run == 2
    assert list(tmp_path.iterdir()) == []


def test_hybrid_second_dag_hook_is_skipped_without_cnv_calling():
    from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline

    pipeline = DNAscopeHybridPipeline()
    pipeline.setup_logging(argparse.Namespace(loglevel="WARNING"))
    pipeline.skip_cnv = True

    assert pipeline.build_second_dag() is None


def test_check_execution_names_the_failed_jobs():
    pipeline = _DummyPipeline()
    pipeline.setup_logging(argparse.Namespace(loglevel="WARNING"))
    job = Job(Pipeline(Command("false")), "broken-step", task_name="test")

    # Jobs report themselves by job_id, not by name.
    with pytest.raises(DagExecutionError, match=r"Job\(broken-step-1\)"):
        pipeline.check_execution(DAG(), _StubExecutor([job]))


class _StubPangenome(SentieonPangenome):
    """A SentieonPangenome with its run stages stubbed out.

    `fail_on` is the 1-based index of the `check_execution` call that
    raises, mimicking a job failure in the first or the second DAG.
    """

    def __init__(self, fail_on=None):
        super().__init__()
        self.fail_on = fail_on
        self.check_calls = 0
        self.dags_built = []

    def validate(self) -> None:
        pass

    def configure(self) -> None:
        pass

    def build_dag(self) -> DAG:
        self.dags_built.append("first")
        self.ploidy_json = self.tmp_dir.joinpath("ploidy.json")
        return DAG()

    def build_second_dag(self) -> Optional[DAG]:
        if not self._needs_second_dag():
            return None
        self.dags_built.append("second")
        return DAG()

    def get_sex(self, ploidy_json) -> None:
        pass

    def run(self, dag):
        return _StubExecutor()

    def check_execution(self, dag, executor):
        self.check_calls += 1
        if self.fail_on == self.check_calls:
            raise DagExecutionError("a job failed")


def _stub_pangenome(monkeypatch, tmp_path, **kwargs):
    """A stubbed pangenome pipeline writing temp dirs into `tmp_path`"""
    monkeypatch.setenv("SENTIEON_TMPDIR", str(tmp_path))
    return _StubPangenome(**kwargs)


def test_pangenome_main_cleans_up_tmpdir_when_the_first_dag_fails(
    tmp_path, monkeypatch
):
    # SentieonPangenome runs the base two-DAG flow; the temp directory
    # has to be removed on the failure path too.
    pipeline = _stub_pangenome(monkeypatch, tmp_path, fail_on=1)

    with pytest.raises(DagExecutionError):
        pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert list(tmp_path.iterdir()) == []  # no leaked temp dir


def test_pangenome_main_cleans_up_tmpdir_when_the_second_dag_fails(
    tmp_path, monkeypatch
):
    pipeline = _stub_pangenome(monkeypatch, tmp_path, fail_on=2)
    pipeline.segdup_caller = []  # request the second DAG

    with pytest.raises(DagExecutionError):
        pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert pipeline.dags_built == ["first", "second"]
    assert list(tmp_path.iterdir()) == []


def test_pangenome_main_runs_the_second_dag_for_cnv_calling(
    tmp_path, monkeypatch
):
    # CNV calling is sex-aware, so it alone triggers the second DAG
    pipeline = _stub_pangenome(monkeypatch, tmp_path)
    pipeline.call_svs = True
    pipeline.has_cnv_model = True

    pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert pipeline.dags_built == ["first", "second"]
    assert pipeline.check_calls == 2


def test_pangenome_main_retains_tmpdir_on_failure_when_requested(
    tmp_path, monkeypatch
):
    pipeline = _stub_pangenome(monkeypatch, tmp_path, fail_on=1)
    pipeline.retain_tmpdir = True

    with pytest.raises(DagExecutionError):
        pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert len(list(tmp_path.iterdir())) == 1  # kept for inspection


def test_pangenome_main_cleans_up_tmpdir_on_success(tmp_path, monkeypatch):
    pipeline = _stub_pangenome(monkeypatch, tmp_path)
    pipeline.segdup_caller = []  # request the second DAG

    pipeline.main(argparse.Namespace(loglevel="WARNING"))

    assert pipeline.dags_built == ["first", "second"]
    assert pipeline.check_calls == 2
    assert list(tmp_path.iterdir()) == []


def test_check_execution_flags_unexecuted_jobs():
    pipeline = _DummyPipeline()
    pipeline.setup_logging(argparse.Namespace(loglevel="WARNING"))
    dag = DAG()
    dag.add_job(
        Job(Pipeline(Command("echo", "x")), "unexecuted", task_name="test")
    )

    with pytest.raises(DagExecutionError, match=r"Job\(unexecuted-1\)"):
        pipeline.check_execution(dag, _StubExecutor())


@pytest.mark.parametrize(
    "pipeline_cls",
    [
        DNAscopePipeline,
        DNAscopeLRPipeline,
        DNAscopeHybridPipeline,
        SentieonPangenome,
        HybridPangenome,
    ],
)
def test_no_pipeline_overrides_main(pipeline_cls):
    """Every pipeline runs through `BasePipeline.main`"""
    assert "main" not in pipeline_cls.__dict__
