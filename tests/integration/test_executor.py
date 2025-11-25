"""
Integration tests for the Executor
"""

import io
import os
import pathlib
import sys
import tempfile
from unittest.mock import patch

import pytest

# Add the parent directory to the path
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.executor import LocalExecutor  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.scheduler import ThreadScheduler  # noqa: E402
from sentieon_cli.shell_pipeline import (  # noqa: E402
    Command,
    InputProcSub,
    Pipeline,
)


def test_local_executor_simple_job():
    """Test LocalExecutor with a single, simple job"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        cmd_out = tmp_dir / "test_out.txt"
        dag = DAG()
        job = Job(Pipeline(Command("echo", "hello executor"), file_output=cmd_out), "echo-job")
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)

        executor.execute()

        cmd_out_str = open(cmd_out).read()
        assert "hello executor" in cmd_out_str
        assert len(executor.jobs_with_errors) == 0


def test_local_executor_pipeline_job():
    """Test LocalExecutor with a job containing a pipeline"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        cmd_out = tmp_dir / "test_out.txt"
        dag = DAG()
        pipeline = Pipeline(Command("echo", "hello pipeline"), Command("cat"), file_output=cmd_out)
        job = Job(pipeline, "pipeline-job")
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)

        executor.execute()

        cmd_out_str = open(cmd_out).read()
        assert "hello pipeline" in cmd_out_str
        assert len(executor.jobs_with_errors) == 0

def test_local_executor_failing_job():
    """Test LocalExecutor with a job that fails"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        cmd_in = tmp_dir / "test_in.txt"
        dag = DAG()
        # This command will fail
        job = Job(Pipeline(Command("cat", str(cmd_in))), "failing-job")
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)

        executor.execute()

        assert len(executor.jobs_with_errors) == 1
        assert executor.jobs_with_errors[0] == job

def test_local_executor_proc_sub_job():
    """Test LocalExecutor with a job using process substitution"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        cmd_out = tmp_dir / "test_out.txt"
        dag = DAG()
        command = Command(
            "diff",
            InputProcSub(Pipeline(Command("echo", "a"))),
            InputProcSub(Pipeline(Command("echo", "b"))),
        )
        job = Job(Pipeline(command, file_output=cmd_out), "proc-sub-job")
        dag.add_job(job)

        # Need at least 2 threads for the two proc subs
        scheduler = ThreadScheduler(dag, 2)
        executor = LocalExecutor(scheduler)

        executor.execute()

        # diff will write to stdout and fail
        output = open(cmd_out).read()
        assert "1c1" in output
        assert "< a" in output
        assert "> b" in output

        # The job should be marked as an error because the diff command fails
        assert len(executor.jobs_with_errors) == 1


if __name__ == "__main__":
    pytest.main([__file__])