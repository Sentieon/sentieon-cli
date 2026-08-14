"""
Integration tests for capturing job output to per-process log files.
"""

import asyncio
import logging
import os
import sys
from typing import List

import pytest

# Add the parent directory to the path
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.executor import LocalExecutor  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.run_logs import RunLogs  # noqa: E402
from sentieon_cli.scheduler import ThreadScheduler  # noqa: E402
from sentieon_cli.shell_pipeline import (  # noqa: E402
    Command,
    Context,
    InputProcSub,
    OutputProcSub,
    Pipeline,
)


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


def _run(dag: DAG, log_dir, cores: int = 2):
    """Execute a DAG with its output captured under ``log_dir``"""
    run_logs = RunLogs(log_dir)
    run_logs.create_dirs()
    executor = LocalExecutor(ThreadScheduler(dag, cores), run_logs=run_logs)
    executor.execute()
    return executor, run_logs


def _task_dir(run_logs: RunLogs, task_name: str):
    return run_logs.task_logs / task_name


def test_each_pipeline_stage_gets_its_own_stderr_log(tmp_path):
    """Every stage's stderr is captured; only uncaptured stdout is."""
    dag = DAG()
    job = Job(
        Pipeline(
            Command("sh", "-c", "echo payload; echo first >&2"),
            Command("sh", "-c", "cat; echo second >&2"),
        ),
        "staged",
        task_name="two-stage",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert executor.jobs_with_errors == []
    task_dir = _task_dir(run_logs, "two-stage")
    assert sorted(p.name for p in task_dir.iterdir()) == [
        "staged-1.0.log",
        "staged-1.1.log",
        "staged-1.1.stdout.log",
    ]
    assert "first" in (task_dir / "staged-1.0.log").read_text()
    assert "second" in (task_dir / "staged-1.1.log").read_text()
    # The last stage's stdout was inherited, so it is captured verbatim
    # (sharing the stage index of that stage's stderr log); the stdout
    # piped between the stages is not.
    assert (task_dir / "staged-1.1.stdout.log").read_text() == "payload\n"


def test_a_file_output_wins_over_the_sink(tmp_path):
    """An explicit destination is never redirected into the log dir."""
    out = tmp_path / "out.txt"
    dag = DAG()
    job = Job(
        Pipeline(Command("echo", "hello"), file_output=out),
        "redirected",
        task_name="file-out",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert executor.jobs_with_errors == []
    assert out.read_text() == "hello\n"
    # The command was silent, so its header-only stderr log was pruned.
    assert list(_task_dir(run_logs, "file-out").iterdir()) == []


def test_input_proc_sub_stderr_is_captured(tmp_path):
    """An inner <(...) command writes to the job's sink, not the terminal."""
    out = tmp_path / "out.txt"
    dag = DAG()
    job = Job(
        Pipeline(
            Command(
                "cat",
                InputProcSub(
                    Pipeline(
                        Command("sh", "-c", "echo inner >&2; echo payload")
                    )
                ),
            ),
            file_output=out,
        ),
        "reader",
        task_name="proc-sub",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert executor.jobs_with_errors == []
    assert out.read_text() == "payload\n"
    logs = list(_task_dir(run_logs, "proc-sub").iterdir())
    assert [p.name for p in logs] == ["reader-1.1.log"]
    assert "inner" in logs[0].read_text()


def test_output_proc_sub_stderr_is_captured(tmp_path):
    """An inner >(...) command writes to the job's sink too."""
    dag = DAG()
    job = Job(
        Pipeline(
            Command(
                "sh",
                "-c",
                'echo payload > "$1"',
                "sh",
                OutputProcSub(
                    Pipeline(Command("sh", "-c", "cat >&2; echo done >&2"))
                ),
            )
        ),
        "writer",
        task_name="proc-sub",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert executor.jobs_with_errors == []
    contents = [
        path.read_text() for path in _task_dir(run_logs, "proc-sub").iterdir()
    ]
    assert any("payload" in text and "done" in text for text in contents)


def test_a_successful_job_prunes_its_empty_logs(tmp_path):
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("true")), "quiet", task_name="cleanup"))
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert executor.jobs_with_errors == []
    assert list(_task_dir(run_logs, "cleanup").iterdir()) == []


def test_a_failing_job_is_reported_with_its_log_and_tail(tmp_path, messages):
    dag = DAG()
    job = Job(
        Pipeline(Command("sh", "-c", "echo boom >&2; exit 3")),
        "boomer",
        task_name="failing",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert job in executor.jobs_with_errors
    log_path = _task_dir(run_logs, "failing") / "boomer-1.0.log"
    assert log_path.is_file()  # nothing is pruned after a failure
    assert "boom" in log_path.read_text()

    report = "\n".join(messages)
    pid = job.shell.nodes[0].proc.pid
    assert f"Failed sub-command of {job}" in report
    assert f"exit code: 3, pid: {pid}" in report
    assert str(log_path) in report
    assert "    boom" in report
    assert f"Task logs for debugging: {run_logs.task_logs}" in report


def test_a_failing_stage_reports_its_own_log(tmp_path, messages):
    """The tail comes from the failing process, not a sibling stage."""
    dag = DAG()
    job = Job(
        Pipeline(
            Command("sh", "-c", "echo upstream-noise >&2"),
            Command("sh", "-c", "cat; echo downstream-boom >&2; exit 4"),
        ),
        "mixed",
        task_name="failing",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert job in executor.jobs_with_errors
    report = "\n".join(messages)
    # Quoted log lines are indented; the sibling's output is not quoted.
    assert "    downstream-boom" in report
    assert "    upstream-noise" not in report
    assert str(_task_dir(run_logs, "failing") / "mixed-1.1.log") in report


def test_a_job_that_cannot_start_points_at_its_partial_logs(
    tmp_path, messages
):
    dag = DAG()
    job = Job(
        Pipeline(
            Command("sh", "-c", "echo started >&2; sleep 5"),
            Command("no_such_cmd_zzz"),
        ),
        "unstartable",
        task_name="launch-failure",
    )
    dag.add_job(job)
    executor, run_logs = _run(dag, tmp_path / "logs")

    assert job in executor.jobs_with_errors
    report = "\n".join(messages)
    # The stage that did spawn keeps its log, and the report names it.
    log_path = _task_dir(run_logs, "launch-failure") / "unstartable-1.0.log"
    assert log_path.is_file()
    assert f"Logs from {job}" in report
    assert str(log_path) in report


@pytest.mark.asyncio
async def test_an_interrupted_job_keeps_and_reports_its_logs(
    tmp_path, messages
):
    """Shutdown after an interrupt prunes nothing and names the logs."""
    run_logs = RunLogs(tmp_path / "logs")
    run_logs.create_dirs()
    job = Job(
        Pipeline(Command("sleep", "30")),
        "sleeper",
        task_name="interrupted",
    )
    sink = run_logs.job_sink(job)
    context = Context(log_sink=sink)
    context.file_handles.append(sink.open_stderr(Command("sleep", "30")))

    async def _done() -> int:
        return 0

    executor = LocalExecutor(ThreadScheduler(DAG(), 1), run_logs=run_logs)
    executor.running = [(job, context, asyncio.create_task(_done()), 0)]
    await executor._shutdown()

    log_path = sink.log_paths()[0]
    assert log_path.is_file()
    assert str(log_path) in "\n".join(messages)


def test_without_run_logs_no_files_are_written(tmp_path):
    """An executor with no log directory behaves exactly as before."""
    out = tmp_path / "out.txt"
    dag = DAG()
    dag.add_job(
        Job(
            Pipeline(Command("echo", "hi"), file_output=out),
            "plain",
            task_name="no-logs",
        )
    )
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()

    assert executor.jobs_with_errors == []
    assert [p.name for p in tmp_path.iterdir()] == ["out.txt"]
