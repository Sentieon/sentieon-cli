"""
Unit tests for the machine-readable metrics files
"""

import os
import pathlib
import sys
from typing import List

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.executor import LocalExecutor  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.run_logs import (  # noqa: E402
    JOB_METRICS_COLUMNS,
    PROCESS_METRICS_COLUMNS,
    RunLogs,
)
from sentieon_cli.scheduler import ThreadScheduler  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402


def _run(log_dir: pathlib.Path, pipeline: Pipeline) -> RunLogs:
    """Run a single job with its metrics collected under ``log_dir``"""
    dag = DAG()
    dag.add_job(Job(pipeline, "metrics", task_name="metrics"))
    run_logs = RunLogs(log_dir)
    run_logs.create_dirs()
    LocalExecutor(ThreadScheduler(dag, 1), run_logs=run_logs).execute()
    return run_logs


def _rows(path: pathlib.Path) -> List[List[str]]:
    """The rows of a metrics file, header included, split on tabs"""
    lines = path.read_text().splitlines()
    return [line.split("\t") for line in lines]


def _row_dict(path: pathlib.Path, index: int = 1) -> dict:
    rows = _rows(path)
    return dict(zip(rows[0], rows[index]))


def test_process_metrics_has_a_row_per_process(tmp_path):
    run_logs = _run(tmp_path / "logs", Pipeline(Command("true")))

    rows = _rows(run_logs.process_metrics)
    assert rows[0] == PROCESS_METRICS_COLUMNS
    assert len(rows) == 2
    assert len(rows[1]) == len(PROCESS_METRICS_COLUMNS)
    row = dict(zip(rows[0], rows[1]))
    assert row["job_id"] == "metrics-1"
    assert row["task_name"] == "metrics"
    assert row["exit_code"] == "0"
    assert int(row["pid"]) > 0
    assert float(row["user_s"]) >= 0
    assert float(row["maxrss_mib"]) > 0
    assert row["command"] == "true"


def test_job_metrics_has_a_row_per_job(tmp_path):
    # `wall_s` is rendered with two decimals, so the job has to take long
    # enough to be visible there.
    run_logs = _run(tmp_path / "logs", Pipeline(Command("sleep", "0.05")))

    rows = _rows(run_logs.job_metrics)
    assert rows[0] == JOB_METRICS_COLUMNS
    assert len(rows) == 2
    row = dict(zip(rows[0], rows[1]))
    assert row["job_id"] == "metrics-1"
    assert row["status"] == "ok"
    assert float(row["wall_s"]) > 0
    assert row["processes"] == "1"
    assert row["threads"] == "1"


def test_process_rows_link_to_job_rows(tmp_path):
    run_logs = _run(
        tmp_path / "logs",
        Pipeline(Command("printf", "x"), Command("cat")),
    )

    process_rows = _rows(run_logs.process_metrics)[1:]
    assert len(process_rows) == 2
    header = PROCESS_METRICS_COLUMNS
    pids = {row[header.index("pid")] for row in process_rows}
    assert len(pids) == 2
    job_row = _row_dict(run_logs.job_metrics)
    assert {row[header.index("job_id")] for row in process_rows} == {
        job_row["job_id"]
    }
    assert job_row["processes"] == "2"


def test_stage_links_rows_to_task_log_files(tmp_path):
    run_logs = _run(tmp_path / "logs", Pipeline(Command("true")))

    # The log file itself may be gone: a successful job prunes logs that
    # hold nothing but their header.
    assert _row_dict(run_logs.process_metrics)["stage"] == "0"


def test_failed_job_still_records_metrics(tmp_path):
    run_logs = _run(tmp_path / "logs", Pipeline(Command("false")))

    process_row = _row_dict(run_logs.process_metrics)
    assert process_row["exit_code"] == "1"
    assert process_row["user_s"] != ""
    assert process_row["maxrss_mib"] != ""
    assert _row_dict(run_logs.job_metrics)["status"] == "failed"


def test_command_field_keeps_the_row_rectangular(tmp_path):
    run_logs = _run(tmp_path / "logs", Pipeline(Command("printf", "a\tb\nc")))

    rows = _rows(run_logs.process_metrics)
    assert len(rows) == 2
    assert len(rows[1]) == len(PROCESS_METRICS_COLUMNS)
    command = rows[1][PROCESS_METRICS_COLUMNS.index("command")]
    assert "\t" not in command
    assert "\n" not in command


def test_rerun_truncates_metrics_files(tmp_path):
    log_dir = tmp_path / "logs"
    _run(log_dir, Pipeline(Command("false")))

    Job.reset_ids()
    run_logs = _run(log_dir, Pipeline(Command("true")))

    assert len(_rows(run_logs.process_metrics)) == 2
    assert _row_dict(run_logs.process_metrics)["exit_code"] == "0"
    assert len(_rows(run_logs.job_metrics)) == 2
    assert _row_dict(run_logs.job_metrics)["status"] == "ok"


def test_metrics_files_start_as_headers_only(tmp_path):
    run_logs = RunLogs(tmp_path / "logs")
    run_logs.create_dirs()

    assert _rows(run_logs.process_metrics) == [PROCESS_METRICS_COLUMNS]
    assert _rows(run_logs.job_metrics) == [JOB_METRICS_COLUMNS]
