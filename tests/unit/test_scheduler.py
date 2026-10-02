"""
Unit tests for the ThreadScheduler.
"""

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
from sentieon_cli.scheduler import ThreadScheduler  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402


def _job(name, threads=1, resources=None):
    return Job(
        Pipeline(Command(name)), name, threads, resources, task_name="test"
    )


def test_schedules_independent_jobs_together():
    dag = DAG()
    j1 = _job("j1", threads=1)
    j2 = _job("j2", threads=1)
    dag.add_job(j1)
    dag.add_job(j2)

    scheduler = ThreadScheduler(dag, threads=4)
    assert scheduler.start() == {j1, j2}


def test_thread_budget_gates_scheduling():
    dag = DAG()
    j1 = _job("j1", threads=2)
    j2 = _job("j2", threads=2)
    dag.add_job(j1)
    dag.add_job(j2)

    # A 2-thread budget only fits one 2-thread job at a time.
    scheduler = ThreadScheduler(dag, threads=2)

    first = scheduler.start()
    assert len(first) == 1
    running = next(iter(first))

    # The second job starts only once the first frees its threads.
    second = scheduler.job_finished(running)
    assert len(second) == 1
    remaining = next(iter(second))

    assert {running, remaining} == {j1, j2}


def test_named_resource_gates_scheduling():
    dag = DAG()
    j1 = _job("j1", threads=1, resources={"node0": 1})
    j2 = _job("j2", threads=1, resources={"node0": 1})
    dag.add_job(j1)
    dag.add_job(j2)

    # Plenty of threads, but only a single node0 token.
    scheduler = ThreadScheduler(dag, threads=8, resources={"node0": 1})

    first = scheduler.start()
    assert len(first) == 1
    running = next(iter(first))

    second = scheduler.job_finished(running)
    assert len(second) == 1
    assert {running, next(iter(second))} == {j1, j2}


def test_resource_requirements_ignored_when_unmanaged():
    # Jobs may request a resource the scheduler does not track; those
    # requirements are simply ignored.
    dag = DAG()
    j1 = _job("j1", threads=1, resources={"gpu": 1})
    j2 = _job("j2", threads=1, resources={"gpu": 1})
    dag.add_job(j1)
    dag.add_job(j2)

    scheduler = ThreadScheduler(dag, threads=8)
    assert scheduler.start() == {j1, j2}


def test_oversized_thread_request_is_rejected():
    dag = DAG()
    dag.add_job(
        Job(Pipeline(Command("echo", "hi")), "big", 100, task_name="test")
    )
    with pytest.raises(DagExecutionError):
        ThreadScheduler(dag, threads=2).start()


def test_oversized_resource_request_is_rejected():
    dag = DAG()
    job = Job(
        Pipeline(Command("echo", "hi")),
        "greedy",
        1,
        {"node0": 3},
        task_name="test",
    )
    dag.add_job(job)
    with pytest.raises(DagExecutionError):
        ThreadScheduler(dag, threads=8, resources={"node0": 1}).start()
