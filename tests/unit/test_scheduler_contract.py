"""
Conformance tests for the BaseScheduler contract.

``drive`` mimics what an executor does: ask the scheduler for the initial
batch, then report each job finished (FIFO) and collect whatever that
unblocks, until the DAG drains. Any BaseScheduler implementation can be
dropped into ``drive`` to check that it honors the contract; the
dependency-order tests run against every scheduler in ``SCHEDULERS``.
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

SCHEDULERS = [
    pytest.param(lambda dag: ThreadScheduler(dag, threads=8), id="thread"),
]


def _job(name, threads=1, resources=None):
    return Job(Pipeline(Command(name)), name, threads, resources)


def drive(scheduler):
    """Drive a scheduler to completion; return the job finish order."""
    order = []
    queue = list(scheduler.start())
    while queue:
        job = queue.pop(0)
        order.append(job)
        queue.extend(scheduler.job_finished(job))
    return order


@pytest.mark.parametrize("make_scheduler", SCHEDULERS)
def test_drives_a_diamond_to_completion(make_scheduler):
    dag = DAG()
    top = _job("top")
    left = _job("left")
    right = _job("right")
    bottom = _job("bottom")
    dag.add_job(top)
    dag.add_job(left, {top})
    dag.add_job(right, {top})
    dag.add_job(bottom, {left, right})

    order = drive(make_scheduler(dag))

    # Every job runs exactly once.
    assert len(order) == len(set(order))
    assert sorted(j.name for j in order) == [
        "bottom",
        "left",
        "right",
        "top",
    ]
    # Dependencies precede their dependents.
    pos = {job: i for i, job in enumerate(order)}
    assert pos[top] < pos[left] < pos[bottom]
    assert pos[top] < pos[right] < pos[bottom]


@pytest.mark.parametrize("make_scheduler", SCHEDULERS)
def test_drives_a_chain_to_completion(make_scheduler):
    dag = DAG()
    a = _job("a")
    b = _job("b")
    c = _job("c")
    dag.add_job(a)
    dag.add_job(b, {a})
    dag.add_job(c, {b})

    order = drive(make_scheduler(dag))
    assert [j.name for j in order] == ["a", "b", "c"]


@pytest.mark.parametrize("make_scheduler", SCHEDULERS)
def test_job_finished_on_unknown_job_raises(make_scheduler):
    dag = DAG()
    dag.add_job(_job("a"))
    scheduler = make_scheduler(dag)
    scheduler.start()

    with pytest.raises(DagExecutionError):
        scheduler.job_finished(_job("never-scheduled"))


def test_budget_constrained_run_still_completes():
    # ThreadScheduler-specific: a 1-thread budget serializes two independent
    # jobs, but driving to completion still runs both exactly once.
    dag = DAG()
    dag.add_job(_job("a", threads=1))
    dag.add_job(_job("b", threads=1))

    order = drive(ThreadScheduler(dag, threads=1))
    assert sorted(j.name for j in order) == ["a", "b"]
