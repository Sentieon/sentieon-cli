"""
Unit tests for Job identity: task_name, job_id, and id numbering.
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
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402
from sentieon_cli.util import sanitize  # noqa: E402


def _job(name, arg, task_name="test"):
    """A job whose identity is its (unique) command argument."""
    return Job(Pipeline(Command("echo", arg)), name, task_name=task_name)


def test_task_name_is_required():
    with pytest.raises(TypeError):
        Job(  # type: ignore[call-arg]
            Pipeline(Command("echo", "x")), "no-task"
        )


def test_ids_are_numbered_per_name():
    jobs = [_job("shard", str(i)) for i in range(3)]
    assert [job.job_id for job in jobs] == ["shard-1", "shard-2", "shard-3"]


def test_each_name_has_its_own_counter():
    first = _job("dedup", "a")
    other = _job("metrics", "b")
    second = _job("dedup", "c")

    assert first.job_id == "dedup-1"
    assert other.job_id == "metrics-1"
    assert second.job_id == "dedup-2"


def test_names_differing_only_in_unsafe_characters_share_a_counter():
    # Log file names are sanitized, so ids that only differ in unsafe
    # characters would name the same file and truncate each other.
    first = _job("a b", "x")
    second = _job("a-b", "y")

    assert (first.job_id, second.job_id) == ("a b-1", "a-b-2")
    assert sanitize(first.job_id) == "a-b-1"
    assert sanitize(second.job_id) == "a-b-2"


def test_reset_ids_restarts_the_sequence():
    assert _job("multiqc", "a").job_id == "multiqc-1"
    assert _job("multiqc", "b").job_id == "multiqc-2"

    Job.reset_ids()

    assert _job("multiqc", "c").job_id == "multiqc-1"


def test_ids_keep_counting_across_two_dags():
    # The pangenome pipeline executes two DAGs in one process; ids must stay
    # unique for the whole run, so they are never reset between DAGs.
    def build_dag(tag):
        dag = DAG()
        for i in range(2):
            dag.add_job(_job("calling", f"{tag}-{i}"))
        return dag

    first = build_dag("a")
    second = build_dag("b")

    ids = sorted(
        job.job_id
        for dag in (first, second)
        for job in list(dag.ready_jobs) + list(dag.waiting_jobs)
    )
    assert ids == ["calling-1", "calling-2", "calling-3", "calling-4"]


def test_task_name_and_job_id_are_not_part_of_identity():
    j1 = Job(Pipeline(Command("echo", "x")), "a", task_name="alignment")
    j2 = Job(Pipeline(Command("echo", "x")), "b", task_name="cleanup")

    assert j1.job_id != j2.job_id
    assert j1.task_name != j2.task_name
    assert j1 == j2
    assert hash(j1) == hash(j2)


def test_jobs_with_different_task_names_still_collide_in_a_dag():
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("echo", "x")), "a", task_name="dedup"))
    with pytest.raises(ValueError):
        dag.add_job(
            Job(Pipeline(Command("echo", "x")), "b", task_name="metrics")
        )


def test_repr_and_str_report_the_job_id():
    job = _job("locuscollector", "x")
    assert repr(job) == "Job(locuscollector-1)"
    assert str(job) == "Job(locuscollector-1)"


def test_task_name_is_stored():
    job = _job("dnascope", "x", task_name="variant-calling")
    assert job.task_name == "variant-calling"
