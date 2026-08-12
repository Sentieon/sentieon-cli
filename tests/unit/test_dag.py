"""
Unit tests for DAG traversal, duplicate detection, and the resume prune.
"""

import os
import sys

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.dag import DAG, has_all_outputs  # noqa: E402
from sentieon_cli.exceptions import DagExecutionError  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402


def _job(name, threads=1, resources=None):
    """A trivial job whose identity is its (unique) command name."""
    return Job(Pipeline(Command(name)), name, threads, resources)


def test_job_identity_is_the_pipeline():
    # Two jobs with the same pipeline are equal even if named
    # differently; a DAG deduplicates on this identity.
    j1 = Job(Pipeline(Command("echo", "x")), "name-a")
    j2 = Job(Pipeline(Command("echo", "x")), "name-b")
    j3 = Job(Pipeline(Command("echo", "y")), "name-c")
    assert j1 == j2
    assert hash(j1) == hash(j2)
    assert j1 != j3


class TestDAGTraversal:
    """Test DAG bookkeeping through mark_finished"""

    def test_add_job_with_unknown_dependency_raises(self):
        dag = DAG()
        known = _job("known")
        dag.add_job(known)

        unknown = _job("unknown")
        with pytest.raises(ValueError):
            dag.add_job(_job("downstream"), {unknown})

    def test_finishing_an_unready_job_raises(self):
        dag = DAG()
        a = _job("a")
        b = _job("b")
        dag.add_job(a)
        dag.add_job(b, {a})

        with pytest.raises(DagExecutionError):
            dag.mark_finished(b)  # b is not ready

    def test_diamond_dependencies(self):
        dag = DAG()
        top = _job("top")
        left = _job("left")
        right = _job("right")
        bottom = _job("bottom")

        dag.add_job(top)
        dag.add_job(left, {top})
        dag.add_job(right, {top})
        dag.add_job(bottom, {left, right})

        assert set(dag.ready_jobs) == {top}

        assert set(dag.mark_finished(top)) == {left, right}

        # bottom needs BOTH left and right
        assert set(dag.mark_finished(left)) == set()
        assert set(dag.mark_finished(right)) == {bottom}

    def test_multilevel_chain(self):
        dag = DAG()
        a = _job("a")
        b = _job("b")
        c = _job("c")

        dag.add_job(a)
        dag.add_job(b, {a})
        dag.add_job(c, {b})

        assert set(dag.ready_jobs) == {a}
        assert set(dag.mark_finished(a)) == {b}
        assert set(dag.mark_finished(b)) == {c}

        # Everything has been consumed and marked finished.
        dag.mark_finished(c)
        assert a in dag.finished_jobs
        assert b in dag.finished_jobs
        assert c in dag.finished_jobs
        assert len(dag.waiting_jobs) == 0
        assert len(dag.ready_jobs) == 0

    def test_add_duplicate_pipeline_raises(self):
        dag = DAG()
        dag.add_job(Job(Pipeline(Command("echo", "x")), "step-1"))
        # A second job with an identical pipeline collides on identity.
        with pytest.raises(ValueError):
            dag.add_job(Job(Pipeline(Command("echo", "x")), "step-2"))

    def test_jobs_differing_only_in_exec_kwargs_are_not_duplicates(self):
        # A command's env/cwd is part of its identity, so two jobs that
        # differ only there are distinct, not rejected as duplicates.
        dag = DAG()
        j1 = Job(Pipeline(Command("run", exec_kwargs={"cwd": "/a"})), "a")
        j2 = Job(Pipeline(Command("run", exec_kwargs={"cwd": "/b"})), "b")
        dag.add_job(j1)
        dag.add_job(j2)  # must not raise "already in the DAG"
        assert j1 in dag.ready_jobs
        assert j2 in dag.ready_jobs

    def test_add_job_accepts_list_dependencies(self):
        # Regression: non-set dependencies used to be silently dropped,
        # leaving the job immediately ready.
        dag = DAG()
        a = _job("a")
        b = _job("b")
        dag.add_job(a)
        dag.add_job(b, [a])

        assert b not in dag.ready_jobs
        assert dag.waiting_jobs[b] == {a}
        assert set(dag.mark_finished(a)) == {b}

    def test_add_job_accepts_any_iterable_dependencies(self):
        # Tuples, frozensets, and single-use generators all work; the
        # generator case also guards against consuming the iterable twice.
        for make in (tuple, frozenset, iter):
            dag = DAG()
            a = _job("a")
            b = _job("b")
            dag.add_job(a)
            dag.add_job(b, make([a]))
            assert dag.waiting_jobs[b] == {a}, make

    def test_add_job_rejects_unknown_dependency_in_list(self):
        dag = DAG()
        dag.add_job(_job("known"))
        with pytest.raises(ValueError):
            dag.add_job(_job("downstream"), [_job("unknown")])

    def test_add_job_with_empty_iterable_is_ready(self):
        for deps in ([], iter([])):
            dag = DAG()
            job = _job("solo")
            dag.add_job(job, deps)
            assert job in dag.ready_jobs

    def test_stored_dependencies_are_independent_of_caller_set(self):
        dag = DAG()
        a = _job("a")
        b = _job("b")
        c = _job("c")
        dag.add_job(a)
        dag.add_job(b)

        deps = {a}
        dag.add_job(c, deps)
        deps.add(b)  # mutating the caller's set must not affect the DAG

        assert dag.waiting_jobs[c] == {a}


class TestSkipSatisfied:
    """Test the pre-run resume prune"""

    def test_satisfied_root_is_marked_finished(self):
        dag = DAG()
        job = _job("a")
        dag.add_job(job)
        skipped = dag.skip_satisfied(lambda j: True)
        assert skipped == [job]
        assert job in dag.finished_jobs
        assert len(dag.ready_jobs) == 0

    def test_unsatisfied_root_stays_ready(self):
        dag = DAG()
        job = _job("a")
        dag.add_job(job)
        skipped = dag.skip_satisfied(lambda j: False)
        assert skipped == []
        assert job in dag.ready_jobs
        assert dag.finished_jobs == []

    def test_cascades_down_a_chain_in_topological_order(self):
        dag = DAG()
        a, b, c = _job("a"), _job("b"), _job("c")
        dag.add_job(a)
        dag.add_job(b, {a})
        dag.add_job(c, {b})
        skipped = dag.skip_satisfied(lambda j: True)
        assert skipped == [a, b, c]
        assert dag.finished_jobs == [a, b, c]
        assert len(dag.ready_jobs) == 0
        assert len(dag.waiting_jobs) == 0

    def test_stops_at_the_first_unsatisfied_job(self):
        dag = DAG()
        a, b, c = _job("a"), _job("b"), _job("c")
        dag.add_job(a)
        dag.add_job(b, {a})
        dag.add_job(c, {b})
        skipped = dag.skip_satisfied(lambda j: j == a)
        assert skipped == [a]
        assert b in dag.ready_jobs
        assert c in dag.waiting_jobs

    def test_blocked_jobs_are_not_offered_to_the_predicate(self):
        dag = DAG()
        a, b, c = _job("a"), _job("b"), _job("c")
        dag.add_job(a)
        dag.add_job(b, {a})
        dag.add_job(c, {b})
        seen = []

        def predicate(job):
            seen.append(job)
            return job == a

        dag.skip_satisfied(predicate)
        # b is offered (a was skipped, unblocking it); c never is -- the
        # ancestor rule: a job behind an unsatisfied job cannot be skipped.
        assert seen == [a, b]

    def test_diamond_with_partial_parents(self):
        dag = DAG()
        top = _job("top")
        left, right = _job("left"), _job("right")
        bottom = _job("bottom")
        dag.add_job(top)
        dag.add_job(left, {top})
        dag.add_job(right, {top})
        dag.add_job(bottom, {left, right})
        skipped = dag.skip_satisfied(lambda j: j in (top, left))
        assert skipped == [top, left]
        assert right in dag.ready_jobs
        assert dag.waiting_jobs[bottom] == {right}

    def test_empty_dag_returns_empty(self):
        assert DAG().skip_satisfied(lambda j: True) == []


class TestHasAllOutputs:
    """Test the default resume predicate"""

    def test_no_declared_outputs_is_never_satisfied(self):
        assert has_all_outputs(_job("a")) is False

    def test_true_when_all_outputs_exist(self, tmp_path):
        full = tmp_path / "full.txt"
        full.write_text("data")
        empty = tmp_path / "empty.txt"
        empty.touch()  # zero-byte files count as present
        job = Job(Pipeline(Command("true")), "j", outputs=[full, empty])
        assert has_all_outputs(job) is True

    def test_false_when_any_output_is_missing(self, tmp_path):
        present = tmp_path / "present.txt"
        present.touch()
        missing = tmp_path / "missing.txt"
        job = Job(Pipeline(Command("true")), "j", outputs=[present, missing])
        assert has_all_outputs(job) is False

    def test_dangling_symlink_counts_as_missing(self, tmp_path):
        target = tmp_path / "target.txt"
        target.touch()
        link = tmp_path / "link.txt"
        link.symlink_to(target)
        target.unlink()  # now dangling
        job = Job(Pipeline(Command("true")), "j", outputs=[link])
        assert has_all_outputs(job) is False
