"""
Integration tests for the executors.
"""

import asyncio
import os
import pathlib
import sys
import tempfile
import threading

import pytest

# Add the parent directory to the path
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.dag import DAG, has_all_outputs  # noqa: E402
from sentieon_cli.exceptions import DagExecutionError  # noqa: E402
from sentieon_cli.executor import (  # noqa: E402
    DryRunExecutor,
    LocalExecutor,
    _cleanup_quietly,
)
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
        cmd_out = pathlib.Path(tmp_dir_str) / "test_out.txt"
        dag = DAG()
        job = Job(
            Pipeline(Command("echo", "hello executor"), file_output=cmd_out),
            "echo-job",
        )
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)
        executor.execute()

        assert "hello executor" in cmd_out.read_text()
        assert len(executor.jobs_with_errors) == 0


def test_local_executor_pipeline_job():
    """Test LocalExecutor with a job containing a pipeline"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        cmd_out = pathlib.Path(tmp_dir_str) / "test_out.txt"
        dag = DAG()
        pipeline = Pipeline(
            Command("echo", "hello pipeline"),
            Command("cat"),
            file_output=cmd_out,
        )
        job = Job(pipeline, "pipeline-job")
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)
        executor.execute()

        assert "hello pipeline" in cmd_out.read_text()
        assert len(executor.jobs_with_errors) == 0


def test_local_executor_failing_job():
    """Test LocalExecutor with a job that fails"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        cmd_in = pathlib.Path(tmp_dir_str) / "test_in.txt"
        dag = DAG()
        # This command will fail: the input file does not exist
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
        cmd_out = pathlib.Path(tmp_dir_str) / "test_out.txt"
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
        output = cmd_out.read_text()
        assert "1c1" in output
        assert "< a" in output
        assert "> b" in output

        # The job is marked as an error because the diff command fails
        assert len(executor.jobs_with_errors) == 1


def test_dry_run_executor_prints_without_running(capsys):
    """DryRunExecutor prints the shell but does not execute it"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        cmd_out = pathlib.Path(tmp_dir_str) / "out.txt"
        dag = DAG()
        job = Job(
            Pipeline(Command("echo", "dryrun"), file_output=cmd_out),
            "dry-job",
        )
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = DryRunExecutor(scheduler)
        executor.execute()

        printed = capsys.readouterr().out
        assert "echo" in printed
        assert "dryrun" in printed
        assert not cmd_out.exists()
        assert len(executor.jobs_with_errors) == 0


def test_dry_run_executor_drains_dependencies(capsys):
    """DryRunExecutor walks the whole DAG, including dependent jobs"""
    dag = DAG()
    a = Job(Pipeline(Command("echo", "a")), "a")
    b = Job(Pipeline(Command("echo", "b")), "b")
    dag.add_job(a)
    dag.add_job(b, {a})

    scheduler = ThreadScheduler(dag, 4)
    executor = DryRunExecutor(scheduler)
    executor.execute()

    printed = capsys.readouterr().out
    assert "echo a" in printed
    assert "echo b" in printed


def test_fail_ok_subcommand_does_not_fail_job():
    """A failing but fail_ok sub-command does not mark the job as errored"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        missing = tmp_dir / "missing.txt"
        out = tmp_dir / "out.txt"

        # `cat missing | wc -c`; cat fails (no file) but is tolerated.
        dag = DAG()
        job = Job(
            Pipeline(
                Command("cat", str(missing), fail_ok=True),
                Command("wc", "-c"),
                file_output=out,
            ),
            "tolerant-job",
        )
        dag.add_job(job)

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)
        executor.execute()

        assert len(executor.jobs_with_errors) == 0


def test_independent_jobs_all_complete():
    """Two independent jobs both run and produce output"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        o1 = tmp_dir / "o1.txt"
        o2 = tmp_dir / "o2.txt"
        dag = DAG()
        dag.add_job(Job(Pipeline(Command("echo", "1"), file_output=o1), "j1"))
        dag.add_job(Job(Pipeline(Command("echo", "2"), file_output=o2), "j2"))

        scheduler = ThreadScheduler(dag, 2)
        executor = LocalExecutor(scheduler)
        executor.execute()

        assert o1.read_text().strip() == "1"
        assert o2.read_text().strip() == "2"
        assert len(executor.jobs_with_errors) == 0


def test_dependent_job_runs_after_dependency():
    """A dependent job runs and both jobs complete successfully"""
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = pathlib.Path(tmp_dir_str)
        first = tmp_dir / "first.txt"
        second = tmp_dir / "second.txt"
        dag = DAG()
        a = Job(Pipeline(Command("echo", "a"), file_output=first), "a")
        b = Job(Pipeline(Command("echo", "b"), file_output=second), "b")
        dag.add_job(a)
        dag.add_job(b, {a})

        scheduler = ThreadScheduler(dag, 1)
        executor = LocalExecutor(scheduler)
        executor.execute()

        assert first.read_text().strip() == "a"
        assert second.read_text().strip() == "b"
        assert len(executor.jobs_with_errors) == 0


def test_job_that_fails_to_start_is_recorded():
    """A job whose command cannot be spawned is recorded as an error."""
    dag = DAG()
    job = Job(Pipeline(Command("no_such_cmd_zzz")), "ghost")
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    executor.execute()
    assert job in executor.jobs_with_errors


def test_launch_loop_stops_after_a_failed_launch():
    """Once a job in a ready batch fails to launch, the loop stops starting
    the rest of the batch instead of launching all of it."""
    dag = DAG()
    for i in range(5):
        # Distinct args -> distinct pipeline identities (else the DAG
        # rejects them as duplicates).
        dag.add_job(Job(Pipeline(Command("true", str(i))), f"j{i}"))
    executor = LocalExecutor(ThreadScheduler(dag, 5))

    launched = []

    async def fake_run_job(job):
        launched.append(job)
        executor.start_new_jobs = False  # simulate a launch failure

    executor.run_job = fake_run_job
    executor.execute()

    # The guard stopped the loop after the first job rather than starting
    # all five.
    assert len(launched) == 1


def test_infeasible_job_raises_instead_of_deadlocking():
    """A job larger than the budget raises rather than stalling silently."""
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("echo", "hi")), "big", 100))
    with pytest.raises(DagExecutionError):
        LocalExecutor(ThreadScheduler(dag, 2)).execute()


def test_dry_run_also_rejects_infeasible_job():
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("echo", "hi")), "big", 100))
    with pytest.raises(DagExecutionError):
        DryRunExecutor(ThreadScheduler(dag, 2)).execute()


def test_non_oserror_on_start_propagates(monkeypatch):
    """A non-OSError while starting a job is a bug and must propagate."""

    async def boom(self, *args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr("sentieon_cli.shell_pipeline.Pipeline.run", boom)
    dag = DAG()
    dag.add_job(Job(Pipeline(Command("echo", "x")), "j"))
    with pytest.raises(RuntimeError):
        LocalExecutor(ThreadScheduler(dag, 1)).execute()


def test_missing_input_fails_before_run(tmp_path):
    """A declared input that is missing fails the job before it runs."""
    missing = tmp_path / "in.txt"
    marker = tmp_path / "ran.txt"
    dag = DAG()
    job = Job(
        Pipeline(Command("touch", str(marker))),
        "j",
        inputs=[missing],
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors
    assert not marker.exists()  # the command never ran


def test_present_inputs_and_outputs_allow_success(tmp_path):
    """Empty (zero-byte) declared inputs/outputs still count as present."""
    inp = tmp_path / "in.txt"
    inp.touch()  # empty input file
    out = tmp_path / "out.txt"
    dag = DAG()
    job = Job(
        Pipeline(Command("touch", str(out))),  # produces an empty output
        "j",
        inputs=[inp],
        outputs=[out],
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert executor.jobs_with_errors == []
    assert out.exists()


def test_missing_output_fails_after_run(tmp_path):
    """A declared output the command does not produce fails the job."""
    missing = tmp_path / "never.txt"
    dag = DAG()
    job = Job(Pipeline(Command("true")), "j", outputs=[missing])
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors


def test_failed_job_outputs_are_removed(tmp_path):
    """A failed job's declared outputs are deleted from disk."""
    out = tmp_path / "out.txt"
    out.write_text("stale")
    dag = DAG()
    job = Job(Pipeline(Command("false")), "j", outputs=[out])
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors
    assert not out.exists()


def test_missing_output_failure_removes_produced_siblings(tmp_path):
    """A missing declared output deletes the outputs that were produced."""
    produced = tmp_path / "out1.txt"
    never = tmp_path / "out2.txt"
    dag = DAG()
    job = Job(
        Pipeline(Command("touch", str(produced))),
        "j",
        outputs=[produced, never],
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors
    assert not produced.exists()


def test_launch_failure_removes_outputs(tmp_path):
    """A job that fails to launch has its declared outputs deleted."""
    out = tmp_path / "out.txt"
    out.write_text("stale")
    dag = DAG()
    job = Job(Pipeline(Command("no_such_cmd_zzz")), "j", outputs=[out])
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors
    assert not out.exists()


def test_missing_input_failure_removes_outputs(tmp_path):
    """Any failed attempt clears the job's declared outputs -- even a
    missing-input failure where the command never ran."""
    missing = tmp_path / "in.txt"
    out = tmp_path / "out.txt"
    out.write_text("stale")
    dag = DAG()
    job = Job(Pipeline(Command("true")), "j", inputs=[missing], outputs=[out])
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors
    assert not out.exists()


def test_successful_job_keeps_outputs(tmp_path):
    """A successful job's outputs are never deleted."""
    out = tmp_path / "out.txt"
    dag = DAG()
    job = Job(Pipeline(Command("touch", str(out))), "j", outputs=[out])
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert executor.jobs_with_errors == []
    assert out.exists()


def test_failed_job_directory_output_is_left_in_place(tmp_path):
    """Directories declared as outputs are never removed on failure."""
    out_dir = tmp_path / "out_dir"
    out_dir.mkdir()
    (out_dir / "keep.txt").write_text("keep")
    dag = DAG()
    job = Job(Pipeline(Command("false")), "j", outputs=[out_dir])
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert job in executor.jobs_with_errors
    assert out_dir.is_dir()
    assert (out_dir / "keep.txt").read_text() == "keep"


def test_local_executor_empty_dag_completes():
    """An empty DAG (e.g. everything skipped by resume) runs cleanly."""
    executor = LocalExecutor(ThreadScheduler(DAG(), 1))
    executor.execute()
    assert executor.jobs_with_errors == []


def test_dry_run_empty_dag_completes(capsys):
    """DryRunExecutor handles an empty DAG without printing anything."""
    executor = DryRunExecutor(ThreadScheduler(DAG(), 1))
    executor.execute()
    assert capsys.readouterr().out == ""
    assert executor.jobs_with_errors == []


def test_resume_after_failure_reruns_only_failed_subtree(tmp_path):
    """After a partial failure, skip_satisfied + has_all_outputs resumes
    the run from the surviving intermediate files."""
    a_out = tmp_path / "a.txt"
    b_out = tmp_path / "b.txt"

    def build_dag(b_cmd):
        dag = DAG()
        a = Job(
            Pipeline(Command("echo", "a"), file_output=a_out),
            "a",
            outputs=[a_out],
        )
        b = Job(Pipeline(b_cmd), "b", outputs=[b_out])
        dag.add_job(a)
        dag.add_job(b, {a})
        return dag, a, b

    # First run: a succeeds, b fails (and b's partial output is deleted).
    b_out.write_text("partial")
    dag, _a, b = build_dag(Command("false"))
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert b in executor.jobs_with_errors
    assert a_out.exists()
    assert not b_out.exists()

    # Second run: rebuild the DAG with b fixed; resume skips a.
    dag, a, b = build_dag(Command("touch", str(b_out)))
    skipped = dag.skip_satisfied(has_all_outputs)
    assert skipped == [a]
    executor = LocalExecutor(ThreadScheduler(dag, 1))
    executor.execute()
    assert executor.jobs_with_errors == []
    assert b_out.exists()
    assert a_out.read_text().strip() == "a"  # a was not rerun over


def _execute_bounded(executor, timeout=15):
    """Run executor.execute() with a hang guard."""
    holder = {}

    def run():
        try:
            executor.execute()
        except BaseException as exc:  # pragma: no cover - regression guard
            holder["exc"] = exc

    thread = threading.Thread(target=run, daemon=True)
    thread.start()
    thread.join(timeout=timeout)
    assert not thread.is_alive(), "executor.execute() hung"
    assert "exc" not in holder, holder.get("exc")


def test_procsub_job_with_early_exiting_outer_does_not_hang():
    """An outer command that exits without opening its proc-sub FIFO must
    not hang the executor (regression: cleanup deadlock)."""
    dag = DAG()
    job = Job(
        Pipeline(
            Command("false", InputProcSub(Pipeline(Command("echo", "x"))))
        ),
        "early-exit",
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)
    assert job in executor.jobs_with_errors


def test_procsub_job_with_unspawnable_outer_does_not_hang():
    """A missing outer executable with a proc-sub argument must not hang
    the executor's error path (regression: cleanup deadlock)."""
    dag = DAG()
    job = Job(
        Pipeline(
            Command(
                "no_such_cmd_zzz",
                InputProcSub(Pipeline(Command("echo", "x"))),
            )
        ),
        "ghost-with-procsub",
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)
    assert job in executor.jobs_with_errors


def test_procsub_inner_that_fails_to_launch_is_recorded():
    """An inner proc-sub command that cannot be spawned fails its job
    instead of crashing the executor (regression: a raw FileNotFoundError
    from the background task propagated out of execute())."""
    dag = DAG()
    job = Job(
        Pipeline(
            Command(
                "cat",
                InputProcSub(Pipeline(Command("no_such_cmd_zzz"))),
            )
        ),
        "inner-fail",
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)
    assert job in executor.jobs_with_errors


def test_procsub_inner_failure_does_not_skip_sibling_bookkeeping():
    """A job whose inner proc-sub fails to launch must not abort the run
    for independent sibling jobs (regression: the exception escaped mid-
    flight and skipped other jobs' bookkeeping)."""
    dag = DAG()
    bad = Job(
        Pipeline(
            Command(
                "cat",
                InputProcSub(Pipeline(Command("no_such_cmd_zzz"))),
            )
        ),
        "inner-fail",
    )
    good = Job(Pipeline(Command("true")), "sibling")
    dag.add_job(bad)
    dag.add_job(good)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)
    assert bad in executor.jobs_with_errors
    # The sibling ran to completion and was recorded as a success, not
    # silently dropped when the bad job's failure surfaced.
    assert good not in executor.jobs_with_errors


def test_mid_pipeline_spawn_failure_terminates_running_stages(monkeypatch):
    """When a later pipeline stage fails to spawn, earlier stages that are
    already running must be terminated, not leaked (regression: a
    Pipeline(sleep 600, no_such_cmd) left sleep alive past the workflow)."""
    import sentieon_cli.executor as executor_mod

    contexts = []
    real_context = executor_mod.Context

    class RecordingContext(real_context):  # type: ignore[valid-type,misc]
        def __init__(self) -> None:
            super().__init__()
            contexts.append(self)

    monkeypatch.setattr(executor_mod, "Context", RecordingContext)

    dag = DAG()
    job = Job(
        Pipeline(
            Command("sleep", "600"),
            Command("no_such_cmd_zzz"),
        ),
        "leak",
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)

    assert job in executor.jobs_with_errors
    # The first stage spawned; every spawned proc must be reaped, not left
    # running to outlive the workflow.
    spawned = [sub.proc for c in contexts for sub in c.commands if sub.proc]
    assert spawned, "expected the first stage to have spawned"
    for proc in spawned:
        assert proc.returncode is not None


def test_sigpipe_producer_does_not_fail_job(tmp_path):
    """A producer killed by SIGPIPE when a consumer exits early (e.g.
    `seq | head -1`) produces the right output and must not fail the job
    (regression: the producer's rc -13 was treated as a failure)."""
    out = tmp_path / "out.txt"
    dag = DAG()
    job = Job(
        Pipeline(
            Command("seq", "1000000"),
            Command("head", "-1"),
            file_output=out,
        ),
        "sigpipe",
    )
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)
    assert job not in executor.jobs_with_errors
    assert out.read_text().strip() == "1"


def test_nonzero_producer_still_fails_job():
    """SIGPIPE forgiveness must not mask a producer that exits non-zero
    for another reason: `false | cat` still fails the job."""
    dag = DAG()
    job = Job(Pipeline(Command("false"), Command("cat")), "real-fail")
    dag.add_job(job)
    executor = LocalExecutor(ThreadScheduler(dag, 2))
    _execute_bounded(executor)
    assert job in executor.jobs_with_errors


class _FakeContext:
    """A stand-in Context whose cleanup() can be told to raise."""

    def __init__(self, cleanup_error=None):
        self.commands = []
        self.cleaned = False
        self._error = cleanup_error

    async def cleanup(self):
        self.cleaned = True
        if self._error is not None:
            raise self._error


@pytest.mark.asyncio
async def test_cleanup_quietly_swallows_oserror():
    """A launch failure surfacing from cleanup() during teardown is logged,
    not raised -- so it cannot escape execute() on the abort/interrupt
    paths (symmetry with _drive's per-job handling)."""
    context = _FakeContext(FileNotFoundError("no_such_inner"))
    await _cleanup_quietly(context)  # must not raise
    assert context.cleaned


@pytest.mark.asyncio
async def test_cleanup_quietly_propagates_non_oserror():
    """A non-OSError from cleanup() is a bug and must still propagate."""
    with pytest.raises(RuntimeError):
        await _cleanup_quietly(_FakeContext(RuntimeError("bug")))


@pytest.mark.asyncio
async def test_shutdown_releases_all_contexts_despite_cleanup_error():
    """A cleanup() OSError during interrupt shutdown is swallowed and the
    remaining contexts are still released, instead of the error escaping
    _shutdown() and skipping their teardown."""
    executor = LocalExecutor(ThreadScheduler(DAG(), 1))
    bad = _FakeContext(FileNotFoundError("no_such_inner"))
    good = _FakeContext()

    async def _done():
        return 0

    executor.running = [
        (
            Job(Pipeline(Command("true")), "j1"),
            bad,
            asyncio.create_task(_done()),
            0,
        ),
        (
            Job(Pipeline(Command("false")), "j2"),
            good,
            asyncio.create_task(_done()),
            0,
        ),
    ]
    await executor._shutdown()  # must not raise

    assert bad.cleaned and good.cleaned
