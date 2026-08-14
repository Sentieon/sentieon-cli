"""
Signal-handling tests for LocalExecutor.

These deliver real signals to the current process. That is safe here because
the executor installs an asyncio loop signal handler that suppresses the
default disposition (e.g. SIGINT -> KeyboardInterrupt) for the duration of the
run, and the timer that raises the signal is always armed *after* the handler
is installed. If these ever prove flaky in a constrained CI, the same behavior
can be exercised by running the workflow in a subprocess and signalling the
child.
"""

import os
import signal
import sys
import threading
import time

# Add the parent directory to the path
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

import sentieon_cli.executor as executor_mod  # noqa: E402
from sentieon_cli.dag import DAG  # noqa: E402
from sentieon_cli.executor import LocalExecutor  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.scheduler import ThreadScheduler  # noqa: E402
from sentieon_cli.shell_pipeline import (  # noqa: E402
    Command,
    InputProcSub,
    Pipeline,
)


def _echo_dag(out):
    dag = DAG()
    job = Job(
        Pipeline(Command("echo", "hi"), file_output=out),
        "echo",
        task_name="test",
    )
    dag.add_job(job)
    return dag


def test_off_main_thread_does_not_crash(tmp_path):
    """LocalExecutor with signal handlers runs fine off the main thread."""
    out = tmp_path / "o.txt"
    holder = {}

    def run():
        try:
            LocalExecutor(
                ThreadScheduler(_echo_dag(out), 1),
                install_signal_handlers=True,
            ).execute()
        except Exception as exc:  # pragma: no cover - regression guard
            holder["exc"] = exc

    thread = threading.Thread(target=run)
    thread.start()
    thread.join(timeout=10)

    assert not thread.is_alive()
    assert "exc" not in holder, holder.get("exc")
    assert out.read_text().strip() == "hi"


def test_default_does_not_touch_signal_handlers(tmp_path):
    """The default executor installs no handlers and leaves them intact."""
    before_int = signal.getsignal(signal.SIGINT)
    before_term = signal.getsignal(signal.SIGTERM)

    out = tmp_path / "o.txt"
    LocalExecutor(ThreadScheduler(_echo_dag(out), 1)).execute()

    assert signal.getsignal(signal.SIGINT) is before_int
    assert signal.getsignal(signal.SIGTERM) is before_term


def test_opt_in_restores_signal_handlers(tmp_path):
    """Opt-in handlers are restored to their prior values after a run."""
    before_int = signal.getsignal(signal.SIGINT)
    before_term = signal.getsignal(signal.SIGTERM)

    out = tmp_path / "o.txt"
    LocalExecutor(
        ThreadScheduler(_echo_dag(out), 1),
        install_signal_handlers=True,
    ).execute()

    assert signal.getsignal(signal.SIGINT) is before_int
    assert signal.getsignal(signal.SIGTERM) is before_term


def test_sigint_triggers_graceful_shutdown(monkeypatch):
    """A delivered SIGINT terminates running jobs and cleans up."""
    temp_dirs = []
    real_context = executor_mod.Context

    class SpyContext(real_context):
        async def cleanup(self):
            temp_dirs.append(self.temp_dir.name)
            await super().cleanup()

    monkeypatch.setattr(executor_mod, "Context", SpyContext)

    dag = DAG()
    dag.add_job(
        Job(Pipeline(Command("sleep", "30")), "sleeper", task_name="test")
    )
    executor = LocalExecutor(
        ThreadScheduler(dag, 1),
        install_signal_handlers=True,
        shutdown_grace_period=1.0,
    )

    timer = threading.Timer(0.5, lambda: os.kill(os.getpid(), signal.SIGINT))
    timer.start()
    try:
        start = time.monotonic()
        executor.execute()
        elapsed = time.monotonic() - start
    finally:
        timer.cancel()

    assert elapsed < 10  # returned well before `sleep 30` would finish
    assert executor.running == []
    assert temp_dirs  # cleanup ran for the interrupted job
    for temp_dir in temp_dirs:
        assert not os.path.exists(temp_dir)  # no leaked temp dir


def test_grace_period_escalates_to_sigkill():
    """A job that ignores SIGTERM is still killed after the grace period."""
    script = (
        "import signal, time; "
        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
        "time.sleep(30)"
    )
    dag = DAG()
    job = Job(
        Pipeline(Command(sys.executable, "-c", script)),
        "stubborn",
        task_name="test",
    )
    dag.add_job(job)
    executor = LocalExecutor(
        ThreadScheduler(dag, 1),
        install_signal_handlers=True,
        shutdown_grace_period=0.5,
    )

    timer = threading.Timer(0.5, lambda: os.kill(os.getpid(), signal.SIGINT))
    timer.start()
    try:
        start = time.monotonic()
        executor.execute()
        elapsed = time.monotonic() - start
    finally:
        timer.cancel()

    assert elapsed < 10  # SIGKILL fired after grace, not after time.sleep(30)
    assert executor.running == []


def test_sigint_with_blocked_procsub_cleans_up():
    """SIGINT while a proc-sub FIFO open is still blocked must not hang
    the shutdown path (regression: cleanup deadlock)."""
    dag = DAG()
    # sh never opens the FIFO in "$1", so the proc sub's open() stays
    # blocked until cleanup unblocks it.
    job = Job(
        Pipeline(
            Command(
                "sh",
                "-c",
                "sleep 30",
                "sh",
                InputProcSub(Pipeline(Command("echo", "x"))),
            )
        ),
        "blocked-procsub",
        task_name="test",
    )
    dag.add_job(job)
    executor = LocalExecutor(
        ThreadScheduler(dag, 2),
        install_signal_handlers=True,
        shutdown_grace_period=1.0,
    )

    timer = threading.Timer(0.5, lambda: os.kill(os.getpid(), signal.SIGINT))
    timer.start()
    try:
        start = time.monotonic()
        executor.execute()
        elapsed = time.monotonic() - start
    finally:
        timer.cancel()

    assert elapsed < 10
    assert executor.running == []
