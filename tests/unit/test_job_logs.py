"""
Unit tests for the per-process job log files
"""

import os
import sys

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.executor import TAIL_LINES, _log_tail  # noqa: E402
from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.run_logs import RunLogs  # noqa: E402
from sentieon_cli.shell_pipeline import (  # noqa: E402
    Command,
    Context,
    Pipeline,
)
from sentieon_cli.util import sanitize  # noqa: E402


def _sink(tmp_path, name: str = "bwa", task_name: str = "alignment"):
    job = Job(Pipeline(Command("true")), name, task_name=task_name)
    return RunLogs(tmp_path / "logs").job_sink(job)


def test_sanitize_keeps_safe_characters():
    assert sanitize("bwa-1.0_x") == "bwa-1.0_x"
    assert sanitize("call vars/chr1:1-2") == "call-vars-chr1-1-2"


def test_stage_indices_follow_spawn_order_across_both_streams(tmp_path):
    sink = _sink(tmp_path)
    for handle in (
        sink.open_stderr(Command("a")),
        sink.open_stdout(Command("b")),
        sink.open_stderr(Command("c")),
    ):
        handle.close()

    assert [path.name for path in sink.log_paths()] == [
        "bwa-1.0.log",
        "bwa-1.1.stdout.log",
        "bwa-1.2.log",
    ]


def test_both_streams_of_one_process_share_a_stage_index(tmp_path):
    sink = _sink(tmp_path)
    first = Command("a")
    second = Command("b")
    for handle in (
        sink.open_stderr(first),
        sink.open_stdout(first),
        sink.open_stderr(second),
        sink.open_stdout(second),
    ):
        handle.close()

    assert [path.name for path in sink.log_paths()] == [
        "bwa-1.0.log",
        "bwa-1.0.stdout.log",
        "bwa-1.1.log",
        "bwa-1.1.stdout.log",
    ]


def test_the_task_directory_is_created_lazily(tmp_path):
    sink = _sink(tmp_path)
    assert not sink.log_dir.exists()

    sink.open_stderr(Command("true")).close()

    assert sink.log_dir == tmp_path / "logs" / "task_logs" / "alignment"
    assert sink.log_dir.is_dir()


def test_path_components_are_sanitized(tmp_path):
    sink = _sink(tmp_path, name="call vars/2", task_name="variant calling")
    sink.open_stderr(Command("true")).close()

    path = sink.log_paths()[0]
    assert path.parent.name == "variant-calling"
    assert path.name == "call-vars-2-1.0.log"


def test_the_stderr_log_starts_with_a_header(tmp_path):
    sink = _sink(tmp_path)
    handle = sink.open_stderr(Command("echo", "hello world"))
    handle.write(b"from the child\n")
    handle.close()

    contents = sink.log_paths()[0].read_text()
    assert "# timestamp: " in contents
    assert "# task: alignment" in contents
    assert "# job: bwa (bwa-1)" in contents
    assert "# stage: 0" in contents
    assert "# command: echo 'hello world'" in contents
    # The header is flushed before the child inherits the FD, so child
    # output always lands after it.
    assert contents.splitlines()[-1] == "from the child"


def test_the_stdout_log_has_no_header(tmp_path):
    sink = _sink(tmp_path)
    handle = sink.open_stdout(Command("echo", "hello"))
    handle.write(b"payload\n")
    handle.close()

    assert sink.log_paths()[0].read_bytes() == b"payload\n"


def test_logs_are_matched_to_commands_by_identity(tmp_path):
    sink = _sink(tmp_path)
    first = Command("echo", "same")
    second = Command("echo", "same")
    assert first == second  # equal, but two distinct processes

    sink.open_stderr(first).close()
    sink.open_stderr(second).close()

    assert sink.stderr_log_for(first).name == "bwa-1.0.log"
    assert sink.stderr_log_for(second).name == "bwa-1.1.log"
    assert sink.stderr_log_for(Command("elsewhere")) is None


def test_stdout_logs_are_not_returned_as_stderr_logs(tmp_path):
    sink = _sink(tmp_path)
    command = Command("echo", "hi")
    sink.open_stdout(command).close()

    assert sink.stderr_log_for(command) is None


def test_finalize_prunes_uninformative_logs_after_a_success(tmp_path):
    sink = _sink(tmp_path)
    header_only = Command("quiet")
    empty_stdout = Command("quiet")
    noisy = Command("loud")
    sink.open_stderr(header_only).close()
    sink.open_stdout(empty_stdout).close()
    handle = sink.open_stderr(noisy)
    handle.write(b"something happened\n")
    handle.close()

    sink.finalize(success=True)

    assert not sink.stderr_log_for(header_only).exists()
    assert not sink.log_paths()[1].exists()
    assert sink.stderr_log_for(noisy).exists()


def test_finalize_keeps_every_log_after_a_failure(tmp_path):
    sink = _sink(tmp_path)
    sink.open_stderr(Command("quiet")).close()
    sink.open_stdout(Command("quiet")).close()

    sink.finalize(success=False)

    assert all(path.exists() for path in sink.log_paths())


def test_a_context_without_a_sink_captures_nothing():
    assert Context().log_sink is None


def test_log_tail_returns_the_last_lines(tmp_path):
    log = tmp_path / "big.log"
    log.write_text("".join(f"line {i}\n" for i in range(100)))

    tail = _log_tail(log)

    assert len(tail) == TAIL_LINES
    assert tail[0] == "line 80"
    assert tail[-1] == "line 99"


def test_log_tail_returns_a_short_file_whole(tmp_path):
    log = tmp_path / "small.log"
    log.write_text("one\ntwo\nthree\n")

    assert _log_tail(log) == ["one", "two", "three"]


def test_log_tail_only_reads_the_end_of_a_large_file(tmp_path):
    log = tmp_path / "huge.log"
    padding = "x" * 200
    with open(log, "w") as handle:
        for i in range(2000):
            handle.write(f"line {i} {padding}\n")
    assert log.stat().st_size > 64 * 1024

    tail = _log_tail(log)

    assert len(tail) == TAIL_LINES
    assert tail[-1].startswith("line 1999 ")
    assert not any(line.startswith("line 0 ") for line in tail)


def test_log_tail_of_a_missing_file_is_empty(tmp_path):
    assert _log_tail(tmp_path / "gone.log") == []
