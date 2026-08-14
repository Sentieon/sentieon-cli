"""
Unit tests for shell_pipeline.py
"""

import asyncio
import fcntl
import os
import pathlib
import signal
import sys
import tempfile

import pytest

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.shell_pipeline import (  # noqa: E402
    Command,
    Context,
    InputProcSub,
    OutputProcSub,
    Pipeline,
)


@pytest.mark.asyncio
async def test_simple_command():
    """Test running a simple command"""
    cmd = Command("echo", "hello world")
    context = Context()

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await cmd.run(context, stdout=stdout_file)
        await proc.async_wait()

        stdout_file.seek(0)
        output = stdout_file.read().strip()

    await context.cleanup()
    os.remove(stdout_file.name)

    assert output == "hello world"
    assert proc.returncode == 0
    # wait4 captured the child's resource usage as it was reaped
    assert proc.rusage is not None
    assert proc.rusage.ru_maxrss > 0


@pytest.mark.asyncio
async def test_a_signalled_command_reports_its_signal_and_its_rusage():
    """A child killed by a signal is still reaped through ``wait4``.

    This pins ``RusagePopen._try_wait``: if a future CPython changes that
    private method's shape, the override stops being called and ``rusage``
    silently stays None while the negative return code still works.
    """
    cmd = Command("sleep", "5")
    context = Context()

    proc = await cmd.run(context)
    proc.send_signal(signal.SIGTERM)

    assert await proc.async_wait() == -signal.SIGTERM
    assert proc.rusage is not None

    await context.cleanup()


@pytest.mark.asyncio
async def test_simple_pipeline():
    """Test a simple pipeline: echo hello | cat"""
    cmd1 = Command("echo", "hello pipeline")
    cmd2 = Command("cat")
    pipeline = Pipeline(cmd1, cmd2)
    context = Context()

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await pipeline.run(context, stdout=stdout_file)
        await proc.async_wait()

        stdout_file.seek(0)
        output = stdout_file.read().strip()

    await context.cleanup()
    os.remove(stdout_file.name)

    assert output == "hello pipeline"
    assert proc.returncode == 0


@pytest.mark.asyncio
async def test_pipeline_with_file_io():
    """Test a pipeline with file input and output"""
    context = Context()

    with tempfile.NamedTemporaryFile(mode="w", delete=False) as infile:
        infile.write("hello from file")
        infile_path = pathlib.Path(infile.name)

    outfile_path = pathlib.Path(context.temp_dir.name) / "output.txt"

    cmd1 = Command("cat")
    cmd2 = Command("wc", "-c")
    pipeline = Pipeline(
        cmd1,
        cmd2,
        file_input=infile_path,
        file_output=outfile_path,
    )

    proc = await pipeline.run(context)
    await proc.async_wait()

    with open(outfile_path, "r") as f:
        output = f.read().strip()

    await context.cleanup()
    os.remove(infile_path)

    assert output == str(len("hello from file"))
    assert proc.returncode == 0


@pytest.mark.asyncio
async def test_input_process_substitution():
    """Test input process substitution: diff <(echo 1) <(echo 2)"""
    context = Context()

    cmd = Command(
        "diff",
        InputProcSub(Pipeline(Command("echo", "1"))),
        InputProcSub(Pipeline(Command("echo", "2"))),
    )

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await cmd.run(context, stdout=stdout_file)
        await proc.async_wait()

        stdout_file.seek(0)
        output = stdout_file.read().strip()

    await context.cleanup()
    os.remove(stdout_file.name)

    # The diff output should indicate a change from 1 to 2
    assert "1c1" in output
    assert "< 1" in output
    assert "---" in output
    assert "> 2" in output
    assert proc.returncode != 0  # diff exits with 1 if inputs differ


@pytest.mark.asyncio
async def test_output_process_substitution():
    """Test output process substitution: tee >(wc -c)"""
    context = Context()

    outfile_path = pathlib.Path(context.temp_dir.name) / "wc_output.txt"

    # echo "hello" | tee >(wc -c > outfile.txt)
    # The main output of the pipeline should be "hello"
    # The output of the process substitution should be the char count.
    wc_cmd = Command("wc", "-c")
    wc_pipeline = Pipeline(wc_cmd, file_output=outfile_path)

    echo_cmd = Command("echo", "hello")
    tee_cmd = Command("tee", OutputProcSub(wc_pipeline))
    main_pipeline = Pipeline(echo_cmd, tee_cmd)

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await main_pipeline.run(context, stdout=stdout_file)
        await proc.async_wait()

        # Wait for background tasks from process substitution to finish
        await asyncio.gather(*context.tasks)

        stdout_file.seek(0)
        main_output = stdout_file.read().strip()

    assert main_output == "hello"

    with open(outfile_path, "r") as f:
        wc_output = f.read().strip()

    # 'echo' adds a newline, so it should be 6 characters for `wc -c`
    assert wc_output == str(len("hello\n"))

    await context.cleanup()
    os.remove(stdout_file.name)
    assert proc.returncode == 0


@pytest.mark.asyncio
async def test_pipe_size_enlarges_internal_pipes(tmp_path):
    """pipe_size raises F_SETPIPE_SZ on the pipeline's internal pipes,
    which the downstream child inherits."""
    if not hasattr(fcntl, "F_SETPIPE_SZ"):
        pytest.skip("F_SETPIPE_SZ not available (non-Linux)")
    try:
        pipe_max = int(open("/proc/sys/fs/pipe-max-size").read())
    except OSError:
        pytest.skip("cannot read fs.pipe-max-size")
    target = 128 * 1024  # a power of two, so the kernel sets it exactly
    if pipe_max < target:
        pytest.skip("fs.pipe-max-size too small for this test")

    report = tmp_path / "size.txt"
    # The first stage's stdout is the internal pipe to `cat`; report that
    # pipe's buffer size (1032 == F_GETPIPE_SZ) rather than writing to it.
    reporter = (
        "import fcntl, sys; "
        "open(sys.argv[1], 'w').write(str(fcntl.fcntl(1, 1032)))"
    )
    context = Context()
    pipeline = Pipeline(
        Command(sys.executable, "-c", reporter, str(report)),
        Command("cat"),
        pipe_size=target,
    )
    proc = await pipeline.run(context)
    await proc.async_wait()
    for sub in context.commands:
        if sub.proc:
            await sub.proc.async_wait()
    await context.cleanup()

    assert int(report.read_text()) == target


def test_command_and_pipeline_hash_identity():
    """Equal commands/pipelines hash equally (DAG relies on this)."""
    c1 = Command("echo", "x")
    c2 = Command("echo", "x")
    c3 = Command("echo", "y")
    assert c1 == c2
    assert hash(c1) == hash(c2)
    assert c1 != c3

    p1 = Pipeline(Command("echo", "x"), Command("cat"))
    p2 = Pipeline(Command("echo", "x"), Command("cat"))
    p3 = Pipeline(Command("echo", "x"))
    assert p1 == p2
    assert hash(p1) == hash(p2)
    assert p1 != p3


def test_procsub_input_and_output_are_distinct():
    """<(x) and >(x) wrap the same pipeline but are not the same node: they
    must be unequal and hash differently."""
    p = Pipeline(Command("echo", "x"))
    assert InputProcSub(p) != OutputProcSub(p)
    assert hash(InputProcSub(p)) != hash(OutputProcSub(p))
    # Each still equals another of its own kind wrapping an equal pipeline.
    other = Pipeline(Command("echo", "x"))
    assert InputProcSub(p) == InputProcSub(other)
    assert hash(InputProcSub(p)) == hash(InputProcSub(other))


def test_command_identity_is_stable_across_running():
    """A Command's identity must not change once it has started: __eq__
    used to mix in self.proc, so running a job mutated its identity."""
    c1 = Command("echo", "x")
    c2 = Command("echo", "x")
    assert c1 == c2 and hash(c1) == hash(c2)
    c1.proc = object()  # stand in for a started asyncio process
    assert c1 == c2
    assert hash(c1) == hash(c2)


def test_command_identity_includes_exec_kwargs():
    """exec_kwargs participates in identity: commands differing only in
    env/cwd are distinct, and equal env hashes the same regardless of the
    dict's insertion order."""
    assert Command("run", exec_kwargs={"cwd": "/a"}) != Command(
        "run", exec_kwargs={"cwd": "/b"}
    )
    assert Command("run", exec_kwargs={"cwd": "/a"}) != Command("run")
    e1 = Command("run", exec_kwargs={"env": {"A": "1", "B": "2"}})
    e2 = Command("run", exec_kwargs={"env": {"B": "2", "A": "1"}})
    assert e1 == e2
    assert hash(e1) == hash(e2)


def test_str_renders_env_and_cwd():
    """str() renders exec_kwargs env/cwd as a conventional bash prefix,
    quoting values that need it."""
    env_cmd = Command("samtools", "sort", exec_kwargs={"env": {"REF": "hg"}})
    assert str(env_cmd) == "REF=hg samtools sort"
    cwd_cmd = Command("samtools", "sort", exec_kwargs={"cwd": "/data"})
    assert str(cwd_cmd) == "(cd /data && samtools sort)"
    both = Command(
        "run", exec_kwargs={"cwd": "/data dir", "env": {"K": "v w"}}
    )
    assert str(both) == "(cd '/data dir' && K='v w' run)"


def test_empty_pipeline_raises():
    with pytest.raises(ValueError):
        Pipeline()


@pytest.mark.asyncio
async def test_file_output_and_stdout_conflict_raises():
    """A file/stdout conflict raises ValueError (not sys.exit)."""
    context = Context()
    out = pathlib.Path(context.temp_dir.name) / "out.txt"
    pipeline = Pipeline(Command("echo", "x"), file_output=out)
    with pytest.raises(ValueError):
        await pipeline.run(context, stdout=1)
    await context.cleanup()


@pytest.mark.asyncio
async def test_file_input_and_stdin_fd0_conflict_raises(tmp_path):
    """file_input plus a stdin of fd 0 conflicts: fd 0 is a valid stream but
    falsy, so the guard must test `is not None`, not truthiness."""
    inp = tmp_path / "in.txt"
    inp.write_text("data")
    context = Context()
    pipeline = Pipeline(Command("cat"), file_input=inp)
    with pytest.raises(ValueError):
        await pipeline.run(context, stdin=0)
    await context.cleanup()


@pytest.mark.asyncio
async def test_file_input_and_stdin_conflict_raises():
    """A file/stdin conflict raises ValueError (not sys.exit)."""
    context = Context()
    infile = pathlib.Path(context.temp_dir.name) / "in.txt"
    pipeline = Pipeline(Command("cat"), file_input=infile)
    with pytest.raises(ValueError):
        await pipeline.run(context, stdin=1)
    await context.cleanup()


@pytest.mark.asyncio
async def test_cleanup_unblocks_unopened_input_procsub():
    """cleanup() must not hang when an <(...) FIFO was never opened."""
    context = Context()
    proc_sub = InputProcSub(Pipeline(Command("echo", "x")))
    fifo_path = await proc_sub.setup(context)
    assert os.path.exists(fifo_path)

    await asyncio.wait_for(context.cleanup(), timeout=10)

    # The inner command never spawned: nothing opened the FIFO.
    assert context.commands == []


@pytest.mark.asyncio
async def test_cleanup_unblocks_unopened_output_procsub():
    """cleanup() must not hang when a >(...) FIFO was never opened."""
    context = Context()
    proc_sub = OutputProcSub(Pipeline(Command("cat")))
    await proc_sub.setup(context)

    await asyncio.wait_for(context.cleanup(), timeout=10)

    assert context.commands == []


@pytest.mark.asyncio
async def test_cleanup_after_outer_exits_without_opening():
    """An outer command that exits without opening its proc-sub FIFO
    (here: false ignores its arguments) must not hang cleanup()."""
    context = Context()
    cmd = Command("false", InputProcSub(Pipeline(Command("echo", "x"))))
    proc = await cmd.run(context)
    await proc.async_wait()

    await asyncio.wait_for(context.cleanup(), timeout=10)

    assert [c.executable for c in context.commands] == ["false"]


@pytest.mark.asyncio
async def test_cleanup_unblocks_multiple_procsubs():
    """Several never-opened proc subs are all unblocked by cleanup()."""
    context = Context()
    cmd = Command(
        "false",
        InputProcSub(Pipeline(Command("echo", "x"))),
        OutputProcSub(Pipeline(Command("cat"))),
    )
    proc = await cmd.run(context)
    await proc.async_wait()

    await asyncio.wait_for(context.cleanup(), timeout=10)

    assert all(task.done() for task in context.tasks)


@pytest.mark.asyncio
async def test_cleanup_still_raises_inner_launch_failure():
    """An inner command that cannot be spawned still propagates from
    cleanup(), but only after all resources have been released."""
    context = Context()
    cmd = Command(
        "cat", InputProcSub(Pipeline(Command("no_such_cmd_zzz")))
    )
    proc = await cmd.run(context)
    await proc.async_wait()

    with pytest.raises(FileNotFoundError):
        await asyncio.wait_for(context.cleanup(), timeout=10)

    assert not os.path.exists(context.temp_dir.name)


@pytest.mark.asyncio
async def test_cleanup_does_not_stall_running_inner_writer():
    """cleanup() must not touch the FIFO of a proc sub whose inner command
    is already running; holding a read end there would prevent the EPIPE
    that lets the inner writer exit."""
    context = Context()
    cmd = Command(
        "sh",
        "-c",
        'sleep 0.5 <"$1"',
        "sh",
        InputProcSub(Pipeline(Command("yes"))),
    )
    proc = await cmd.run(context)
    # Wait until the inner command has spawned (sh + yes).
    async with asyncio.timeout(5):
        while len(context.commands) < 2:
            await asyncio.sleep(0.01)

    await asyncio.wait_for(context.cleanup(), timeout=10)
    assert await proc.async_wait() == 0


@pytest.mark.asyncio
async def test_cleanup_twice_is_safe():
    """A second cleanup() call is a harmless no-op."""
    context = Context()
    cmd = Command("cat", InputProcSub(Pipeline(Command("echo", "hello"))))
    proc = await cmd.run(context, stdout=asyncio.subprocess.DEVNULL)
    await proc.async_wait()

    await asyncio.wait_for(context.cleanup(), timeout=10)
    await asyncio.wait_for(context.cleanup(), timeout=10)
    assert proc.returncode == 0
