"""
Unit tests for shell_pipeline.py
"""

import asyncio
import os
import pathlib
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

    # Use a file for stdout
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await cmd.run(context, stdout=stdout_file)
        await proc.wait()

        stdout_file.seek(0)
        output = stdout_file.read().strip()

    await context.cleanup()
    os.remove(stdout_file.name)

    assert output == "hello world"
    assert proc.returncode == 0


@pytest.mark.asyncio
async def test_simple_pipeline():
    """Test a simple pipeline: echo hello | cat"""
    cmd1 = Command("echo", "hello pipeline")
    cmd2 = Command("cat")
    pipeline = Pipeline(cmd1, cmd2)
    context = Context()

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await pipeline.run(context, stdout=stdout_file)
        await proc.wait()

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

    # Create input and output files
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as infile:
        infile.write("hello from file")
        infile_path = pathlib.Path(infile.name)

    outfile_path = pathlib.Path(context.temp_dir.name) / "output.txt"

    cmd1 = Command("cat")
    cmd2 = Command("wc", "-c")
    pipeline = Pipeline(cmd1, cmd2, file_input=infile_path, file_output=outfile_path)

    proc = await pipeline.run(context)
    await proc.wait()

    # Verify output file content
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
        await proc.wait()

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

    # We will write the output of wc -c to a file to check it
    outfile_path = pathlib.Path(context.temp_dir.name) / "wc_output.txt"

    # The pipeline will be: echo "hello" | tee >(wc -c > outfile.txt)
    # The main output of the pipeline should be "hello"
    # The output of the process substitution should be the character count "6"

    # Inner pipeline for process substitution
    wc_cmd = Command("wc", "-c")
    wc_pipeline = Pipeline(wc_cmd, file_output=outfile_path)

    # Main pipeline
    echo_cmd = Command("echo", "hello")
    tee_cmd = Command("tee", OutputProcSub(wc_pipeline))
    main_pipeline = Pipeline(echo_cmd, tee_cmd)

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as stdout_file:
        proc = await main_pipeline.run(context, stdout=stdout_file)
        await proc.wait()

        # Wait for background tasks from process substitution to finish
        await asyncio.gather(*context.tasks)

        stdout_file.seek(0)
        main_output = stdout_file.read().strip()

    # Verify main output
    assert main_output == "hello"

    # Verify process substitution output
    with open(outfile_path, "r") as f:
        wc_output = f.read().strip()

    # 'echo' adds a newline, so it should be 6 characters for `wc -c`
    assert wc_output == str(len("hello\n"))

    await context.cleanup()
    os.remove(stdout_file.name)
    assert proc.returncode == 0
