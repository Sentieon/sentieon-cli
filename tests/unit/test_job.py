"""
Unit tests for the Job model.
"""

import os
import pathlib
import sys

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

from sentieon_cli.job import Job  # noqa: E402
from sentieon_cli.shell_pipeline import Command, Pipeline  # noqa: E402


def _pipe(text):
    return Pipeline(Command("echo", text))


def test_new_fields_default_empty():
    job = Job(_pipe("a"), "a")
    assert job.inputs == []
    assert job.outputs == []
    assert job.image is None


def test_new_fields_are_stored():
    ins = [pathlib.Path("in1"), pathlib.Path("in2")]
    outs = [pathlib.Path("out")]
    job = Job(_pipe("a"), "a", inputs=ins, outputs=outs, image="org/img:1")
    assert job.inputs == ins
    assert job.outputs == outs
    assert job.image == "org/img:1"


def test_input_output_sequences_are_snapshotted():
    ins = [pathlib.Path("in1")]
    job = Job(_pipe("a"), "a", inputs=ins)
    ins.append(pathlib.Path("in2"))  # mutating the caller's list...
    assert job.inputs == [pathlib.Path("in1")]  # ...doesn't change the job


def test_identity_ignores_metadata_fields():
    # Identity is keyed on the shell pipeline; the other fields don't
    # affect it.
    j1 = Job(_pipe("same"), "j1", inputs=[pathlib.Path("x")], image="a")
    j2 = Job(_pipe("same"), "j2", outputs=[pathlib.Path("y")], image="b")
    assert j1 == j2
    assert hash(j1) == hash(j2)
