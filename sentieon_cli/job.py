"""
Job objects
"""

import pathlib
from typing import Dict, Iterable, List, Optional

from .shell_pipeline import Pipeline


class Job:
    """A unit of work: a shell pipeline plus its execution metadata.

    Fields:

    * ``shell`` -- the pipeline to run.
    * ``name`` -- a human-readable label (not part of identity).
    * ``threads`` -- CPU threads the job needs (a local scheduling budget).
    * ``resources`` -- named resource counts (e.g. NUMA-node tokens).
    * ``inputs`` / ``outputs`` -- local files the job reads / writes; when
      declared, the executor verifies them around the run and removes the
      outputs of a failed attempt.
    * ``image`` -- container image (or environment) to run in on a backend.

    A job's identity is keyed only on ``shell``: two jobs with the same
    pipeline are equal (and collide in a DAG) regardless of the other fields.
    """

    def __init__(
        self,
        pipeline: Pipeline,
        name: str,
        threads: int = 1,
        resources: Optional[Dict[str, int]] = None,
        *,
        inputs: Optional[Iterable[pathlib.Path]] = None,
        outputs: Optional[Iterable[pathlib.Path]] = None,
        image: Optional[str] = None,
    ) -> None:
        self.shell = pipeline
        self.name = name
        self.threads = threads
        self.resources = {} if resources is None else resources
        self.inputs: List[pathlib.Path] = (
            [] if inputs is None else list(inputs)
        )
        self.outputs: List[pathlib.Path] = (
            [] if outputs is None else list(outputs)
        )
        self.image = image

    def __hash__(self) -> int:
        return hash(self.shell)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Job):
            return self.shell == other.shell
        return False

    def __ne__(self, other: object) -> bool:
        return not self == other

    def __repr__(self) -> str:
        return f"Job({self.name})"

    def __str__(self) -> str:
        return f"Job({self.name})"
