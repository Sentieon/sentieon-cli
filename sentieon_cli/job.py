"""
Job objects
"""

from typing import Dict, Optional

from .shell_pipeline import Pipeline


class Job:
    """A unit of work: a shell pipeline plus its execution metadata.

    Fields:

    * ``shell`` -- the pipeline to run.
    * ``name`` -- a human-readable label (not part of identity).
    * ``threads`` -- CPU threads the job needs (a local scheduling budget).
    * ``resources`` -- named resource counts (e.g. NUMA-node tokens).

    A job's identity is keyed only on ``shell``: two jobs with the same
    pipeline are equal (and collide in a DAG) regardless of the other fields.
    """

    def __init__(
        self,
        pipeline: Pipeline,
        name: str,
        threads: int = 1,
        resources: Optional[Dict[str, int]] = None,
    ) -> None:
        self.shell = pipeline
        self.name = name
        self.threads = threads
        self.resources = {} if resources is None else resources

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
