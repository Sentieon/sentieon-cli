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
    * ``task_name`` -- the pipeline stage this job belongs to; it groups the
      job's log files, so every shard of an operation shares one value.
    * ``job_id`` -- ``{name}-{n}``, unique across a run (not part of
      identity).

    A job's identity is keyed only on ``shell``: two jobs with the same
    pipeline are equal (and collide in a DAG) regardless of the other fields.
    So two Job objects built from identical pipelines are "the same job" to
    the DAG even though each carries its own ``job_id``.
    """

    # Per-name counters backing ``job_id``; see ``reset_ids``.
    _id_counters: Dict[str, int] = {}

    def __init__(
        self,
        pipeline: Pipeline,
        name: str,
        threads: int = 1,
        resources: Optional[Dict[str, int]] = None,
        *,
        task_name: str,
    ) -> None:
        self.shell = pipeline
        self.name = name
        self.threads = threads
        self.resources = {} if resources is None else resources
        self.task_name = task_name
        count = Job._id_counters.get(name, 0) + 1
        Job._id_counters[name] = count
        self.job_id = f"{name}-{count}"

    @classmethod
    def reset_ids(cls) -> None:
        """Restart id numbering.

        Ids must stay unique for a whole run, which can execute more than one
        DAG, so this is called once per CLI invocation -- never between DAGs.
        """
        cls._id_counters.clear()

    def __hash__(self) -> int:
        return hash(self.shell)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Job):
            return self.shell == other.shell
        return False

    def __ne__(self, other: object) -> bool:
        return not self == other

    def __repr__(self) -> str:
        return f"Job({self.job_id})"

    def __str__(self) -> str:
        return f"Job({self.job_id})"
