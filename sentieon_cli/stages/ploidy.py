"""
Sample ploidy and sex estimation
"""

from dataclasses import dataclass
from importlib.resources import files
import pathlib
from typing import Iterable, List, Optional

from .. import command_strings as cmds
from ..dag import DAG
from ..job import Job
from ..shard import ploidy_contigs_for_build
from .base import Stage, StageResult


@dataclass(kw_only=True)
class PloidyResult(StageResult):
    """The job and output of `PloidyStage`"""

    job: Job
    ploidy_json: pathlib.Path


@dataclass(kw_only=True)
class PloidyStage(Stage):
    """Estimate the sample ploidy and sex from an alignment.

    The JSON report is written next to the output VCF and is read after
    the DAG has run, by the pipelines that call variants with a
    sex-aware caller.
    """

    inputs: List[pathlib.Path]
    reference_build: Optional[str] = None

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> PloidyResult:
        deps = set(upstream)
        ploidy_json = pathlib.Path(
            str(self.ctx.output_vcf).replace(".vcf.gz", "_ploidy.json")
        )
        estimate_ploidy = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("estimate_ploidy.py"))
        ).resolve()
        ploidy_contigs = ploidy_contigs_for_build(self.reference_build)
        job = Job(
            cmds.cmd_estimate_ploidy(
                ploidy_json,
                self.inputs,
                estimate_ploidy,
                contigs=ploidy_contigs.contigs,
                autosomes=ploidy_contigs.autosomes,
                x_contig=ploidy_contigs.x_contig,
                y_contig=ploidy_contigs.y_contig,
            ),
            "estimate-ploidy",
            0,
            task_name="ploidy",
        )
        dag.add_job(job, deps)

        return PloidyResult(
            jobs=[job],
            terminal={job},
            job=job,
            ploidy_json=ploidy_json,
        )
