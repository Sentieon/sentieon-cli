"""
Repeat expansion genotyping with ExpansionHunter
"""

from dataclasses import dataclass
import pathlib
from typing import Dict, Iterable, Optional

import packaging.version

from .. import command_strings as cmds
from ..dag import DAG
from ..job import Job
from ..util import SampleSex, caller_sex_arg
from .base import Stage, StageResult

EXPANSION_MIN_VERSIONS: Dict[str, Optional[packaging.version.Version]] = {
    "ExpansionHunter": None,
}


@dataclass(kw_only=True)
class ExpansionHunterResult(StageResult):
    """The job and output of `ExpansionHunterStage`"""

    job: Job
    output_prefix: pathlib.Path


@dataclass(kw_only=True)
class ExpansionHunterStage(Stage):
    """Genotype repeat expansions with ExpansionHunter.

    ExpansionHunter is sex-aware, so callers run this stage after ploidy
    estimation. ``output_prefix`` is the `--output-prefix` the tool
    writes its files under.
    """

    alignment: pathlib.Path
    variant_catalog: pathlib.Path
    output_prefix: pathlib.Path
    sample_sex: Optional[SampleSex] = None
    name: str = "expansion-hunter"
    task_name: str = "expansion-hunter"

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> ExpansionHunterResult:
        deps = set(upstream)

        job = Job(
            cmds.cmd_expansion_hunter(
                self.output_prefix,
                self.alignment,
                reference=self.ctx.reference,
                variant_catalog=self.variant_catalog,
                sex=caller_sex_arg(self.sample_sex),
                threads=self.ctx.cores,
            ),
            self.name,
            self.ctx.cores,
            task_name=self.task_name,
        )
        dag.add_job(job, deps)

        return ExpansionHunterResult(
            jobs=[job],
            terminal={job},
            job=job,
            output_prefix=self.output_prefix,
        )
