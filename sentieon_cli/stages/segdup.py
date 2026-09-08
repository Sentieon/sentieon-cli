"""
Variant calling in segmental duplications with segdup-caller
"""

from dataclasses import dataclass
import pathlib
from typing import Dict, Iterable, Optional, Sequence

import packaging.version

from .. import command_strings as cmds
from ..dag import DAG
from ..job import Job
from ..util import SampleSex, caller_sex_arg
from .base import Stage, StageResult

SEGDUP_MIN_VERSIONS: Dict[str, Optional[packaging.version.Version]] = {
    # 0.7.0 adds the `--set` argument
    "segdup-caller": packaging.version.Version("0.7.0"),
}


@dataclass(kw_only=True)
class SegdupResult(StageResult):
    """The job and output of `SegdupStage`"""

    job: Job
    output_dir: pathlib.Path


@dataclass(kw_only=True)
class SegdupStage(Stage):
    """Call variants in difficult segmental duplications.

    segdup-caller is sex-aware, so callers run this stage after ploidy
    estimation. Long reads are optional and need both ``lr_alignment``
    and ``lr_bundle``. ``genes`` is a comma-separated gene list; without
    one the caller's default gene set is used. ``overrides`` become
    ``--set KEY=VALUE`` arguments.
    """

    sr_alignment: pathlib.Path
    sr_bundle: pathlib.Path
    output_dir: pathlib.Path
    input_vcf: Optional[pathlib.Path] = None
    sample_sex: Optional[SampleSex] = None
    genes: Optional[str] = None
    lr_alignment: Optional[pathlib.Path] = None
    lr_bundle: Optional[pathlib.Path] = None
    overrides: Sequence[str] = ()
    name: str = "segdup-caller"
    task_name: str = "segdup"

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> SegdupResult:
        deps = set(upstream)

        job = Job(
            cmds.cmd_segdup_caller(
                self.output_dir,
                self.sr_alignment,
                reference=self.ctx.reference,
                sr_bundle=self.sr_bundle,
                input_vcf=self.input_vcf,
                sex=caller_sex_arg(self.sample_sex),
                genes=self.genes,
                overrides=list(self.overrides),
                lr_alignments=self.lr_alignment,
                lr_bundle=self.lr_bundle,
            ),
            self.name,
            self.ctx.cores,
            task_name=self.task_name,
        )
        dag.add_job(job, deps)

        return SegdupResult(
            jobs=[job],
            terminal={job},
            job=job,
            output_dir=self.output_dir,
        )
