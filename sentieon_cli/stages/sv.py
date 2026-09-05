"""
Structural variant calling with LongReadSV
"""

from dataclasses import dataclass
import pathlib
from typing import Iterable, List, Optional

import packaging.version

from ..dag import DAG
from ..driver import LongReadSV
from ..job import Job
from .base import Stage, StageResult, driver_job

LONGREADSV_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}


@dataclass(kw_only=True)
class LongReadSVResult(StageResult):
    """The job and output of `LongReadSVStage`"""

    job: Job
    sv_vcf: pathlib.Path


@dataclass(kw_only=True)
class LongReadSVStage(Stage):
    """Call SVs from long reads with Sentieon LongReadSV.

    ``output`` names the SV VCF; when it is None the path is derived from
    the run's output VCF, so the caller only supplies the alignment and
    the model. The remaining fields are LongReadSV algo options, emitted
    only when set.
    """

    inputs: List[pathlib.Path]
    model: pathlib.Path
    output: Optional[pathlib.Path] = None
    min_map_qual: Optional[int] = None
    min_sv_size: Optional[int] = None
    min_dp: Optional[int] = None
    min_af: Optional[float] = None
    interval: Optional[pathlib.Path] = None
    replace_rg: Optional[List[List[str]]] = None
    name: str = "LongReadSV"
    task_name: str = "sv-calling"

    def build(self) -> LongReadSVResult:
        """Construct this stage's job without touching a DAG"""
        sv_vcf = self.output
        if sv_vcf is None:
            sv_vcf = pathlib.Path(
                str(self.ctx.output_vcf).replace(".vcf.gz", ".sv.vcf.gz")
            )
        job = driver_job(
            self.ctx,
            [
                LongReadSV(
                    sv_vcf,
                    model=self.model,
                    min_map_qual=self.min_map_qual,
                    min_sv_size=self.min_sv_size,
                    min_dp=self.min_dp,
                    min_af=self.min_af,
                )
            ],
            inputs=self.inputs,
            interval=self.interval,
            replace_rg=self.replace_rg,
            name=self.name,
            task_name=self.task_name,
        )
        return LongReadSVResult(
            jobs=[job],
            terminal={job},
            job=job,
            sv_vcf=sv_vcf,
        )

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> LongReadSVResult:
        result = self.build()
        dag.add_job(result.job, set(upstream))
        return result
