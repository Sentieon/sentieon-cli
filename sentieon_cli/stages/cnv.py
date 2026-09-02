"""
CNV calling with CNVscope
"""

from dataclasses import dataclass
import pathlib
from typing import Iterable, List, Optional

import packaging.version

from ..dag import DAG
from ..driver import CNVModelApply, CNVscope
from ..job import Job
from ..util import SampleSex, cnvscope_sex_args
from .base import Stage, StageResult, driver_job

CNV_MIN_VERSIONS = {
    # 202503.04 adds the CNVscope `--sex` and `--par` arguments
    "sentieon driver": packaging.version.Version("202503.04"),
}


@dataclass(kw_only=True)
class CNVResult(StageResult):
    """The jobs and output of `CNVscopeStage`"""

    cnvscope_job: Job
    apply_job: Job
    cnv_vcf: pathlib.Path


@dataclass(kw_only=True)
class CNVscopeStage(Stage):
    """Call CNVs with CNVscope, then filter them with CNVModelApply.

    The sample sex has to be known before the jobs are built, so callers
    run this stage after ploidy estimation. Version checks stay with the
    pipelines, which run them during validation.
    """

    inputs: List[pathlib.Path]
    model: pathlib.Path
    cnvscope_vcf: pathlib.Path
    cnv_vcf: pathlib.Path
    sample_sex: Optional[SampleSex] = None
    par_bed: Optional[pathlib.Path] = None
    interval: Optional[pathlib.Path] = None
    replace_rg: Optional[List[List[str]]] = None
    name: str = "cnvscope"
    apply_name: str = "cnv-model-apply"
    task_name: str = "cnv"

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> CNVResult:
        deps = set(upstream)
        sex, par = cnvscope_sex_args(self.sample_sex, self.par_bed)

        cnvscope = driver_job(
            self.ctx,
            [
                CNVscope(
                    self.cnvscope_vcf,
                    self.model,
                    sex=sex,
                    par=par,
                )
            ],
            inputs=self.inputs,
            interval=self.interval,
            replace_rg=self.replace_rg,
            name=self.name,
            task_name=self.task_name,
        )
        apply_job = driver_job(
            self.ctx,
            [
                CNVModelApply(
                    self.cnv_vcf,
                    self.model,
                    vcf=self.cnvscope_vcf,
                )
            ],
            name=self.apply_name,
            task_name=self.task_name,
        )

        dag.add_job(cnvscope, deps)
        dag.add_job(apply_job, {cnvscope})

        return CNVResult(
            jobs=[cnvscope, apply_job],
            terminal={apply_job},
            cnvscope_job=cnvscope,
            apply_job=apply_job,
            cnv_vcf=self.cnv_vcf,
        )
