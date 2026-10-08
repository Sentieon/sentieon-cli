"""
CNV calling with CNVscope
"""

from dataclasses import dataclass
from importlib.resources import files
import pathlib
from typing import Iterable, List, Optional

import packaging.version

from .. import command_strings as cmds
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
    # The CNVscope `--sex` and `--par` arguments, None when omitted
    sex: Optional[str] = None
    par: Optional[pathlib.Path] = None


@dataclass(kw_only=True)
class CNVscopeStage(Stage):
    """Call CNVs with CNVscope, then filter them with CNVModelApply.

    The sample sex has to be known before the jobs are built, so callers
    run this stage after ploidy estimation. Version checks stay with the
    pipelines, which run them during validation. ``dump_probes`` is the
    path CNVscope writes its per-probe debugging output to.
    """

    inputs: List[pathlib.Path]
    model: pathlib.Path
    cnvscope_vcf: pathlib.Path
    cnv_vcf: pathlib.Path
    sample_sex: Optional[SampleSex] = None
    par_bed: Optional[pathlib.Path] = None
    dump_probes: Optional[pathlib.Path] = None
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
                    dump_probes=self.dump_probes,
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
            sex=sex,
            par=par,
        )


@dataclass(kw_only=True)
class PangenomeCNVResult(StageResult):
    """The jobs and output of `PangenomeCNVStage`"""

    cnvscope_job: Job
    apply_job: Job
    indel2cnv_job: Job
    combine_job: Job
    cnv_vcf: pathlib.Path


@dataclass(kw_only=True)
class PangenomeCNVStage(Stage):
    """CNV calling for the pangenome pipelines.

    CNVscope and CNVModelApply run on the short reads, the PangenomeSV
    INDELs are converted to CNVs with `indel2cnv.py`, and the two call
    sets are combined by `combine_cnv.py` into `<output>_cnv.vcf.gz`.
    Both inputs -- the alignments and the SV VCF -- come from an earlier
    DAG, and the sample sex has to be known before the jobs are built.
    `combine_cnv.py` gets the same sex and PAR BED file as CNVscope and
    the CNVscope output as its raw segments. `combine_preset` selects its
    rules: `PE` for paired-end CNVscope models and `SE` for single-end
    models.
    """

    inputs: List[pathlib.Path]
    model: pathlib.Path
    sv_vcf: pathlib.Path
    sample_sex: Optional[SampleSex] = None
    par_bed: Optional[pathlib.Path] = None
    combine_preset: str = "PE"
    interval: Optional[pathlib.Path] = None
    replace_rg: Optional[List[List[str]]] = None
    task_name: str = "cnv"

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> PangenomeCNVResult:
        deps = set(upstream)

        cnvscope_vcf = self.ctx.tmp_dir.joinpath("sample-cnvscope.vcf.gz")
        apply_vcf = self.ctx.tmp_dir.joinpath("sample-cnv_model_apply.vcf.gz")
        converted_vcf = self.ctx.tmp_dir.joinpath("sample-sv_cnv.vcf.gz")
        cnv_vcf = pathlib.Path(
            str(self.ctx.output_vcf).replace(".vcf.gz", "_cnv.vcf.gz")
        )
        cnv_probes = pathlib.Path(
            str(self.ctx.output_vcf).replace(".vcf.gz", "_cnv.probes")
        )

        cnv_result = CNVscopeStage(
            ctx=self.ctx,
            inputs=self.inputs,
            model=self.model,
            cnvscope_vcf=cnvscope_vcf,
            cnv_vcf=apply_vcf,
            sample_sex=self.sample_sex,
            par_bed=self.par_bed,
            dump_probes=cnv_probes,
            interval=self.interval,
            replace_rg=self.replace_rg,
            task_name=self.task_name,
        ).add_to(dag, deps)

        # Convert the PangenomeSV output to CNVs
        indel2cnv_script = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("indel2cnv.py"))
        )
        indel2cnv_job = Job(
            cmds.cmd_pyexec_indel2cnv(
                converted_vcf,
                self.sv_vcf,
                self.ctx.reference,
                indel2cnv_script,
                self.ctx.cores,
            ),
            "indel2cnv",
            0,
            task_name=self.task_name,
        )
        dag.add_job(indel2cnv_job, deps)

        # Combine the CNVModelApply output with the converted SVs. The
        # CNVscope output, written before CNVModelApply runs, is the raw
        # segmentation that `combine_cnv.py` checks long converted gains
        # against
        combine_script = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("combine_cnv.py"))
        )
        combine_job = Job(
            cmds.cmd_pyexec_combine_cnv(
                cnv_vcf,
                apply_vcf,
                converted_vcf,
                combine_script,
                raw_vcf=cnvscope_vcf,
                preset=self.combine_preset,
                sex=cnv_result.sex,
                par=cnv_result.par,
            ),
            "combine-cnv",
            0,
            task_name=self.task_name,
        )
        dag.add_job(combine_job, {cnv_result.apply_job, indel2cnv_job})

        return PangenomeCNVResult(
            jobs=[*cnv_result.jobs, indel2cnv_job, combine_job],
            terminal={combine_job},
            cnvscope_job=cnv_result.cnvscope_job,
            apply_job=cnv_result.apply_job,
            indel2cnv_job=indel2cnv_job,
            combine_job=combine_job,
            cnv_vcf=cnv_vcf,
        )
