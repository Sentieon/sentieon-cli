"""
Unit tests for the CNVscope stage
"""

import dataclasses
import importlib.util
from importlib.resources import files
import json
import logging
import pathlib
import sys
from typing import List

import packaging.version
import pytest

from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.cnv import (
    CNV_MIN_VERSIONS,
    CNVscopeStage,
    PangenomeCNVStage,
)
from sentieon_cli.util import SampleSex
from tests.utils.test_helpers import write_bundle

PACKAGE_LOGGER = "sentieon_cli"


class _RecordingHandler(logging.Handler):
    """Collects formatted messages from the package logger."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


@pytest.fixture
def messages():
    """`caplog` does not see the package logger, which does not propagate"""
    logger = logging.getLogger(PACKAGE_LOGGER)
    handler = _RecordingHandler()
    logger.addHandler(handler)
    try:
        yield handler.messages
    finally:
        logger.removeHandler(handler)


def make_ctx(tmp_path: pathlib.Path, cores: int = 4) -> StageContext:
    """A StageContext over a temporary directory"""
    return StageContext(
        reference=tmp_path / "ref.fa",
        output_vcf=tmp_path / "output.vcf.gz",
        tmp_dir=tmp_path,
        cores=cores,
        dry_run=True,
        skip_version_check=True,
    )


def make_stage(tmp_path: pathlib.Path, **kwargs) -> CNVscopeStage:
    """A CNVscopeStage with the paths a caller supplies"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.cram"],
        model=tmp_path / "bundle" / "cnv.model",
        cnvscope_vcf=tmp_path / "cnvscope.vcf.gz",
        cnv_vcf=tmp_path / "output.cnv.vcf.gz",
    )
    defaults.update(kwargs)
    return CNVscopeStage(**defaults)  # type: ignore[arg-type]


class TestMinVersions:
    """The version constant moved here with the stage"""

    def test_cnvscope_sex_arguments_need_202503_04(self):
        assert CNV_MIN_VERSIONS["sentieon driver"] == (
            packaging.version.Version("202503.04")
        )


class TestCommands:
    """The `sentieon driver` commands the stage builds"""

    def test_cnvscope_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert str(result.cnvscope_job.shell) == (
            f"sentieon driver --input {tmp_path}/sample.cram "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            f"--algo CNVscope --model {tmp_path}/bundle/cnv.model "
            f"--blocks {tmp_path}/cnvscope-blocks.bed "
            f"{tmp_path}/cnvscope.vcf.gz"
        )

    def test_cnv_model_apply_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        # CNVModelApply reads the CNVscope VCF; it takes no alignment.
        # `--emit_conf 0` writes the rejected segments too
        assert str(result.apply_job.shell) == (
            f"sentieon driver --reference {tmp_path}/ref.fa "
            f"--thread_count 4 --algo CNVModelApply "
            f"--model {tmp_path}/bundle/cnv.model "
            f"--vcf {tmp_path}/cnvscope.vcf.gz "
            "--emit_conf 0 "
            f"{tmp_path}/output.cnv.vcf.gz"
        )

    def test_interval_is_passed_to_cnvscope_only(self, tmp_path):
        bed = tmp_path / "regions.bed"
        result = make_stage(tmp_path, interval=bed).add_to(DAG())

        assert f"--interval {bed}" in str(result.cnvscope_job.shell)
        assert "--interval" not in str(result.apply_job.shell)

    def test_replace_rg_is_passed_to_cnvscope_only(self, tmp_path):
        result = make_stage(
            tmp_path,
            replace_rg=[["rg1=ID:rg1\\tSM:sample"]],
        ).add_to(DAG())

        assert "--replace_rg 'rg1=ID:rg1\\tSM:sample'" in str(
            result.cnvscope_job.shell
        )
        assert "--replace_rg" not in str(result.apply_job.shell)

    def test_dump_probes_is_passed_to_cnvscope_only(self, tmp_path):
        probes = tmp_path / "output.cnv.probes"
        result = make_stage(tmp_path, dump_probes=probes).add_to(DAG())

        assert f"--dump_probes {probes}" in str(result.cnvscope_job.shell)
        assert "--dump_probes" not in str(result.apply_job.shell)

    def test_multiple_inputs(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.bam", tmp_path / "two.bam"],
        ).add_to(DAG())

        shell = str(result.cnvscope_job.shell)
        assert f"--input {tmp_path}/one.bam" in shell
        assert f"--input {tmp_path}/two.bam" in shell


class TestSexArguments:
    """`--sex` and `--par` follow the sample sex"""

    def test_male_sample(self, tmp_path):
        par_bed = tmp_path / "par.bed"
        result = make_stage(
            tmp_path, sample_sex=SampleSex.MALE, par_bed=par_bed
        ).add_to(DAG())

        shell = str(result.cnvscope_job.shell)
        assert "--sex M" in shell
        assert f"--par {par_bed}" in shell
        # The result records the ploidy CNVscope was called with
        assert (result.sex, result.par) == ("M", par_bed)

    def test_female_sample(self, tmp_path):
        # A female sample does not need the PAR regions
        result = make_stage(
            tmp_path,
            sample_sex=SampleSex.FEMALE,
            par_bed=tmp_path / "par.bed",
        ).add_to(DAG())

        shell = str(result.cnvscope_job.shell)
        assert "--sex F" in shell
        assert "--par" not in shell

    @pytest.mark.parametrize("sample_sex", [None, SampleSex.UNKNOWN])
    def test_without_a_known_sex(self, tmp_path, sample_sex, messages):
        result = make_stage(tmp_path, sample_sex=sample_sex).add_to(DAG())

        shell = str(result.cnvscope_job.shell)
        blocks = tmp_path / "cnvscope-blocks.bed"
        assert "--sex" not in shell
        assert "--par" not in shell
        assert f"--blocks {blocks}" in shell
        assert result.blocks == blocks
        # A dry run does not write the file
        assert not blocks.exists()
        assert any("Sex estimation failed" in msg for msg in messages)

    def test_without_a_known_sex_writes_the_blocks(self, tmp_path):
        (tmp_path / "ref.fa.fai").write_text(
            "chr20\t64444167\t7\t60\t61\n"
            "chrX\t156040895\t65518245\t60\t61\n"
            "chrY\t57227415\t224159863\t60\t61\n"
        )
        par_bed = tmp_path / "par.bed"
        par_bed.write_text(
            "chrX\t10000\t2781479\nchrY\t10000\t2781479\n"
            "chrX\t155701382\t156030895\nchrY\t56887902\t57217415\n"
        )
        ctx = dataclasses.replace(make_ctx(tmp_path), dry_run=False)

        result = make_stage(
            tmp_path,
            ctx=ctx,
            sample_sex=SampleSex.UNKNOWN,
            par_bed=par_bed,
        ).add_to(DAG())

        assert result.blocks is not None
        assert result.blocks.read_text() == (
            "chrX\t0\t156040895\t2\n"
            "chrY\t0\t10000\t1\n"
            "chrY\t2781479\t56887902\t1\n"
            "chrY\t57217415\t57227415\t1\n"
        )


class TestJobMetadata:
    """Job names, task name and thread counts"""

    def test_default_names(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.cnvscope_job.name == "cnvscope"
        assert result.apply_job.name == "cnv-model-apply"

    def test_overridden_names(self, tmp_path):
        result = make_stage(
            tmp_path, name="CNVscope", apply_name="CNVModelApply"
        ).add_to(DAG())

        assert result.cnvscope_job.name == "CNVscope"
        assert result.apply_job.name == "CNVModelApply"

    def test_task_name_and_threads(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.jobs:
            assert job.task_name == "cnv"
            assert job.threads == 4


class TestDagWiring:
    """The edges the stage inserts"""

    def test_upstream_lands_on_cnvscope_only(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.cnvscope_job] == {upstream}
        assert dag.waiting_jobs[result.apply_job] == {result.cnvscope_job}

    def test_cnvscope_is_a_root_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert result.cnvscope_job in dag.ready_jobs
        assert dag.waiting_jobs[result.apply_job] == {result.cnvscope_job}

    def test_result_fields(self, tmp_path):
        cnv_vcf = tmp_path / "output.cnv.vcf.gz"
        result = make_stage(tmp_path, cnv_vcf=cnv_vcf).add_to(DAG())

        assert result.jobs == [result.cnvscope_job, result.apply_job]
        assert result.terminal == {result.apply_job}
        assert result.cnv_vcf == cnv_vcf


def make_pangenome_stage(
    tmp_path: pathlib.Path, **kwargs
) -> PangenomeCNVStage:
    """A PangenomeCNVStage over the short-read alignment and the SV VCF"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.cram"],
        model=tmp_path / "bundle" / "cnv.model",
        sv_vcf=tmp_path / "output_sv.vcf.gz",
    )
    defaults.update(kwargs)
    return PangenomeCNVStage(**defaults)  # type: ignore[arg-type]


def script_path(name: str) -> str:
    """The packaged helper script, unresolved as the stage builds it"""
    return str(files("sentieon_cli.scripts").joinpath(name))


class TestPangenomeCNVStage:
    """CNVscope, CNVModelApply, indel2cnv and combine_cnv together"""

    def test_job_names_threads_and_task(self, tmp_path):
        result = make_pangenome_stage(tmp_path).add_to(DAG())

        assert [job.name for job in result.jobs] == [
            "cnvscope",
            "cnv-model-apply",
            "indel2cnv",
            "combine-cnv",
        ]
        assert [job.threads for job in result.jobs] == [4, 4, 0, 0]
        for job in result.jobs:
            assert job.task_name == "cnv"

    def test_derived_paths(self, tmp_path):
        result = make_pangenome_stage(tmp_path).add_to(DAG())

        # The intermediates live in the temporary directory; the CNV VCF
        # and the probes file sit next to the output VCF
        cnvscope_shell = str(result.cnvscope_job.shell)
        assert f"{tmp_path}/sample-cnvscope.vcf.gz" in cnvscope_shell
        assert f"--dump_probes {tmp_path}/output_cnv.probes" in cnvscope_shell

        apply_shell = str(result.apply_job.shell)
        assert f"--vcf {tmp_path}/sample-cnvscope.vcf.gz" in apply_shell
        assert f"{tmp_path}/sample-cnv_model_apply.vcf.gz" in apply_shell

        assert result.cnv_vcf == tmp_path / "output_cnv.vcf.gz"

    def test_cnv_model_apply_writes_every_segment(self, tmp_path):
        """combine_cnv.py reads the rejected segments from CNVModelApply"""
        result = make_pangenome_stage(tmp_path).add_to(DAG())

        assert "--emit_conf 0 " in str(result.apply_job.shell)
        # So it needs no raw CNVscope VCF
        combine_shell = str(result.combine_job.shell)
        assert "--raw" not in combine_shell
        assert "sample-cnvscope.vcf.gz" not in combine_shell

    def test_indel2cnv_command(self, tmp_path):
        result = make_pangenome_stage(tmp_path).add_to(DAG())

        assert str(result.indel2cnv_job.shell) == (
            f"{sys.executable} {script_path('indel2cnv.py')} "
            f"{tmp_path}/ref.fa {tmp_path}/output_sv.vcf.gz "
            f"{tmp_path}/sample-sv_cnv.vcf.gz -t 4"
        )

    def test_combine_cnv_command(self, tmp_path):
        """Without a preset combine_cnv.py uses its defaults"""
        result = make_pangenome_stage(tmp_path).add_to(DAG())

        assert str(result.combine_job.shell) == (
            f"{sys.executable} {script_path('combine_cnv.py')} "
            f"--cnv {tmp_path}/sample-cnv_model_apply.vcf.gz "
            f"--converted {tmp_path}/sample-sv_cnv.vcf.gz "
            f"--blocks {tmp_path}/cnvscope-blocks.bed "
            f"-o {tmp_path}/output_cnv.vcf.gz"
        )

    def test_combine_cnv_preset(self, tmp_path):
        """The preset, a model bundle member, is passed as it is"""
        preset = tmp_path / "model.bundle" / "combine_cnv.json"
        stage = make_pangenome_stage(tmp_path, combine_preset=preset)
        result = stage.add_to(DAG())

        assert str(result.combine_job.shell) == (
            f"{sys.executable} {script_path('combine_cnv.py')} "
            f"--cnv {tmp_path}/sample-cnv_model_apply.vcf.gz "
            f"--converted {tmp_path}/sample-sv_cnv.vcf.gz "
            f"--preset {preset} "
            f"--blocks {tmp_path}/cnvscope-blocks.bed "
            f"-o {tmp_path}/output_cnv.vcf.gz"
        )

    def test_combine_cnv_male_sample(self, tmp_path):
        """combine_cnv.py gets the same sex and PAR BED file as CNVscope"""
        par_bed = tmp_path / "par.bed"
        result = make_pangenome_stage(
            tmp_path, sample_sex=SampleSex.MALE, par_bed=par_bed
        ).add_to(DAG())

        assert str(result.combine_job.shell) == (
            f"{sys.executable} {script_path('combine_cnv.py')} "
            f"--cnv {tmp_path}/sample-cnv_model_apply.vcf.gz "
            f"--converted {tmp_path}/sample-sv_cnv.vcf.gz "
            f"--sex M --par {par_bed} "
            f"-o {tmp_path}/output_cnv.vcf.gz"
        )

    def test_combine_cnv_female_sample(self, tmp_path):
        result = make_pangenome_stage(
            tmp_path,
            sample_sex=SampleSex.FEMALE,
            par_bed=tmp_path / "par.bed",
        ).add_to(DAG())

        shell = str(result.combine_job.shell)
        assert "--sex F" in shell
        assert "--par" not in shell

    def test_combine_cnv_unknown_sex(self, tmp_path):
        """Without a known sex both tools get the same assumed ploidies"""
        result = make_pangenome_stage(
            tmp_path,
            sample_sex=SampleSex.UNKNOWN,
            par_bed=tmp_path / "par.bed",
        ).add_to(DAG())

        blocks = tmp_path / "cnvscope-blocks.bed"
        for job in (result.cnvscope_job, result.combine_job):
            shell = str(job.shell)
            assert "--sex" not in shell
            assert "--par" not in shell
            assert f"--blocks {blocks}" in shell

    def test_upstream_lands_on_the_entry_jobs(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_pangenome_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.cnvscope_job] == {upstream}
        assert dag.waiting_jobs[result.indel2cnv_job] == {upstream}
        assert dag.waiting_jobs[result.apply_job] == {result.cnvscope_job}
        assert dag.waiting_jobs[result.combine_job] == {
            result.apply_job,
            result.indel2cnv_job,
        }

    def test_entry_jobs_are_roots_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_pangenome_stage(tmp_path).add_to(dag)

        assert result.cnvscope_job in dag.ready_jobs
        assert result.indel2cnv_job in dag.ready_jobs

    def test_interval_and_replace_rg_reach_cnvscope_only(self, tmp_path):
        bed = tmp_path / "regions.bed"
        result = make_pangenome_stage(
            tmp_path,
            interval=bed,
            replace_rg=[["rg1=ID:rg1\\tSM:sample"]],
        ).add_to(DAG())

        cnvscope_shell = str(result.cnvscope_job.shell)
        assert f"--interval {bed}" in cnvscope_shell
        assert "--replace_rg 'rg1=ID:rg1\\tSM:sample'" in cnvscope_shell
        for job in (
            result.apply_job,
            result.indel2cnv_job,
            result.combine_job,
        ):
            assert "--interval" not in str(job.shell)
            assert "--replace_rg" not in str(job.shell)

    def test_sex_arguments_reach_cnvscope(self, tmp_path):
        par_bed = tmp_path / "par.bed"
        result = make_pangenome_stage(
            tmp_path, sample_sex=SampleSex.MALE, par_bed=par_bed
        ).add_to(DAG())

        shell = str(result.cnvscope_job.shell)
        assert "--sex M" in shell
        assert f"--par {par_bed}" in shell

    def test_female_sample_needs_no_par(self, tmp_path):
        result = make_pangenome_stage(
            tmp_path,
            sample_sex=SampleSex.FEMALE,
            par_bed=tmp_path / "par.bed",
        ).add_to(DAG())

        shell = str(result.cnvscope_job.shell)
        assert "--sex F" in shell
        assert "--par" not in shell

    def test_result_fields(self, tmp_path):
        result = make_pangenome_stage(tmp_path).add_to(DAG())

        assert result.jobs == [
            result.cnvscope_job,
            result.apply_job,
            result.indel2cnv_job,
            result.combine_job,
        ]
        assert result.terminal == {result.combine_job}
        assert result.cnv_vcf == tmp_path / "output_cnv.vcf.gz"


def load_combine_cnv():
    """The packaged combine_cnv.py, imported as a module"""
    spec = importlib.util.spec_from_file_location(
        "combine_cnv", script_path("combine_cnv.py")
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestCombinePreset:
    """combine_cnv.py reads the preset the pipelines pass"""

    def test_bundle_preset_loads_with_combine_cnv(self, tmp_path):
        """`--preset <bundle>/combine_cnv.json` names a bundle member"""
        preset = {
            "description": "Single-end CNVscope models",
            "max_converted_loss": 3000,
            "loss_support_size": 10000,
        }
        bundle = write_bundle(
            tmp_path / "model.bundle",
            {
                "cnv.model": b"",
                "combine_cnv.json": json.dumps(preset).encode(),
            },
        )
        combine_cnv = load_combine_cnv()

        rules = combine_cnv.load_preset(bundle / "combine_cnv.json")
        assert rules == dict(
            combine_cnv.DEFAULTS,
            max_converted_loss=3000,
            loss_support_size=10000,
        )
