"""
Unit tests for the segdup-caller stage
"""

import pathlib

import packaging.version
import pytest

from sentieon_cli import command_strings as cmds
from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.segdup import (
    SEGDUP_MIN_VERSIONS,
    SegdupStage,
)
from sentieon_cli.util import SampleSex


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


def make_stage(tmp_path: pathlib.Path, **kwargs) -> SegdupStage:
    """A SegdupStage with the paths a caller supplies"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        sr_alignment=tmp_path / "sample.cram",
        sr_bundle=tmp_path / "bundle",
        output_dir=tmp_path / "output_segdups",
        input_vcf=tmp_path / "output.vcf.gz",
        sample_sex=SampleSex.MALE,
    )
    defaults.update(kwargs)
    return SegdupStage(**defaults)  # type: ignore[arg-type]


class TestMinVersions:
    """The version constant lives next to the stage"""

    def test_set_argument_needs_0_7_0(self):
        assert SEGDUP_MIN_VERSIONS["segdup-caller"] == (
            packaging.version.Version("0.7.0")
        )


class TestCommands:
    """The `segdup-caller` commands the stage builds"""

    def test_short_read_command(self, tmp_path):
        # The command dnascope-pangenome builds today
        result = make_stage(tmp_path, genes="CFH,CYP2D6,SMN1").add_to(DAG())

        assert str(result.job.shell) == (
            f"segdup-caller --short {tmp_path}/sample.cram "
            f"--reference {tmp_path}/ref.fa "
            f"--sr_model {tmp_path}/bundle "
            f"--input_vcf {tmp_path}/output.vcf.gz "
            "--sex male --genes CFH,CYP2D6,SMN1 "
            f"--outdir {tmp_path}/output_segdups"
        )

    def test_long_read_command(self, tmp_path):
        result = make_stage(
            tmp_path,
            lr_alignment=tmp_path / "long.cram",
            lr_bundle=tmp_path / "bundle",
        ).add_to(DAG())

        assert str(result.job.shell) == (
            f"segdup-caller --short {tmp_path}/sample.cram "
            f"--long {tmp_path}/long.cram "
            f"--reference {tmp_path}/ref.fa "
            f"--sr_model {tmp_path}/bundle "
            f"--lr_model {tmp_path}/bundle "
            f"--input_vcf {tmp_path}/output.vcf.gz "
            "--sex male "
            f"--outdir {tmp_path}/output_segdups"
        )

    def test_overrides_become_set_arguments(self, tmp_path):
        result = make_stage(
            tmp_path,
            overrides=["main.min_map_qual=30", "main.other=1"],
        ).add_to(DAG())

        assert "--set main.min_map_qual=30 --set main.other=1 --outdir" in str(
            result.job.shell
        )

    def test_without_an_input_vcf_or_genes(self, tmp_path):
        result = make_stage(tmp_path, input_vcf=None).add_to(DAG())

        shell = str(result.job.shell)
        assert "--input_vcf" not in shell
        assert "--genes" not in shell
        assert "--set" not in shell

    @pytest.mark.parametrize(
        "sample_sex,expected",
        [
            (SampleSex.MALE, "male"),
            (SampleSex.FEMALE, "female"),
            (SampleSex.UNKNOWN, "female"),
            (None, "female"),
        ],
    )
    def test_sex_argument(self, tmp_path, sample_sex, expected):
        result = make_stage(tmp_path, sample_sex=sample_sex).add_to(DAG())

        assert f"--sex {expected} " in str(result.job.shell)


class TestLongReadValidation:
    """Long-read calling needs both the alignment and the model bundle"""

    def test_alignment_without_a_bundle(self, tmp_path):
        with pytest.raises(ValueError):
            cmds.cmd_segdup_caller(
                tmp_path / "out",
                tmp_path / "sample.cram",
                reference=tmp_path / "ref.fa",
                sr_bundle=tmp_path / "bundle",
                lr_alignments=tmp_path / "long.cram",
            )

    def test_bundle_without_an_alignment(self, tmp_path):
        with pytest.raises(ValueError):
            cmds.cmd_segdup_caller(
                tmp_path / "out",
                tmp_path / "sample.cram",
                reference=tmp_path / "ref.fa",
                sr_bundle=tmp_path / "bundle",
                lr_bundle=tmp_path / "bundle",
            )

    def test_the_stage_propagates_the_error(self, tmp_path):
        with pytest.raises(ValueError):
            make_stage(tmp_path, lr_bundle=tmp_path / "bundle").add_to(DAG())


class TestJobMetadata:
    """Job name, task name and thread count"""

    def test_default_names_and_threads(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.job.name == "segdup-caller"
        assert result.job.task_name == "segdup"
        assert result.job.threads == 4

    def test_overridden_names(self, tmp_path):
        result = make_stage(
            tmp_path, name="segdups", task_name="segdups"
        ).add_to(DAG())

        assert result.job.name == "segdups"
        assert result.job.task_name == "segdups"


class TestDagWiring:
    """The edges the stage inserts"""

    def test_upstream_becomes_the_dependencies(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.job] == {upstream}

    def test_root_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert result.job in dag.ready_jobs

    def test_result_fields(self, tmp_path):
        output_dir = tmp_path / "output_segdups"
        result = make_stage(tmp_path, output_dir=output_dir).add_to(DAG())

        assert result.jobs == [result.job]
        assert result.terminal == {result.job}
        assert result.output_dir == output_dir
