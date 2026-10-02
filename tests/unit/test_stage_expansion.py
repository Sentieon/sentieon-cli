"""
Unit tests for the ExpansionHunter stage
"""

import pathlib

import pytest

from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.expansion import (
    EXPANSION_MIN_VERSIONS,
    ExpansionHunterStage,
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


def make_stage(tmp_path: pathlib.Path, **kwargs) -> ExpansionHunterStage:
    """An ExpansionHunterStage with the paths a caller supplies"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        alignment=tmp_path / "sample.cram",
        variant_catalog=tmp_path / "catalog.json",
        output_prefix=tmp_path / "output_expansion",
    )
    defaults.update(kwargs)
    return ExpansionHunterStage(**defaults)  # type: ignore[arg-type]


class TestMinVersions:
    """ExpansionHunter has no minimum version"""

    def test_no_minimum_version(self):
        assert EXPANSION_MIN_VERSIONS == {"ExpansionHunter": None}


class TestCommands:
    """The `ExpansionHunter` command the stage builds"""

    def test_expansion_hunter_command(self, tmp_path):
        result = make_stage(tmp_path, sample_sex=SampleSex.MALE).add_to(DAG())

        assert str(result.job.shell) == (
            f"ExpansionHunter --reads {tmp_path}/sample.cram "
            f"--reference {tmp_path}/ref.fa "
            f"--variant-catalog {tmp_path}/catalog.json "
            "--sex male --threads 4 "
            f"--output-prefix {tmp_path}/output_expansion"
        )

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

    def test_threads_follow_the_run_cores(self, tmp_path):
        result = make_stage(tmp_path, ctx=make_ctx(tmp_path, cores=16)).add_to(
            DAG()
        )

        assert "--threads 16" in str(result.job.shell)


class TestJobMetadata:
    """Job name, task name and thread count"""

    def test_default_names_and_threads(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.job.name == "expansion-hunter"
        assert result.job.task_name == "expansion-hunter"
        assert result.job.threads == 4

    def test_overridden_names(self, tmp_path):
        result = make_stage(
            tmp_path, name="expansions", task_name="expansions"
        ).add_to(DAG())

        assert result.job.name == "expansions"
        assert result.job.task_name == "expansions"


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
        output_prefix = tmp_path / "output_expansion"
        result = make_stage(tmp_path, output_prefix=output_prefix).add_to(
            DAG()
        )

        assert result.jobs == [result.job]
        assert result.terminal == {result.job}
        assert result.output_prefix == output_prefix
