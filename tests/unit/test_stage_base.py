"""
Unit tests for the pipeline stage building blocks
"""

import argparse
import dataclasses
import pathlib
from typing import Iterable
from unittest.mock import patch

import packaging.version
import pytest

from sentieon_cli import util
from sentieon_cli.dag import DAG
from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.driver import LocusCollector
from sentieon_cli.job import Job
from sentieon_cli.stages.base import (
    Stage,
    StageContext,
    StageResult,
    driver_job,
    rm_job,
)
from sentieon_cli.util import require_versions, versions_available


def create_mock_args(loglevel: str = "WARNING") -> argparse.Namespace:
    """The argparse namespace `setup_logging` needs"""
    return argparse.Namespace(loglevel=loglevel)


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


@dataclasses.dataclass(kw_only=True)
class ToyStage(Stage):
    """Two jobs: `b` depends on `a`, `a` is the only entry job"""

    tag: str = "toy"

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> StageResult:
        deps = set(upstream)
        job_a = rm_job(
            [self.ctx.tmp_dir / f"{self.tag}-a"], f"{self.tag}-first"
        )
        job_b = rm_job(
            [self.ctx.tmp_dir / f"{self.tag}-b"], f"{self.tag}-second"
        )
        dag.add_job(job_a, deps)
        dag.add_job(job_b, {job_a})
        return StageResult(jobs=[job_a, job_b], terminal={job_b})


class TestStageContext:
    """The run-wide settings object"""

    def test_fields(self, tmp_path):
        ctx = make_ctx(tmp_path)
        assert ctx.reference == tmp_path / "ref.fa"
        assert ctx.output_vcf == tmp_path / "output.vcf.gz"
        assert ctx.tmp_dir == tmp_path
        assert ctx.cores == 4
        assert ctx.dry_run is True
        assert ctx.skip_version_check is True

    def test_frozen(self, tmp_path):
        ctx = make_ctx(tmp_path)
        with pytest.raises(dataclasses.FrozenInstanceError):
            ctx.cores = 8  # type: ignore[misc]


class TestStage:
    """The `add_to` contract"""

    def test_upstream_lands_on_entry_jobs_only(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = ToyStage(ctx=make_ctx(tmp_path)).add_to(dag, [upstream])

        job_a, job_b = result.jobs
        assert result.jobs == [job_a, job_b]
        assert result.terminal == {job_b}
        # The entry job gets the upstream dependencies...
        assert dag.waiting_jobs[job_a] == {upstream}
        # ...and the internal edge is added on top of them.
        assert dag.waiting_jobs[job_b] == {job_a}

    def test_no_upstream_makes_the_entry_job_ready(self, tmp_path):
        dag = DAG()
        result = ToyStage(ctx=make_ctx(tmp_path)).add_to(dag)

        job_a, job_b = result.jobs
        assert job_a in dag.ready_jobs
        assert dag.waiting_jobs[job_b] == {job_a}

    def test_upstream_iterable_consumed_once(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = ToyStage(ctx=make_ctx(tmp_path)).add_to(dag, iter([upstream]))
        assert dag.waiting_jobs[result.jobs[0]] == {upstream}

    def test_add_to_twice_raises(self, tmp_path):
        dag = DAG()
        stage = ToyStage(ctx=make_ctx(tmp_path))
        stage.add_to(dag)
        # Job identity is the shell pipeline, so the same jobs collide
        with pytest.raises(ValueError):
            stage.add_to(dag)

    def test_constructing_a_stage_adds_nothing(self, tmp_path):
        dag = DAG()
        ToyStage(ctx=make_ctx(tmp_path))
        assert not dag.waiting_jobs
        assert not dag.ready_jobs


class TestDriverJob:
    """The `sentieon driver` job helper"""

    def test_basic_command(self, tmp_path):
        ctx = make_ctx(tmp_path, cores=8)
        job = driver_job(
            ctx,
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="locuscollector",
            task_name="dedup",
            inputs=[tmp_path / "sample.bam"],
        )
        expected = (
            f"sentieon driver --input {tmp_path}/sample.bam "
            f"--reference {tmp_path}/ref.fa --thread_count 8 "
            f"--algo LocusCollector {tmp_path}/score.txt.gz"
        )
        assert str(job.shell) == expected
        assert job.name == "locuscollector"
        assert job.task_name == "dedup"
        assert job.threads == 8

    def test_no_interval_padding_by_default(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="locuscollector",
            task_name="dedup",
        )
        assert "--interval_padding" not in str(job.shell)
        assert "--interval" not in str(job.shell)
        assert "--read_filter" not in str(job.shell)
        assert "--input" not in str(job.shell)

    def test_interval_and_padding(self, tmp_path):
        bed = tmp_path / "regions.bed"
        job = driver_job(
            make_ctx(tmp_path),
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="locuscollector",
            task_name="dedup",
            interval=bed,
            interval_padding=5,
        )
        assert f"--interval {bed} --interval_padding 5" in str(job.shell)

    def test_interval_as_a_string(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="readwriter",
            task_name="extract",
            interval="chr6:28510020-33480577",
        )
        assert "--interval chr6:28510020-33480577" in str(job.shell)

    def test_threads_zero(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="locuscollector",
            task_name="dedup",
            threads=0,
        )
        # `threads` is the scheduling budget only; the driver still gets
        # the run's core count.
        assert job.threads == 0
        assert "--thread_count 4" in str(job.shell)

    def test_read_filter_repeats_the_flag(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="locuscollector",
            task_name="dedup",
            read_filter=["QualCalFilter,table=t", "IndelLeftAlign,rgid=rg"],
        )
        assert (
            "--read_filter QualCalFilter,table=t "
            "--read_filter IndelLeftAlign,rgid=rg" in str(job.shell)
        )

    def test_replace_rg_interleaves_with_input(self, tmp_path):
        bam1 = tmp_path / "one.bam"
        bam2 = tmp_path / "two.bam"
        job = driver_job(
            make_ctx(tmp_path),
            [LocusCollector(tmp_path / "score.txt.gz")],
            name="locuscollector",
            task_name="dedup",
            inputs=[bam1, bam2],
            replace_rg=[["rg1=@RG\\tID:rg1"], ["rg2=@RG\\tID:rg2"]],
        )
        shell = str(job.shell)
        assert (
            f"--replace_rg 'rg1=@RG\\tID:rg1' --input {bam1} "
            f"--replace_rg 'rg2=@RG\\tID:rg2' --input {bam2}" in shell
        )

    def test_algos_are_added_in_order(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [
                LocusCollector(tmp_path / "first.txt.gz"),
                LocusCollector(tmp_path / "second.txt.gz"),
            ],
            name="locuscollector",
            task_name="dedup",
        )
        shell = str(job.shell)
        assert shell.index("first.txt.gz") < shell.index("second.txt.gz")


class TestRmJob:
    """The cleanup job helper"""

    def test_command_and_metadata(self, tmp_path):
        vcf = tmp_path / "tmp.vcf.gz"
        job = rm_job([vcf, str(vcf) + ".tbi"], "rm-tmp-vcf")
        assert str(job.shell) == f"rm {vcf} {vcf}.tbi"
        assert job.name == "rm-tmp-vcf"
        assert job.task_name == "cleanup"
        assert job.threads == 0
        assert job.shell.nodes[0].fail_ok is True


class TestRequireVersions:
    """`require_versions` and `versions_available`"""

    min_versions = {"sentieon driver": packaging.version.Version("202308")}

    def test_require_skips(self):
        with patch.object(util, "check_version") as check:
            require_versions(self.min_versions, skip=True)
        check.assert_not_called()

    def test_require_passes(self):
        with patch.object(util, "check_version", return_value=True) as check:
            require_versions(self.min_versions)
        check.assert_called_once_with(
            "sentieon driver", packaging.version.Version("202308")
        )

    def test_require_exits_on_failure(self):
        with patch.object(util, "check_version", return_value=False):
            with pytest.raises(SystemExit) as excinfo:
                require_versions(self.min_versions)
        assert excinfo.value.code == 2

    def test_require_exits_on_the_first_failure(self):
        versions = {
            "samtools": packaging.version.Version("1.16"),
            "bcftools": packaging.version.Version("1.10"),
        }
        with patch.object(util, "check_version", return_value=False) as check:
            with pytest.raises(SystemExit):
                require_versions(versions)
        assert check.call_count == 1

    def test_available_skips(self):
        with patch.object(util, "check_version") as check:
            assert versions_available(self.min_versions, skip=True) is True
        check.assert_not_called()

    def test_available_true(self):
        with patch.object(util, "check_version", return_value=True):
            assert versions_available(self.min_versions) is True

    def test_available_false(self):
        with patch.object(util, "check_version", return_value=False):
            assert versions_available(self.min_versions) is False


class TestPipelineStageContext:
    """`BasePipeline.stage_context` and `BasePipeline.required`"""

    def make_pipeline(self, tmp_path) -> DNAscopePipeline:
        pipeline = DNAscopePipeline()
        pipeline.setup_logging(create_mock_args())
        pipeline.reference = tmp_path / "ref.fa"
        pipeline.output_vcf = tmp_path / "output.vcf.gz"
        pipeline.cores = 3
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = tmp_path
        return pipeline

    def test_stage_context(self, tmp_path):
        ctx = self.make_pipeline(tmp_path).stage_context()
        assert ctx == StageContext(
            reference=tmp_path / "ref.fa",
            output_vcf=tmp_path / "output.vcf.gz",
            tmp_dir=tmp_path,
            cores=3,
            dry_run=True,
            skip_version_check=True,
        )

    def test_missing_reference(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path)
        pipeline.reference = None
        with pytest.raises(SystemExit) as excinfo:
            pipeline.stage_context()
        assert excinfo.value.code == 2

    def test_missing_output_vcf(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path)
        pipeline.output_vcf = None
        with pytest.raises(SystemExit) as excinfo:
            pipeline.stage_context()
        assert excinfo.value.code == 2

    def test_missing_tmp_dir(self, tmp_path):
        pipeline = DNAscopePipeline()
        pipeline.setup_logging(create_mock_args())
        pipeline.reference = tmp_path / "ref.fa"
        pipeline.output_vcf = tmp_path / "output.vcf.gz"
        assert not hasattr(pipeline, "tmp_dir")
        with pytest.raises(SystemExit) as excinfo:
            pipeline.stage_context()
        assert excinfo.value.code == 2

    def test_required_returns_the_value(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path)
        assert pipeline.required(pipeline.output_vcf, "output_vcf") == (
            tmp_path / "output.vcf.gz"
        )

    def test_required_exits_when_missing(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path)
        with pytest.raises(SystemExit) as excinfo:
            pipeline.required(None, "model_bundle")
        assert excinfo.value.code == 2
