"""
Unit tests for the annotation-transfer stage
"""

import pathlib
from types import SimpleNamespace

import pytest

from sentieon_cli.dag import DAG
from sentieon_cli.shard import Shard
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.transfer import TransferConfig, TransferStage


def make_ctx(tmp_path: pathlib.Path, cores: int = 4) -> StageContext:
    """A StageContext over a temporary directory"""
    return StageContext(
        reference=tmp_path / "ref.fa",
        output_vcf=tmp_path / "output.vcf.gz",
        tmp_dir=tmp_path,
        cores=cores,
        # `build_transfer_jobs` shells out to bcftools unless this is set
        dry_run=True,
        skip_version_check=True,
    )


def make_config(tmp_path: pathlib.Path) -> TransferConfig:
    """Two shards on a contig the pop VCF has, one on a contig it lacks"""
    return TransferConfig(
        pop_vcf=tmp_path / "pop.vcf.gz",
        shards=[
            Shard("chr1", 0, 1000),
            Shard("chr1", 1000, 2000),
            Shard("chrUn_extra", 0, 500),
        ],
        pop_vcf_contigs={"chr1": 2000},
        fai_data={
            "chr1": {"length": 2000},
            "chrUn_extra": {"length": 500},
        },
    )


def make_stage(tmp_path: pathlib.Path, **kwargs) -> TransferStage:
    """A TransferStage over `make_config`"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        config=make_config(tmp_path),
        raw_vcf=tmp_path / "raw.vcf.gz",
        out_vcf=tmp_path / "transfer.vcf.gz",
    )
    defaults.update(kwargs)
    return TransferStage(**defaults)  # type: ignore[arg-type]


class TestJobs:
    """The jobs the fan-out produces"""

    def test_one_job_per_shard_plus_the_extra_contig(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert [job.name for job in result.shard_jobs] == [
            "merge-trim-0",
            "merge-trim-1",
            "merge-trim-extra",
        ]
        assert result.concat_job.name == "merge-trim-concat"

    def test_task_name_and_threads(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.jobs:
            assert job.task_name == "annotation-transfer"
        for job in result.shard_jobs:
            assert job.threads == 1
        # The concat runs with every core
        assert result.concat_job.threads == 4

    def test_shard_jobs_read_the_raw_vcf(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.shard_jobs:
            assert str(tmp_path / "raw.vcf.gz") in str(job.shell)

    def test_concat_writes_the_output(self, tmp_path):
        out_vcf = tmp_path / "annotated.vcf.gz"
        result = make_stage(tmp_path, out_vcf=out_vcf).add_to(DAG())

        assert str(out_vcf) in str(result.concat_job.shell)
        assert result.out_vcf == out_vcf


class TestTag:
    """`tag` disambiguates two fan-outs in one run"""

    def test_tagged_names(self, tmp_path):
        result = make_stage(tmp_path, tag="hp").add_to(DAG())

        assert [job.name for job in result.shard_jobs] == [
            "merge-trim-hp-0",
            "merge-trim-hp-1",
            "merge-trim-hp-extra",
        ]
        assert result.concat_job.name == "merge-trim-hp-concat"


class TestDagWiring:
    """The edges the stage inserts"""

    def test_shards_depend_on_upstream_and_concat_on_the_shards(
        self, tmp_path
    ):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        for job in result.shard_jobs:
            assert dag.waiting_jobs[job] == {upstream}
        assert dag.waiting_jobs[result.concat_job] == set(result.shard_jobs)

    def test_shards_are_roots_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        for job in result.shard_jobs:
            assert job in dag.ready_jobs

    def test_result_fields(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.jobs == [*result.shard_jobs, result.concat_job]
        assert result.terminal == {result.concat_job}


class TestFromPipeline:
    """`TransferConfig.from_pipeline` collecting a pipeline's attributes"""

    @staticmethod
    def make_pipeline(tmp_path: pathlib.Path, **kwargs) -> SimpleNamespace:
        """A stand-in carrying the four `TransferInputs` attributes"""
        config = make_config(tmp_path)
        defaults = dict(
            pop_vcf=config.pop_vcf,
            shards=config.shards,
            pop_vcf_contigs=config.pop_vcf_contigs,
            fai_data=config.fai_data,
        )
        defaults.update(kwargs)
        return SimpleNamespace(**defaults)

    def test_collects_the_four_values(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path)

        config = TransferConfig.from_pipeline(pipeline)

        assert config.pop_vcf == pipeline.pop_vcf
        assert config.shards == pipeline.shards
        assert config.pop_vcf_contigs == pipeline.pop_vcf_contigs
        assert config.fai_data == pipeline.fai_data

    def test_without_a_pop_vcf(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path, pop_vcf=None)

        with pytest.raises(ValueError, match="population VCF"):
            TransferConfig.from_pipeline(pipeline)
