"""
Unit tests for the long-read SV calling stage
"""

import pathlib

from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.sv import LongReadSVStage


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


def make_stage(tmp_path: pathlib.Path, **kwargs) -> LongReadSVStage:
    """A LongReadSVStage over one alignment"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.bam"],
        model=tmp_path / "bundle" / "longreadsv.model",
    )
    defaults.update(kwargs)
    return LongReadSVStage(**defaults)  # type: ignore[arg-type]


class TestCommand:
    """The `sentieon driver` command the stage builds"""

    def test_longreadsv_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert str(result.job.shell) == (
            f"sentieon driver --input {tmp_path}/sample.bam "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            f"--algo LongReadSV "
            f"--model {tmp_path}/bundle/longreadsv.model "
            f"{tmp_path}/output.sv.vcf.gz"
        )

    def test_interval(self, tmp_path):
        bed = tmp_path / "regions.bed"
        result = make_stage(tmp_path, interval=bed).add_to(DAG())

        assert f"--interval {bed}" in str(result.job.shell)

    def test_without_an_interval(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert "--interval" not in str(result.job.shell)

    def test_multiple_inputs(self, tmp_path):
        result = make_stage(
            tmp_path, inputs=[tmp_path / "one.bam", tmp_path / "two.bam"]
        ).add_to(DAG())

        shell = str(result.job.shell)
        assert f"--input {tmp_path}/one.bam" in shell
        assert f"--input {tmp_path}/two.bam" in shell

    def test_replace_rg_interleaves_with_the_inputs(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.bam", tmp_path / "two.bam"],
            replace_rg=[["rg1=@RG\\tID:rg1"], ["rg2=@RG\\tID:rg2"]],
        ).add_to(DAG())

        shell = str(result.job.shell)
        assert (
            "--replace_rg 'rg1=@RG\\tID:rg1' "
            f"--input {tmp_path}/one.bam "
            "--replace_rg 'rg2=@RG\\tID:rg2' "
            f"--input {tmp_path}/two.bam"
        ) in shell

    def test_without_replace_rg(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert "--replace_rg" not in str(result.job.shell)


class TestOutput:
    """The SV VCF is derived from the run's output VCF"""

    def test_sv_vcf(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.sv_vcf == tmp_path / "output.sv.vcf.gz"
        assert str(result.sv_vcf) in str(result.job.shell)


class TestJobMetadata:
    """Job name, task name and thread count"""

    def test_defaults(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.job.name == "LongReadSV"
        assert result.job.task_name == "sv-calling"
        assert result.job.threads == 4

    def test_custom_name_and_task(self, tmp_path):
        result = make_stage(
            tmp_path, name="longreadsv-lr", task_name="sv"
        ).add_to(DAG())

        assert result.job.name == "longreadsv-lr"
        assert result.job.task_name == "sv"


class TestDagWiring:
    """`build` touches no DAG; `add_to` inserts the job"""

    def test_build_inserts_nothing(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).build()

        assert result.jobs == [result.job]
        assert result.terminal == {result.job}
        assert not dag.ready_jobs
        assert not dag.waiting_jobs

    def test_upstream_lands_on_the_job(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.job] == {upstream}

    def test_a_root_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert result.job in dag.ready_jobs
