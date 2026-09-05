"""
Unit tests for the ploidy estimation stage
"""

import pathlib

from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.ploidy import PloidyStage


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


def make_stage(tmp_path: pathlib.Path, **kwargs) -> PloidyStage:
    """A PloidyStage over one alignment"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.cram"],
    )
    defaults.update(kwargs)
    return PloidyStage(**defaults)  # type: ignore[arg-type]


class TestOutput:
    """The ploidy JSON sits next to the output VCF"""

    def test_ploidy_json_derivation(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.ploidy_json == tmp_path / "output_ploidy.json"
        assert str(result.ploidy_json) in str(result.job.shell)


class TestCommand:
    """The `estimate_ploidy.py` command line"""

    def test_script_and_inputs(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.bam", tmp_path / "two.bam"],
        ).add_to(DAG())

        cmd = str(result.job.shell)
        assert cmd.startswith("python3 ")
        assert "estimate_ploidy.py" in cmd
        assert f"-i {tmp_path}/one.bam {tmp_path}/two.bam" in cmd
        assert f"--outfile {tmp_path}/output_ploidy.json" in cmd

    def test_the_script_defaults_apply_without_a_build(self, tmp_path):
        # The `estimate_ploidy.py` defaults are the chr-prefixed names
        cmd = str(make_stage(tmp_path).add_to(DAG()).job.shell)

        assert "--contigs" not in cmd
        assert "--autosomes" not in cmd
        assert "--x_contig" not in cmd
        assert "--y_contig" not in cmd

    def test_an_unrecognized_build_uses_the_defaults(self, tmp_path):
        cmd = str(
            make_stage(tmp_path, reference_build="hg38")
            .add_to(DAG())
            .job.shell
        )

        assert "--contigs" not in cmd
        assert "--x_contig" not in cmd

    def test_b37_contig_names(self, tmp_path):
        cmd = str(
            make_stage(tmp_path, reference_build="b37").add_to(DAG()).job.shell
        )

        assert "--contigs 1 2 " in cmd
        assert "--autosomes 1 2 " in cmd
        # The sex chromosomes are in `--contigs` but not in `--autosomes`
        assert (
            "--contigs "
            + " ".join([str(i) for i in range(1, 23)] + ["X", "Y"])
            in cmd
        )
        assert "--autosomes " + " ".join([str(i) for i in range(1, 23)]) in cmd
        assert "--x_contig X" in cmd
        assert "--y_contig Y" in cmd


class TestJobMetadata:
    """Job name, task name and thread count"""

    def test_metadata(self, tmp_path):
        job = make_stage(tmp_path).add_to(DAG()).job

        assert job.name == "estimate-ploidy"
        assert job.task_name == "ploidy"
        # The script is single-threaded and does not reserve a core
        assert job.threads == 0


class TestDagWiring:
    """The edges the stage inserts"""

    def test_upstream_edge(self, tmp_path):
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
        result = make_stage(tmp_path).add_to(DAG())

        assert result.jobs == [result.job]
        assert result.terminal == {result.job}
