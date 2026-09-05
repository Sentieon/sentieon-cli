"""
Unit tests for the duplicate-marking stage
"""

import pathlib

from sentieon_cli.dag import DAG
from sentieon_cli.driver import AlignmentStat, GCBias
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.dedup import DedupStage


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


def make_stage(tmp_path: pathlib.Path, **kwargs) -> DedupStage:
    """A DedupStage with the paths a caller supplies"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.bam"],
        output=tmp_path / "deduped.cram",
        score_file=tmp_path / "score.txt.gz",
    )
    defaults.update(kwargs)
    return DedupStage(**defaults)  # type: ignore[arg-type]


class TestCommands:
    """The `sentieon driver` commands the stage builds"""

    def test_locuscollector_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert str(result.lc_job.shell) == (
            f"sentieon driver --input {tmp_path}/sample.bam "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            f"--algo LocusCollector {tmp_path}/score.txt.gz"
        )

    def test_dedup_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        # The defaults produce the shortest possible Dedup command
        assert str(result.dedup_job.shell) == (
            f"sentieon driver --input {tmp_path}/sample.bam "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            f"--algo Dedup --score_info {tmp_path}/score.txt.gz "
            f"{tmp_path}/deduped.cram"
        )

    def test_consensus_reads(self, tmp_path):
        result = make_stage(tmp_path, consensus=True).add_to(DAG())

        assert "--algo LocusCollector --consensus" in str(result.lc_job.shell)
        assert "--consensus" not in str(result.dedup_job.shell)

    def test_rmdup_and_cram_write_options(self, tmp_path):
        result = make_stage(
            tmp_path,
            rmdup=True,
            cram_write_options="version=3.0,compressor=rans",
        ).add_to(DAG())

        shell = str(result.dedup_job.shell)
        assert "--cram_write_options version=3.0,compressor=rans" in shell
        assert "--rmdup" in shell
        # The options belong to Dedup alone
        assert "--rmdup" not in str(result.lc_job.shell)

    def test_dedup_metrics(self, tmp_path):
        metrics = tmp_path / "sample.dedup_metrics.txt"
        result = make_stage(tmp_path, dedup_metrics=metrics).add_to(DAG())

        assert f"--metrics {metrics}" in str(result.dedup_job.shell)

    def test_without_optional_dedup_arguments(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        shell = str(result.dedup_job.shell)
        for flag in ("--cram_write_options", "--metrics", "--rmdup"):
            assert flag not in shell

    def test_read_filters_reach_both_passes(self, tmp_path):
        read_filter = "IndelLeftAlignReadTransform,rgid=rg1-mm2"
        result = make_stage(tmp_path, read_filters=[read_filter]).add_to(DAG())

        for job in (result.lc_job, result.dedup_job):
            assert f"--read_filter {read_filter}" in str(job.shell)

    def test_no_read_filters_emits_nothing(self, tmp_path):
        result = make_stage(tmp_path, read_filters=[]).add_to(DAG())

        for job in (result.lc_job, result.dedup_job):
            assert "--read_filter" not in str(job.shell)

    def test_multiple_inputs(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.bam", tmp_path / "two.bam"],
        ).add_to(DAG())

        for job in (result.lc_job, result.dedup_job):
            shell = str(job.shell)
            assert f"--input {tmp_path}/one.bam" in shell
            assert f"--input {tmp_path}/two.bam" in shell


class TestExtraAlgos:
    """Metrics that ride along on the LocusCollector pass"""

    def test_extra_algos_follow_locuscollector(self, tmp_path):
        result = make_stage(
            tmp_path,
            lc_extra_algos=[
                AlignmentStat(tmp_path / "as.txt"),
                GCBias(tmp_path / "gc.txt", summary=tmp_path / "sum.txt"),
            ],
        ).add_to(DAG())

        shell = str(result.lc_job.shell)
        assert shell.index("--algo LocusCollector") < shell.index(
            "--algo AlignmentStat"
        )
        assert shell.index("--algo AlignmentStat") < shell.index(
            "--algo GCBias"
        )
        # They stay off the Dedup pass
        assert "--algo AlignmentStat" not in str(result.dedup_job.shell)


class TestJobMetadata:
    """Job names, task name and thread counts"""

    def test_default_names(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.lc_job.name == "locuscollector"
        assert result.dedup_job.name == "dedup"

    def test_tagged_names(self, tmp_path):
        result = make_stage(tmp_path, tag="bwa").add_to(DAG())

        assert result.lc_job.name == "locuscollector-bwa"
        assert result.dedup_job.name == "dedup-bwa"

    def test_threads_default_to_the_core_count(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.jobs:
            assert job.threads == 4
            assert job.task_name == "dedup"

    def test_overridden_threads_and_task_name(self, tmp_path):
        result = make_stage(
            tmp_path, threads=1, task_name="preprocess"
        ).add_to(DAG())

        for job in result.jobs:
            assert job.threads == 1
            assert job.task_name == "preprocess"
        # The driver still uses every core
        assert "--thread_count 4" in str(result.lc_job.shell)


class TestDagWiring:
    """The edges the stage inserts"""

    def test_upstream_lands_on_locuscollector_only(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.lc_job] == {upstream}
        assert dag.waiting_jobs[result.dedup_job] == {result.lc_job}

    def test_locuscollector_is_a_root_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert result.lc_job in dag.ready_jobs
        assert dag.waiting_jobs[result.dedup_job] == {result.lc_job}

    def test_result_fields(self, tmp_path):
        output = tmp_path / "deduped.cram"
        result = make_stage(tmp_path, output=output).add_to(DAG())

        assert result.jobs == [result.lc_job, result.dedup_job]
        assert result.terminal == {result.dedup_job}
        assert result.output == output
