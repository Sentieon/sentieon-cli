"""
Unit tests for the metrics stage and the metrics paths
"""

import pathlib
import sys

from sentieon_cli.dag import DAG
from sentieon_cli.driver import (
    AlignmentStat,
    InsertSizeMetricAlgo,
    WgsMetricsAlgo,
)
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.metrics import MetricsPaths, MetricsStage


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


def make_stage(tmp_path: pathlib.Path, **kwargs) -> MetricsStage:
    """A MetricsStage with the paths a caller supplies"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "deduped.cram"],
        algos=[AlignmentStat(tmp_path / "as.txt")],
    )
    defaults.update(kwargs)
    return MetricsStage(**defaults)  # type: ignore[arg-type]


class TestMetricsPaths:
    """Every metrics file is derived from the output VCF"""

    def test_directory_and_base(self):
        paths = MetricsPaths.from_output_vcf(
            pathlib.Path("/data/run/sample.vcf.gz")
        )

        assert paths.metrics_dir == pathlib.Path("/data/run/sample_metrics")
        assert paths.metric_base == "sample.txt"

    def test_metric_files(self):
        paths = MetricsPaths.from_output_vcf(pathlib.Path("/data/out.vcf.gz"))
        base = pathlib.Path("/data/out_metrics")

        assert paths.score == base / "out.txt.score.txt.gz"
        assert paths.insert_size == base / "out.txt.insert_size.txt"
        assert (
            paths.mean_qual_by_cycle == base / "out.txt.mean_qual_by_cycle.txt"
        )
        assert paths.base_distribution_by_cycle == (
            base / "out.txt.base_distribution_by_cycle.txt"
        )
        assert paths.qual_distribution == (
            base / "out.txt.qual_distribution.txt"
        )
        assert paths.alignment_stat == base / "out.txt.alignment_stat.txt"
        assert paths.hybrid_selection == base / "out.txt.hybrid-selection.txt"
        assert paths.wgs == base / "out.txt.wgs.txt"
        assert paths.gc_bias == base / "out.txt.gc_bias.txt"
        assert paths.gc_bias_summary == base / "out.txt.gc_bias_summary.txt"
        assert paths.dedup_metrics == base / "out.txt.dedup_metrics.txt"

    def test_coverage_is_not_named_after_the_vcf(self):
        paths = MetricsPaths.from_output_vcf(pathlib.Path("/data/out.vcf.gz"))

        assert paths.coverage == pathlib.Path("/data/out_metrics/coverage")

    def test_ensure_dir_creates_the_directory(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        paths.ensure_dir(dry_run=False)

        assert paths.metrics_dir.is_dir()

    def test_ensure_dir_creates_nothing_on_a_dry_run(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        paths.ensure_dir(dry_run=True)

        assert not paths.metrics_dir.exists()

    def test_ensure_dir_tolerates_an_existing_directory(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        paths.metrics_dir.mkdir()
        paths.ensure_dir(dry_run=False)

        assert paths.metrics_dir.is_dir()


class TestCommands:
    """The commands the stage builds"""

    def test_metrics_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert str(result.metrics_job.shell) == (
            f"sentieon driver --input {tmp_path}/deduped.cram "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            f"--algo AlignmentStat --adapter_seq '' {tmp_path}/as.txt"
        )

    def test_algo_order_is_the_callers(self, tmp_path):
        result = make_stage(
            tmp_path,
            algos=[
                InsertSizeMetricAlgo(tmp_path / "is.txt"),
                WgsMetricsAlgo(tmp_path / "wgs.txt", include_unpaired="true"),
                AlignmentStat(tmp_path / "as.txt"),
            ],
        ).add_to(DAG())

        shell = str(result.metrics_job.shell)
        assert shell.index("--algo InsertSizeMetricAlgo") < shell.index(
            "--algo WgsMetricsAlgo"
        )
        assert shell.index("--algo WgsMetricsAlgo") < shell.index(
            "--algo AlignmentStat"
        )

    def test_interval(self, tmp_path):
        bed = tmp_path / "regions.bed"
        result = make_stage(tmp_path, interval=bed).add_to(DAG())

        assert f"--interval {bed}" in str(result.metrics_job.shell)

    def test_without_an_interval(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert "--interval" not in str(result.metrics_job.shell)

    def test_multiple_inputs(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.cram", tmp_path / "two.cram"],
        ).add_to(DAG())

        shell = str(result.metrics_job.shell)
        assert f"--input {tmp_path}/one.cram" in shell
        assert f"--input {tmp_path}/two.cram" in shell


class TestRehead:
    """The optional WGS-metrics rehead job"""

    def test_rehead_command(self, tmp_path):
        wgs = tmp_path / "output.txt.wgs.txt"
        result = make_stage(tmp_path, rehead_metrics=wgs).add_to(DAG())

        assert result.rehead_job is not None
        shell = str(result.rehead_job.shell)
        assert shell.startswith(sys.executable)
        assert "rehead_wgs_metrics.py" in shell
        assert f"--metrics_file {wgs}" in shell

    def test_rehead_job_metadata(self, tmp_path):
        result = make_stage(
            tmp_path, rehead_metrics=tmp_path / "wgs.txt"
        ).add_to(DAG())

        assert result.rehead_job is not None
        assert result.rehead_job.name == "rehead-metrics"
        assert result.rehead_job.threads == 0
        assert result.rehead_job.task_name == "metrics"

    def test_no_rehead_job_by_default(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.rehead_job is None
        assert result.jobs == [result.metrics_job]
        assert result.terminal == {result.metrics_job}

    def test_terminal_is_the_rehead_job(self, tmp_path):
        result = make_stage(
            tmp_path, rehead_metrics=tmp_path / "wgs.txt"
        ).add_to(DAG())

        assert result.jobs == [result.metrics_job, result.rehead_job]
        assert result.terminal == {result.rehead_job}


class TestJobMetadata:
    """Job name, task name and thread count"""

    def test_defaults(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.metrics_job.name == "metrics"
        assert result.metrics_job.task_name == "metrics"
        # Metrics run in the background by default
        assert result.metrics_job.threads == 0

    def test_custom_name_task_and_threads(self, tmp_path):
        result = make_stage(
            tmp_path,
            name="sr-metrics",
            task_name="qc",
            threads=4,
        ).add_to(DAG())

        assert result.metrics_job.name == "sr-metrics"
        assert result.metrics_job.task_name == "qc"
        assert result.metrics_job.threads == 4

    def test_rehead_keeps_the_metrics_task(self, tmp_path):
        result = make_stage(
            tmp_path,
            name="sr-metrics",
            task_name="qc",
            rehead_metrics=tmp_path / "wgs.txt",
        ).add_to(DAG())

        assert result.rehead_job is not None
        assert result.rehead_job.task_name == "metrics"


class TestDagWiring:
    """The edges the stage inserts"""

    def test_upstream_lands_on_the_metrics_job(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(
            tmp_path, rehead_metrics=tmp_path / "wgs.txt"
        ).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.metrics_job] == {upstream}
        assert dag.waiting_jobs[result.rehead_job] == {result.metrics_job}

    def test_metrics_is_a_root_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert result.metrics_job in dag.ready_jobs
