"""
Unit tests for the short-read preprocessing stage
"""

import pathlib

from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.metrics import MetricsPaths
from sentieon_cli.stages.preprocessing import (
    ShortReadPreprocessingStage,
    pre_dedup_metrics_algos,
)

# The pre-dedup metrics, in the order the stage collects them
PRE_DEDUP_ALGOS = [
    "MeanQualityByCycle",
    "BaseDistributionByCycle",
    "QualDistribution",
    "AlignmentStat",
]


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


def make_stage(
    tmp_path: pathlib.Path, **kwargs
) -> ShortReadPreprocessingStage:
    """A preprocessing stage over one input alignment"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.bam"],
    )
    defaults.update(kwargs)
    return ShortReadPreprocessingStage(**defaults)  # type: ignore[arg-type]


def algo_order(shell: str) -> list:
    """The algos a driver command runs, in command-line order"""
    tokens = shell.split()
    return [
        tokens[i + 1]
        for i, token in enumerate(tokens[:-1])
        if token == "--algo"
    ]


class TestPreDedupAlgos:
    """The metrics that do not need duplicate-marked reads"""

    def test_wgs_adds_gc_bias(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        algos = pre_dedup_metrics_algos(paths, "WGS")

        assert [algo.name for algo in algos] == PRE_DEDUP_ALGOS + ["GCBias"]

    def test_wes_omits_gc_bias(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        algos = pre_dedup_metrics_algos(paths, "WES")

        assert [algo.name for algo in algos] == PRE_DEDUP_ALGOS

    def test_the_outputs_come_from_the_metrics_paths(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        algos = pre_dedup_metrics_algos(paths, "WGS")

        assert algos[0].output == paths.mean_qual_by_cycle
        assert algos[1].output == paths.base_distribution_by_cycle
        assert algos[2].output == paths.qual_distribution
        assert algos[3].output == paths.alignment_stat
        assert algos[4].output == paths.gc_bias
        assert algos[4].summary == paths.gc_bias_summary


class TestNoDuplicateMarking:
    """`duplicate_marking="none"` collects everything in one pass"""

    def test_skip_metrics_adds_no_jobs(self, tmp_path):
        dag = DAG()
        inputs = [tmp_path / "sample.bam"]
        result = make_stage(
            tmp_path,
            inputs=inputs,
            duplicate_marking="none",
            skip_metrics=True,
        ).add_to(dag)

        assert result.jobs == []
        assert result.terminal == set()
        assert result.qc_jobs == []
        assert result.dedup_job is None
        assert result.deduped == inputs
        assert not dag.ready_jobs
        assert not dag.waiting_jobs

    def test_one_metrics_job(self, tmp_path):
        inputs = [tmp_path / "sample.bam"]
        result = make_stage(
            tmp_path, inputs=inputs, duplicate_marking="none"
        ).add_to(DAG())

        (metrics,) = result.jobs
        assert metrics.name == "metrics"
        assert metrics.task_name == "metrics"
        # This pass reads the whole alignment, so it takes every core
        assert metrics.threads == 4
        assert result.dedup_job is None
        assert result.deduped == inputs
        assert result.qc_jobs == [metrics]

    def test_wgs_algo_order(self, tmp_path):
        result = make_stage(tmp_path, duplicate_marking="none").add_to(DAG())

        assert algo_order(str(result.jobs[0].shell)) == (
            ["InsertSizeMetricAlgo"] + PRE_DEDUP_ALGOS + ["GCBias"]
        )

    def test_wes_omits_gc_bias(self, tmp_path):
        result = make_stage(
            tmp_path, duplicate_marking="none", assay="WES"
        ).add_to(DAG())

        assert algo_order(str(result.jobs[0].shell)) == (
            ["InsertSizeMetricAlgo"] + PRE_DEDUP_ALGOS
        )

    def test_the_metrics_job_reads_the_inputs(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.bam", tmp_path / "two.bam"],
            duplicate_marking="none",
        ).add_to(DAG())

        shell = str(result.jobs[0].shell)
        assert f"--input {tmp_path}/one.bam" in shell
        assert f"--input {tmp_path}/two.bam" in shell

    def test_upstream_lands_on_the_metrics_job(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path, duplicate_marking="none").add_to(
            dag, [upstream]
        )

        assert dag.waiting_jobs[result.jobs[0]] == {upstream}
        # Nothing should wait on this leaf
        assert result.terminal == set()


class TestMarkdupWgs:
    """The default: duplicate marking plus WGS metrics"""

    def test_result_fields(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        lc, dedup, metrics, rehead = result.jobs
        assert lc.name == "locuscollector"
        assert dedup.name == "dedup"
        assert metrics.name == "metrics"
        assert rehead.name == "rehead-metrics"
        assert result.dedup_job is dedup
        assert result.terminal == {dedup}
        assert result.qc_jobs == [lc, metrics, rehead]
        assert result.deduped == [tmp_path / "output_deduped.cram"]

    def test_locuscollector_carries_the_pre_dedup_metrics(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert algo_order(str(result.jobs[0].shell)) == (
            ["LocusCollector"] + PRE_DEDUP_ALGOS + ["GCBias"]
        )

    def test_dedup_writes_metrics_and_cram_options(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        result = make_stage(tmp_path).add_to(DAG())

        shell = str(result.dedup_job.shell)
        assert f"--metrics {paths.dedup_metrics}" in shell
        assert "--cram_write_options version=3.0,compressor=rans" in shell
        assert "--rmdup" not in shell

    def test_post_dedup_metrics(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        metrics = result.jobs[2]
        assert algo_order(str(metrics.shell)) == [
            "InsertSizeMetricAlgo",
            "WgsMetricsAlgo",
            "CoverageMetrics",
        ]
        assert "--include_unpaired true" in str(metrics.shell)
        # Metrics run in the background
        assert metrics.threads == 0
        assert f"--input {tmp_path}/output_deduped.cram" in str(metrics.shell)

    def test_the_post_job_waits_on_dedup(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert dag.waiting_jobs[result.jobs[2]] == {result.dedup_job}
        assert dag.waiting_jobs[result.jobs[3]] == {result.jobs[2]}

    def test_upstream_lands_on_locuscollector(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.jobs[0]] == {upstream}

    def test_insert_size_is_not_on_the_locuscollector_pass(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert "InsertSizeMetricAlgo" not in algo_order(
            str(result.jobs[0].shell)
        )


class TestWes:
    """WES runs collect hybrid-selection metrics instead"""

    def test_with_a_bed(self, tmp_path):
        bed = tmp_path / "regions.bed"
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        result = make_stage(tmp_path, assay="WES", bed=bed).add_to(DAG())

        lc, dedup, metrics = result.jobs
        assert algo_order(str(metrics.shell)) == [
            "HsMetricAlgo",
            "InsertSizeMetricAlgo",
        ]
        assert f"--interval {bed}" in str(metrics.shell)
        assert f"--targets_list {bed}" in str(metrics.shell)
        assert str(paths.hybrid_selection) in str(metrics.shell)
        # No WGS metrics means no rehead job
        assert result.qc_jobs == [lc, metrics]
        assert result.terminal == {dedup}

    def test_insert_size_is_not_on_the_locuscollector_pass_with_a_bed(
        self, tmp_path
    ):
        result = make_stage(
            tmp_path, assay="WES", bed=tmp_path / "regions.bed"
        ).add_to(DAG())

        assert algo_order(str(result.jobs[0].shell)) == (
            ["LocusCollector"] + PRE_DEDUP_ALGOS
        )

    def test_without_a_bed(self, tmp_path):
        result = make_stage(tmp_path, assay="WES").add_to(DAG())

        lc, dedup = result.jobs
        # InsertSize moves onto the LocusCollector pass
        assert algo_order(str(lc.shell)) == (
            ["LocusCollector", "InsertSizeMetricAlgo"] + PRE_DEDUP_ALGOS
        )
        assert result.qc_jobs == [lc]
        assert result.terminal == {dedup}


class TestOptions:
    """The flags a caller passes through"""

    def test_rmdup(self, tmp_path):
        result = make_stage(tmp_path, duplicate_marking="rmdup").add_to(DAG())

        assert "--rmdup" in str(result.dedup_job.shell)

    def test_consensus(self, tmp_path):
        result = make_stage(tmp_path, consensus=True).add_to(DAG())

        assert "--algo LocusCollector --consensus" in str(result.jobs[0].shell)

    def test_bam_format(self, tmp_path):
        result = make_stage(tmp_path, bam_format=True).add_to(DAG())

        assert result.deduped == [tmp_path / "output_deduped.bam"]

    def test_skip_metrics_keeps_duplicate_marking(self, tmp_path):
        paths = MetricsPaths.from_output_vcf(tmp_path / "output.vcf.gz")
        result = make_stage(tmp_path, skip_metrics=True).add_to(DAG())

        lc, dedup = result.jobs
        assert algo_order(str(lc.shell)) == ["LocusCollector"]
        # Dedup still writes its own metrics
        assert f"--metrics {paths.dedup_metrics}" in str(dedup.shell)
        assert result.qc_jobs == [lc]
        assert result.terminal == {dedup}


class TestMetricsDirectory:
    """The metrics directory is created before the jobs are built"""

    def test_created_when_not_a_dry_run(self, tmp_path):
        ctx = StageContext(
            reference=tmp_path / "ref.fa",
            output_vcf=tmp_path / "output.vcf.gz",
            tmp_dir=tmp_path,
            cores=4,
            dry_run=False,
            skip_version_check=True,
        )
        make_stage(tmp_path, ctx=ctx).add_to(DAG())

        assert (tmp_path / "output_metrics").is_dir()

    def test_not_created_on_a_dry_run(self, tmp_path):
        make_stage(tmp_path).add_to(DAG())

        assert not (tmp_path / "output_metrics").exists()

    def test_created_even_when_no_jobs_are_added(self, tmp_path):
        ctx = StageContext(
            reference=tmp_path / "ref.fa",
            output_vcf=tmp_path / "output.vcf.gz",
            tmp_dir=tmp_path,
            cores=4,
            dry_run=False,
            skip_version_check=True,
        )
        make_stage(
            tmp_path,
            ctx=ctx,
            duplicate_marking="none",
            skip_metrics=True,
        ).add_to(DAG())

        assert (tmp_path / "output_metrics").is_dir()
