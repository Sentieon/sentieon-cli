"""
Alignment metrics collection
"""

from dataclasses import dataclass
from importlib.resources import files
import pathlib
import sys
from typing import Iterable, List, Optional

from ..dag import DAG
from ..driver import BaseAlgo
from ..job import Job
from ..shell_pipeline import Command, Pipeline
from .base import Stage, StageResult, driver_job


@dataclass(frozen=True)
class MetricsPaths:
    """The metrics directory and the files the pipelines write into it.

    Every metrics file is named after the output VCF, so the whole set is
    derived once, from :meth:`from_output_vcf`, and passed around instead
    of being rebuilt at each call site.
    """

    metrics_dir: pathlib.Path
    metric_base: str

    @classmethod
    def from_output_vcf(cls, output_vcf: pathlib.Path) -> "MetricsPaths":
        """The metrics paths belonging to an output VCF"""
        return cls(
            metrics_dir=pathlib.Path(
                str(output_vcf).replace(".vcf.gz", "_metrics")
            ),
            metric_base=output_vcf.name.replace(".vcf.gz", "") + ".txt",
        )

    def ensure_dir(self, dry_run: bool) -> None:
        """Create the metrics directory, unless this is a dry run"""
        if not dry_run:
            self.metrics_dir.mkdir(exist_ok=True)

    def _metric(self, suffix: str) -> pathlib.Path:
        return self.metrics_dir.joinpath(self.metric_base + suffix)

    @property
    def score(self) -> pathlib.Path:
        """The LocusCollector score file"""
        return self._metric(".score.txt.gz")

    @property
    def insert_size(self) -> pathlib.Path:
        """InsertSizeMetricAlgo output"""
        return self._metric(".insert_size.txt")

    @property
    def mean_qual_by_cycle(self) -> pathlib.Path:
        """MeanQualityByCycle output"""
        return self._metric(".mean_qual_by_cycle.txt")

    @property
    def base_distribution_by_cycle(self) -> pathlib.Path:
        """BaseDistributionByCycle output"""
        return self._metric(".base_distribution_by_cycle.txt")

    @property
    def qual_distribution(self) -> pathlib.Path:
        """QualDistribution output"""
        return self._metric(".qual_distribution.txt")

    @property
    def alignment_stat(self) -> pathlib.Path:
        """AlignmentStat output"""
        return self._metric(".alignment_stat.txt")

    @property
    def coverage(self) -> pathlib.Path:
        """CoverageMetrics output prefix (not named after the VCF)"""
        return self.metrics_dir.joinpath("coverage")

    @property
    def hybrid_selection(self) -> pathlib.Path:
        """HsMetricAlgo output, for WES runs"""
        return self._metric(".hybrid-selection.txt")

    @property
    def wgs(self) -> pathlib.Path:
        """WgsMetricsAlgo output, for WGS runs"""
        return self._metric(".wgs.txt")

    @property
    def gc_bias(self) -> pathlib.Path:
        """GCBias output"""
        return self._metric(".gc_bias.txt")

    @property
    def gc_bias_summary(self) -> pathlib.Path:
        """The GCBias summary output"""
        return self._metric(".gc_bias_summary.txt")

    @property
    def dedup_metrics(self) -> pathlib.Path:
        """The Dedup metrics output"""
        return self._metric(".dedup_metrics.txt")


@dataclass(kw_only=True)
class MetricsResult(StageResult):
    """The jobs of `MetricsStage`"""

    metrics_job: Job
    rehead_job: Optional[Job]


@dataclass(kw_only=True)
class MetricsStage(Stage):
    """Collect alignment metrics with a single `sentieon driver` pass.

    The caller supplies the ordered algo list, since which metrics are
    worth collecting depends on the assay and on whether the alignment
    has been deduplicated yet. When ``rehead_metrics`` names a WGS
    metrics file, a second job rewrites its header so MultiQC recognizes
    it.
    """

    inputs: List[pathlib.Path]
    algos: List[BaseAlgo]
    interval: Optional[pathlib.Path] = None
    rehead_metrics: Optional[pathlib.Path] = None
    name: str = "metrics"
    task_name: str = "metrics"
    threads: int = 0

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> MetricsResult:
        deps = set(upstream)

        metrics = driver_job(
            self.ctx,
            self.algos,
            inputs=self.inputs,
            interval=self.interval,
            name=self.name,
            task_name=self.task_name,
            threads=self.threads,
        )
        dag.add_job(metrics, deps)
        jobs = [metrics]

        rehead: Optional[Job] = None
        if self.rehead_metrics is not None:
            rehead_script = pathlib.Path(
                str(
                    files("sentieon_cli.scripts").joinpath(
                        "rehead_wgs_metrics.py"
                    )
                )
            )
            rehead = Job(
                Pipeline(
                    Command(
                        sys.executable,
                        str(rehead_script),
                        "--metrics_file",
                        str(self.rehead_metrics),
                    )
                ),
                "rehead-metrics",
                0,
                task_name="metrics",
            )
            dag.add_job(rehead, {metrics})
            jobs.append(rehead)

        return MetricsResult(
            jobs=jobs,
            terminal={rehead} if rehead else {metrics},
            metrics_job=metrics,
            rehead_job=rehead,
        )
