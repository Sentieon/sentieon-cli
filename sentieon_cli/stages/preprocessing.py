"""
Short-read preprocessing: duplicate marking and metrics collection
"""

from dataclasses import dataclass
import pathlib
from typing import Iterable, List, Optional

from ..dag import DAG
from ..driver import (
    AlignmentStat,
    BaseAlgo,
    BaseDistributionByCycle,
    CoverageMetrics,
    GCBias,
    HsMetricAlgo,
    InsertSizeMetricAlgo,
    MeanQualityByCycle,
    QualDistribution,
    WgsMetricsAlgo,
)
from ..job import Job
from .base import Stage, StageResult
from .dedup import DedupStage
from .metrics import MetricsPaths, MetricsStage


def pre_dedup_metrics_algos(paths: MetricsPaths, assay: str) -> List[BaseAlgo]:
    """The metrics that do not need duplicate-marked reads"""
    algos: List[BaseAlgo] = [
        MeanQualityByCycle(paths.mean_qual_by_cycle),
        BaseDistributionByCycle(paths.base_distribution_by_cycle),
        QualDistribution(paths.qual_distribution),
        AlignmentStat(paths.alignment_stat),
    ]
    if assay == "WGS":
        algos.append(GCBias(paths.gc_bias, summary=paths.gc_bias_summary))
    return algos


@dataclass(kw_only=True)
class SrPreprocessingResult(StageResult):
    """The short-read dedup and metrics jobs added to a DAG.

    * ``deduped`` -- the alignment to call variants from; the input files
      unchanged when duplicate marking is off.
    * ``dedup_job`` -- the Dedup job, absent when duplicate marking is off.
    * ``qc_jobs`` -- every metrics-producing job, for MultiQC to wait on.
    """

    deduped: List[pathlib.Path]
    dedup_job: Optional[Job]
    qc_jobs: List[Job]


@dataclass(kw_only=True)
class ShortReadPreprocessingStage(Stage):
    """Mark duplicates and collect metrics from the short reads.

    Which metrics are collected, and whether they ride along on the
    LocusCollector pass or run after duplicate marking, depends on the
    assay and on whether a BED file restricts the run.

    ``terminal`` is the Dedup job, so downstream stages read the
    deduplicated alignment. With ``duplicate_marking="none"`` there is no
    Dedup job and the stage adds only a leaf metrics job that nothing
    should wait on, so ``terminal`` is empty and callers depend on the
    jobs that produced ``inputs`` instead.
    """

    inputs: List[pathlib.Path]
    duplicate_marking: str = "markdup"
    consensus: bool = False
    skip_metrics: bool = False
    assay: str = "WGS"
    bed: Optional[pathlib.Path] = None
    bam_format: bool = False

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> SrPreprocessingResult:
        ctx = self.ctx
        suffix = "bam" if self.bam_format else "cram"

        paths = MetricsPaths.from_output_vcf(ctx.output_vcf)
        paths.ensure_dir(ctx.dry_run)

        if self.duplicate_marking == "none":
            # Without duplicate marking there is nothing to collect after
            # the fact, so every metric comes from one pass over the input
            if self.skip_metrics:
                return SrPreprocessingResult(
                    jobs=[],
                    terminal=set(),
                    deduped=self.inputs,
                    dedup_job=None,
                    qc_jobs=[],
                )
            metrics_result = MetricsStage(
                ctx=ctx,
                inputs=self.inputs,
                algos=[
                    InsertSizeMetricAlgo(paths.insert_size),
                    *pre_dedup_metrics_algos(paths, self.assay),
                ],
                name="metrics",
                task_name="metrics",
                threads=ctx.cores,
            ).add_to(dag, upstream)
            return SrPreprocessingResult(
                jobs=list(metrics_result.jobs),
                terminal=set(),
                deduped=self.inputs,
                dedup_job=None,
                qc_jobs=list(metrics_result.jobs),
            )

        # Metrics that ride along on the LocusCollector pass
        lc_extra_algos: List[BaseAlgo] = []
        if not self.skip_metrics:
            # Prefer to run InsertSizeMetricAlgo after duplicate marking
            if self.assay == "WES" and not self.bed:
                lc_extra_algos.append(InsertSizeMetricAlgo(paths.insert_size))
            lc_extra_algos.extend(pre_dedup_metrics_algos(paths, self.assay))

        deduped = pathlib.Path(
            str(ctx.output_vcf).replace(".vcf.gz", f"_deduped.{suffix}")
        )
        dedup_result = DedupStage(
            ctx=ctx,
            inputs=self.inputs,
            output=deduped,
            score_file=paths.score,
            consensus=self.consensus,
            rmdup=(self.duplicate_marking == "rmdup"),
            cram_write_options="version=3.0,compressor=rans",
            dedup_metrics=paths.dedup_metrics,
            lc_extra_algos=lc_extra_algos,
        ).add_to(dag, upstream)
        jobs: List[Job] = list(dedup_result.jobs)
        qc_jobs: List[Job] = [dedup_result.lc_job]

        # HsMetricAlgo and WgsMetricsAlgo run after duplicate marking to
        # account for duplicate reads
        post_algos: List[BaseAlgo] = []
        rehead_metrics: Optional[pathlib.Path] = None
        if not self.skip_metrics:
            if self.assay == "WES" and self.bed:
                post_algos = [
                    HsMetricAlgo(paths.hybrid_selection, self.bed, self.bed),
                    InsertSizeMetricAlgo(paths.insert_size),
                ]
            elif self.assay == "WGS":
                post_algos = [
                    InsertSizeMetricAlgo(paths.insert_size),
                    WgsMetricsAlgo(paths.wgs, include_unpaired="true"),
                    CoverageMetrics(paths.coverage),
                ]
                # Rehead WGS metrics so they are recognized by MultiQC
                rehead_metrics = paths.wgs
        if post_algos:
            metrics_result = MetricsStage(
                ctx=ctx,
                inputs=[deduped],
                algos=post_algos,
                interval=self.bed,
                rehead_metrics=rehead_metrics,
                threads=0,  # Run metrics in the background
            ).add_to(dag, {dedup_result.dedup_job})
            jobs.extend(metrics_result.jobs)
            qc_jobs.extend(metrics_result.jobs)

        return SrPreprocessingResult(
            jobs=jobs,
            terminal={dedup_result.dedup_job},
            deduped=[deduped],
            dedup_job=dedup_result.dedup_job,
            qc_jobs=qc_jobs,
        )
