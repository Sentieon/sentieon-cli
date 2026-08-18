"""
A base class for pangenome pipelines
"""

import copy
from enum import Enum
import json
import pathlib
import sys
from typing import List, Optional, Tuple

from importlib.resources import files

from . import command_strings as cmds
from .driver import (
    AlignmentStat,
    BaseDistributionByCycle,
    CoverageMetrics,
    Dedup,
    Driver,
    GCBias,
    InsertSizeMetricAlgo,
    LocusCollector,
    MeanQualityByCycle,
    QualDistribution,
    WgsMetricsAlgo,
)
from .job import Job
from .pipeline import BasePipeline
from .shell_pipeline import Command, Pipeline
from .util import path_arg


class SampleSex(Enum):
    FEMALE = 1
    MALE = 2
    UNKNOWN = 3


class BasePangenome(BasePipeline):
    """A pipeline base class for short reads"""

    params = copy.deepcopy(BasePipeline.params)
    params.update(
        {
            # Required arguments
            "gbz": {
                "help": "The pangenome graph file in GBZ format.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "hapl": {
                "help": "The haplotype file.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "model_bundle": {
                "flags": ["-m", "--model_bundle"],
                "help": "The model bundle file.",
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "r1_fastq": {
                "nargs": "*",
                "help": "Sample R1 fastq files.",
                "type": path_arg(exists=True, is_file=True),
            },
            "r2_fastq": {
                "nargs": "*",
                "help": "Sample R2 fastq files.",
                "type": path_arg(exists=True, is_file=True),
            },
            # Additional arguments
            "bam_format": {
                "help": (
                    "Use the BAM format instead of CRAM for output aligned "
                    "files."
                ),
                "action": "store_true",
            },
            "dbsnp": {
                "flags": ["-d", "--dbsnp"],
                "help": (
                    "dbSNP vcf file Supplying this file will annotate "
                    "variants with their dbSNP refSNP ID numbers."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "kmer_memory": {
                "help": "Memory limit for KMC in GB.",
                "default": 30,
                "type": int,
            },
            "pcr_free": {
                "help": "Use arguments for PCR-free data processing",
                "action": "store_true",
            },
        }
    )

    positionals = BasePipeline.positionals

    def __init__(self) -> None:
        super().__init__()
        self.gbz: Optional[pathlib.Path] = None
        self.hapl: Optional[pathlib.Path] = None
        self.model_bundle: Optional[pathlib.Path] = None
        self.r1_fastq: List[pathlib.Path] = []
        self.r2_fastq: List[pathlib.Path] = []
        self.bam_format = False
        self.dbsnp: Optional[pathlib.Path] = None
        self.kmer_memory = 30
        self.pcr_free = False

    def build_kmc_job(
        self, kmer_prefix: pathlib.Path, job_threads: int
    ) -> Job:
        """Build KMC k-mer counting jobs"""
        # Create file list for KMC
        file_list = pathlib.Path(str(kmer_prefix) + ".paths")
        all_fastqs = []

        # Add R1 files
        all_fastqs.extend(self.r1_fastq)

        # Add R2 files if present
        if self.r2_fastq:
            all_fastqs.extend(self.r2_fastq)

        # Write file list
        if not self.dry_run:
            with open(file_list, "w") as f:
                for fq in all_fastqs:
                    f.write(f"{fq}\n")

        # Create KMC job
        kmc_job = Job(
            cmds.cmd_kmc(
                kmer_prefix,
                file_list,
                self.tmp_dir,
                memory=self.kmer_memory,
                threads=self.cores,
            ),
            "kmc",
            job_threads,
            task_name="kmer-counting",
        )

        return kmc_job

    def build_ploidy_job(
        self,
        ploidy_json: pathlib.Path,
        deduped_bam: List[pathlib.Path],
    ) -> Job:
        """Estimate sample ploidy and sex"""
        estimate_ploidy = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("estimate_ploidy.py"))
        ).resolve()
        ploidy_job = Job(
            cmds.cmd_estimate_ploidy(
                ploidy_json,
                deduped_bam,
                estimate_ploidy,
            ),
            "estimate-ploidy",
            0,
            task_name="ploidy",
        )
        return ploidy_job

    def build_dedup_job(
        self,
        output_bam,
        input_bam: List[pathlib.Path],
        tag: str,
        left_align_rgid: Optional[str] = None,
        metrics: Optional[pathlib.Path] = None,
    ) -> Tuple[Job, Job]:
        """Build deduplication job"""
        score_file = self.tmp_dir.joinpath(f"sample-{tag}-score.txt.gz")

        read_filters = []
        if left_align_rgid:
            read_filters.append(
                f"IndelLeftAlignReadTransform,rgid={left_align_rgid}"
            )

        # LocusCollector + Dedup
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bam,
            read_filter=read_filters,
        )
        driver.add_algo(LocusCollector(score_file))

        lc_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            f"locuscollector-{tag}",
            self.cores,
            task_name="dedup",
        )

        driver2 = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bam,
            read_filter=read_filters,
        )
        driver2.add_algo(Dedup(output_bam, score_file, metrics=metrics))

        dedup_job = Job(
            Pipeline(Command(*driver2.build_cmd())),
            f"dedup-{tag}",
            self.cores,
            task_name="dedup",
        )

        return lc_job, dedup_job

    def build_metrics_job(
        self,
        sample_input: List[pathlib.Path],
    ) -> Tuple[Job, Job]:
        """Build a metrics job"""
        if not self.output_vcf:
            self.logger.error("output_vcf is required")
            sys.exit(2)

        # Create the metrics directory
        sample_name = self.output_vcf.name.replace(".vcf.gz", "")
        metric_base = sample_name + ".txt"
        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )
        if not self.dry_run:
            metrics_dir.mkdir(exist_ok=True)

        is_metrics = metrics_dir.joinpath(metric_base + ".insert_size.txt")
        mqbc_metrics = metrics_dir.joinpath(
            metric_base + ".mean_qual_by_cycle.txt"
        )
        bdbc_metrics = metrics_dir.joinpath(
            metric_base + ".base_distribution_by_cycle.txt"
        )
        qualdist_metrics = metrics_dir.joinpath(
            metric_base + ".qual_distribution.txt"
        )
        as_metrics = metrics_dir.joinpath(metric_base + ".alignment_stat.txt")
        coverage_metrics = metrics_dir.joinpath("coverage")

        # WGS metrics
        wgs_metrics = metrics_dir.joinpath(metric_base + ".wgs.txt")
        gc_metrics = metrics_dir.joinpath(metric_base + ".gc_bias.txt")
        gc_summary = metrics_dir.joinpath(metric_base + ".gc_bias_summary.txt")

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=sample_input,
        )

        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        driver.add_algo(MeanQualityByCycle(mqbc_metrics))
        driver.add_algo(BaseDistributionByCycle(bdbc_metrics))
        driver.add_algo(QualDistribution(qualdist_metrics))
        driver.add_algo(AlignmentStat(as_metrics))
        driver.add_algo(GCBias(gc_metrics, summary=gc_summary))
        driver.add_algo(WgsMetricsAlgo(wgs_metrics, include_unpaired="true"))
        driver.add_algo(CoverageMetrics(coverage_metrics))

        metrics_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "metrics",
            0,
            task_name="metrics",
        )

        rehead_script = pathlib.Path(
            str(
                files("sentieon_cli.scripts").joinpath("rehead_wgs_metrics.py")
            )
        )
        rehead_job = Job(
            Pipeline(
                Command(
                    sys.executable,
                    str(rehead_script),
                    "--metrics_file",
                    str(wgs_metrics),
                )
            ),
            "Rehead metrics",
            0,
            task_name="metrics",
        )
        return (metrics_job, rehead_job)

    def get_sex(self, ploidy_json: pathlib.Path) -> None:
        """Retrieve the sample sex"""
        if self.dry_run:
            self.logger.info("Setting sample sex to MALE for dry-run")
            self.sample_sex = SampleSex.MALE
            return
        with open(ploidy_json) as fh:
            data = json.load(fh)
            sex = data["sex"]
            if sex == "female":
                self.sample_sex = SampleSex.FEMALE
            elif sex == "male":
                self.sample_sex = SampleSex.MALE
            else:
                self.sample_sex = SampleSex.UNKNOWN
