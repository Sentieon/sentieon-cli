"""
DNAscope alignment and variant calling
"""

import argparse
import itertools
import multiprocessing as mp
import os
import pathlib
import shlex
import shutil
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import packaging.version

from argh import arg, CommandError

from . import command_strings as cmds
from .dag import DAG
from .driver import (
    AlignmentStat,
    BaseDistributionByCycle,
    CoverageMetrics,
    Dedup,
    DNAModelApply,
    DNAscope,
    Driver,
    GCBias,
    GVCFtyper,
    HsMetricAlgo,
    InsertSizeMetricAlgo,
    MeanQualityByCycle,
    LocusCollector,
    QualDistribution,
    SVSolver,
    WgsMetricsAlgo,
)
from .executor import DryRunExecutor, LocalExecutor
from .job import Job
from .logging import get_logger
from .scheduler import ThreadScheduler
from .util import (
    __version__,
    check_version,
    find_numa_nodes,
    library_preloaded,
    path_arg,
    split_numa_nodes,
    tmp,
    total_memory,
)

logger = get_logger(__name__)


ALN_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
    "samtools": packaging.version.Version("1.16"),
}

FQ_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}

VARIANTS_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}

MULTIQC_MIN_VERSION = {
    "multiqc": packaging.version.Version("1.18"),
}


def set_bwt_max_mem(
    bwt_max_mem_arg: Optional[str],
    total_input_size: Optional[int],
    n_alignment_jobs: int = 1,
):
    """Set the bwt_max_mem environment variable"""
    if bwt_max_mem_arg:
        os.environ["bwt_max_mem"] = bwt_max_mem_arg
        return

    total_input_size = total_input_size if total_input_size else 0
    total_mem = total_memory()
    total_mem_gb = total_mem / (1024.0**3)
    align_mem_gb = (
        total_mem_gb - 4 - total_input_size / (1024.0**3) * 2.3
    )  # some memory for other system processes
    bwa_mem_gb = max(
        int((align_mem_gb / n_alignment_jobs) - 6), 0
    )  # some memory for other alignment processes
    logger.debug("Setting bwt_max_mem to: %sG", bwa_mem_gb)
    os.environ["bwt_max_mem"] = f"{bwa_mem_gb}G"


def check_shm(
    sample_input: Optional[List[pathlib.Path]],
    r1_fastq: Optional[List[pathlib.Path]],
    r2_fastq: Optional[List[pathlib.Path]],
    n_alignment_jobs: int,
) -> Optional[int]:
    """Check if /dev/shm can be used for temporary files"""
    # Find the size of the largest input
    total_input_size = 0
    if sample_input:
        total_input_size = sum([x.stat().st_size for x in sample_input])

    r1_fastq_l = r1_fastq if r1_fastq else []
    r2_fastq_l = r2_fastq if r2_fastq else []

    for r1, r2 in itertools.zip_longest(r1_fastq_l, r2_fastq_l):
        total = sum(
            [
                fq.stat().st_size
                for fq in (r1, r2)
                if isinstance(fq, pathlib.Path)
            ]
        )
        total_input_size += total

    shm_free = shutil.disk_usage("/dev/shm").free
    total_mem = total_memory()

    if (
        shm_free > total_input_size * 2.3
        and total_mem > total_input_size + 40 * 1024**3 * n_alignment_jobs
    ):
        logger.debug("Using /dev/shm for temporary files")
        return total_input_size
    return None


def align_inputs(
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    cores: int = mp.cpu_count(),
    duplicate_marking: str = "markdup",
    dry_run: bool = False,
    collate_align: bool = False,
    skip_version_check: bool = False,
    bam_format: bool = False,
    bwa_args: str = "",
    bwa_k_arg: str = "100000000",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    input_ref: Optional[pathlib.Path] = None,
    **_kwargs: Any,
) -> Tuple[List[pathlib.Path], Set[Job], Job]:
    """Align input BAM/CRAM/uBAM/uCRAM files with bwa"""
    if not skip_version_check:
        for cmd, min_version in ALN_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    res: List[pathlib.Path] = []
    suffix = "bam" if bam_format else "cram"
    if duplicate_marking != "none":
        suffix = "bam"
        util_sort_args += " --bam_compression 1 "
    jobs = set()
    align_outputs: List[pathlib.Path] = []
    for i, input_aln in enumerate(sample_input):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_bwa_sorted_{i}.{suffix}")
        )
        if duplicate_marking != "none":
            out_aln = tmp_dir.joinpath(f"bwa_sorted_{i}.{suffix}")
        rg_lines = cmds.get_rg_lines(input_aln, dry_run)
        rg_header = tmp_dir.joinpath(f"input_{i}.hdr")
        with open(rg_header, "w", encoding="utf-8") as rg_fh:
            for line in rg_lines:
                print(line, file=rg_fh)
        job = Job(
            cmds.cmd_samtools_fastq_bwa(
                out_aln,
                input_aln,
                reference,
                model_bundle,
                cores,
                rg_header,
                input_ref,
                collate=collate_align,
                bwa_args=bwa_args,
                bwa_k=bwa_k_arg,
                util_sort_args=util_sort_args,
            ),
            f"bam-align-{i}",
            cores,
        )
        res.append(out_aln)
        jobs.add(job)
        align_outputs.append(out_aln)
        align_outputs.append(pathlib.Path(str(out_aln) + ".bai"))
        if not bam_format:
            align_outputs.append(pathlib.Path(str(out_aln) + ".crai"))

    # Create an unscheduled job to remove the aligned inputs
    rm_job = Job(
        shlex.join(["rm"] + [str(x) for x in align_outputs]),
        "rm-bam-aln",
        0,
        True,
    )

    return (res, jobs, rm_job)


def align_fastq(
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    numa_nodes: List[str],
    cores: int = mp.cpu_count(),
    duplicate_marking: str = "markdup",
    r1_fastq: Optional[List[pathlib.Path]] = None,
    r2_fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    skip_version_check: bool = False,
    bam_format: bool = False,
    bwa_args: str = "",
    bwa_k_arg: str = "100000000",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    **_kwargs: Any,
) -> Tuple[List[pathlib.Path], Set[Job], Optional[Job]]:
    """Align fastq files to the reference genome using bwa"""
    res: List[pathlib.Path] = []
    jobs: Set[Job] = set()
    if r1_fastq is None and readgroups is None:
        return (res, jobs, None)
    if (not r1_fastq or not readgroups) or (len(r1_fastq) != len(readgroups)):
        logger.error(
            "The number of readgroups does not equal the number of fastq files"
        )
        sys.exit(1)

    if not skip_version_check:
        for cmd, min_version in FQ_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    unzip = "igzip"
    if not shutil.which(unzip):
        logger.warning(
            "igzip is recommended for decompression, but is not available. "
            "Falling back to gzip."
        )
        unzip = "gzip"

    n_alignment_jobs = max(1, len(numa_nodes))
    suffix = "bam" if bam_format else "cram"
    if duplicate_marking != "none":
        suffix = "bam"
        util_sort_args += " --bam_compression 1 "
    r2_fastq = [] if r2_fastq is None else r2_fastq
    align_outputs: List[pathlib.Path] = []
    for i, (r1, r2, rg) in enumerate(
        itertools.zip_longest(r1_fastq, r2_fastq, readgroups)
    ):
        for j in range(n_alignment_jobs):
            out_aln = pathlib.Path(
                str(output_vcf).replace(
                    ".vcf.gz", f"_bwa_sorted_fq_{i}_{j}.{suffix}"
                )
            )
            if duplicate_marking != "none":
                out_aln = tmp_dir.joinpath(f"bwa_sorted_fq_{i}_{j}.{suffix}")
            numa = numa_nodes[j] if len(numa_nodes) > 1 else None
            split = f"{j}/{n_alignment_jobs}" if len(numa_nodes) > 1 else None
            split_cores = max(1, int(cores / n_alignment_jobs))
            job = Job(
                cmds.cmd_fastq_bwa(
                    out_aln,
                    r1,
                    r2,
                    rg,
                    reference,
                    model_bundle,
                    split_cores,
                    unzip,
                    bwa_args,
                    bwa_k_arg,
                    util_sort_args,
                    numa,
                    split,
                ),
                f"bam-align-{i}-{j}",
                split_cores,
                resources={f"node{j}": 1},
            )
            res.append(out_aln)
            jobs.add(job)
            align_outputs.append(out_aln)
            align_outputs.append(pathlib.Path(str(out_aln) + ".bai"))
            if not bam_format:
                align_outputs.append(pathlib.Path(str(out_aln) + ".crai"))

    # Create an unscheduled job to remove the aligned inputs
    rm_job = Job(
        shlex.join(["rm"] + [str(x) for x in align_outputs]),
        "rm-fq-aln",
        0,
        True,
    )

    return (res, jobs, rm_job)


def dedup_and_metrics(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    bed: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),
    duplicate_marking: str = "markdup",
    assay: str = "WGS",
    consensus: bool = False,
    dry_run: bool = False,
    bam_format: bool = False,
    cram_write_options: str = "version=3.0,compressor=rans",
    **_kwargs: Any,
) -> Tuple[
    List[pathlib.Path], Job, Optional[Job], Optional[Job], Optional[Job]
]:
    """Perform dedup and metrics collection"""
    suffix = "bam" if bam_format else "cram"

    # Create the metrics directory
    sample_name = output_vcf.name.replace(".vcf.gz", "")
    metric_base = sample_name + ".txt"
    metrics_dir = pathlib.Path(str(output_vcf).replace(".vcf.gz", "_metrics"))
    if not dry_run:
        metrics_dir.mkdir(exist_ok=True)

    # LocusCollector and Metrics
    out_score = metrics_dir.joinpath(metric_base + ".score.txt.gz")
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

    # WES metrics
    hs_metrics = metrics_dir.joinpath(metric_base + ".hybrid-selection.txt")

    # WGS metrics
    wgs_metrics = metrics_dir.joinpath(metric_base + ".wgs.txt")
    wgs_metrics_tmp = metrics_dir.joinpath(metric_base + ".wgs.txt.tmp")
    gc_metrics = metrics_dir.joinpath(metric_base + ".gc_bias.txt")
    gc_summary = metrics_dir.joinpath(metric_base + ".gc_bias_summary.txt")

    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
    )
    if duplicate_marking != "none":
        driver.add_algo(
            LocusCollector(
                out_score,
                consensus=consensus,
            )
        )

    # Prefer to run InsertSizeMetricAlgo after duplicate marking
    if (assay == "WES" and not bed) or duplicate_marking == "none":
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))

    driver.add_algo(MeanQualityByCycle(mqbc_metrics))
    driver.add_algo(BaseDistributionByCycle(bdbc_metrics))
    driver.add_algo(QualDistribution(qualdist_metrics))
    driver.add_algo(AlignmentStat(as_metrics))
    if assay == "WGS":
        driver.add_algo(GCBias(gc_metrics, summary=gc_summary))

    lc_job = Job(shlex.join(driver.build_cmd()), "locuscollector", cores)

    if duplicate_marking == "none":
        return (sample_input, lc_job, None, None, None)

    # Dedup
    deduped = pathlib.Path(
        str(output_vcf).replace(".vcf.gz", f"_deduped.{suffix}")
    )
    dedup_metrics = metrics_dir.joinpath(metric_base + ".dedup_metrics.txt")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
    )
    driver.add_algo(
        Dedup(
            deduped,
            out_score,
            cram_write_options=cram_write_options,
            metrics=dedup_metrics,
            rmdup=(duplicate_marking == "rmdup"),
        )
    )
    dedup_job = Job(shlex.join(driver.build_cmd()), "dedup", cores)

    # Run HsMetricAlgo after duplicate marking to account for duplicate reads
    metrics_job = None
    rehead_job = None
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=[deduped] if deduped else sample_input,
        interval=bed,
    )
    if assay == "WES" and bed:
        driver.add_algo(HsMetricAlgo(hs_metrics, bed, bed))
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        metrics_job = Job(
            shlex.join(driver.build_cmd()), "metrics", 0
        )  # Run metrics in the background

    # Run WgsMetricsAlgo after duplicate marking to account for duplicate reads
    if assay == "WGS":
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        driver.add_algo(WgsMetricsAlgo(wgs_metrics, include_unpaired="true"))
        driver.add_algo(CoverageMetrics(coverage_metrics))
        metrics_job = Job(
            shlex.join(driver.build_cmd()), "metrics", 0
        )  # Run metrics in the background

        # Rehead WGS metrics so they are recognized by MultiQC
        rehead_job = Job(
            cmds.rehead_wgsmetrics(wgs_metrics, wgs_metrics_tmp),
            "Rehead metrics",
            0,
        )
    return ([deduped], lc_job, dedup_job, metrics_job, rehead_job)


def multiqc(
    output_vcf: pathlib.Path,
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> Optional[Job]:
    """Run MultiQC on the metrics files"""

    if not skip_version_check:
        if not all(
            [
                check_version(cmd, min_version)
                for (cmd, min_version) in MULTIQC_MIN_VERSION.items()
            ]
        ):
            logger.warning(
                "Skipping MultiQC. MultiQC version %s or later not found",
                MULTIQC_MIN_VERSION["multiqc"],
            )
            return None

    metrics_dir = pathlib.Path(str(output_vcf).replace(".vcf.gz", "_metrics"))
    multiqc_job = Job(
        cmds.cmd_multiqc(
            metrics_dir,
            metrics_dir,
            f"Generated by the Sentieon-CLI version {__version__}",
        ),
        "multiqc",
        0,
    )
    return multiqc_job


def call_variants(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    deduped: List[pathlib.Path],
    model_bundle: pathlib.Path,
    dbsnp: Optional[pathlib.Path] = None,
    bed: Optional[pathlib.Path] = None,
    interval_padding: int = 0,
    cores: int = mp.cpu_count(),
    pcr_free: bool = False,
    gvcf: bool = False,
    skip_svs: bool = False,
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> Tuple[Job, Job, Job, Optional[Job], Optional[Job], Optional[Job]]:
    """Call SNVs, indels, and SVs using DNAscope"""
    if not skip_version_check:
        for cmd, min_version in VARIANTS_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    out_gvcf = pathlib.Path(str(output_vcf).replace(".vcf.gz", ".g.vcf.gz"))
    out_svs_tmp = pathlib.Path(
        str(output_vcf).replace(".vcf.gz", "_svs_tmp.vcf.gz")
    )
    out_svs = pathlib.Path(str(output_vcf).replace(".vcf.gz", "_svs.vcf.gz"))
    emit_mode = "gvcf" if gvcf else "variant"
    pcr_indel_model = "NONE" if pcr_free else "CONSERVATIVE"
    model = model_bundle.joinpath("dnascope.model")

    # Call variants with DNAscope
    emit_mode = "variant"
    ds_out = output_vcf
    tmp_vcf = pathlib.Path(str(output_vcf).replace(".vcf.gz", "_tmp.vcf.gz"))
    if gvcf:
        emit_mode = "gvcf"
        ds_out = out_gvcf
        tmp_vcf = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", "_tmp.g.vcf.gz")
        )
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=deduped,
        interval=bed,
        interval_padding=interval_padding,
    )
    driver.add_algo(
        DNAscope(
            tmp_vcf,
            dbsnp=dbsnp,
            emit_mode=emit_mode,
            pcr_indel_model=pcr_indel_model,
            model=model,
        )
    )
    if not skip_svs:
        driver.add_algo(
            DNAscope(
                out_svs_tmp,
                dbsnp=dbsnp,
                var_type="BND",
            )
        )
    call_job = Job(shlex.join(driver.build_cmd()), "variant-calling", cores)

    # Genotyping and filtering with DNAModelApply
    driver = Driver(
        reference=reference,
        thread_count=cores,
        interval=bed,
    )
    driver.add_algo(
        DNAModelApply(
            model,
            tmp_vcf,
            ds_out,
        )
    )
    apply_job = Job(shlex.join(driver.build_cmd()), "model-apply", cores)

    # Remove the tmp_vcf
    rm_cmd = ["rm", str(tmp_vcf), str(tmp_vcf) + ".tbi"]
    rm_job = Job(shlex.join(rm_cmd), "rm-tmp-vcf", 0, True)

    # Genotype gVCFs
    gvcftyper_job = None
    if gvcf:
        driver = Driver(
            reference=reference,
            thread_count=cores,
            interval=bed,
        )
        driver.add_algo(
            GVCFtyper(
                output=output_vcf,
                vcf=out_gvcf,
            )
        )
        gvcftyper_job = Job(shlex.join(driver.build_cmd()), "gvcftyper", cores)

    # Call SVs
    svsolver_job = None
    sv_rm_job = None
    if not skip_svs:
        driver = Driver(
            reference=reference,
            thread_count=cores,
            interval=bed,
        )
        driver.add_algo(
            SVSolver(
                output=out_svs,
                vcf=out_svs_tmp,
            )
        )
        svsolver_job = Job(shlex.join(driver.build_cmd()), "svsolver")
        sv_rm_job = Job(
            shlex.join(["rm", str(out_svs_tmp), str(out_svs_tmp) + ".tbi"]),
            "rm-tmp-sv",
            0,
            True,
        )

    return (
        call_job,
        apply_job,
        rm_job,
        gvcftyper_job,
        svsolver_job,
        sv_rm_job,
    )


@arg(
    "-r",
    "--reference",
    required=True,
    help="fasta for reference genome",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "-i",
    "--sample-input",
    nargs="*",
    help="sample BAM or CRAM file",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--r1-fastq",
    nargs="*",
    help="Sample R1 fastq files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--r2-fastq",
    nargs="*",
    help="Sample R2 fastq files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--readgroups",
    nargs="*",
    help="Readgroup information for the fastq files",
)
@arg(
    "-m",
    "--model-bundle",
    help="The model bundle file",
    required=True,
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "output-vcf",
    help="Output VCF File. The file name must end in .vcf.gz",
    type=path_arg(),
)
@arg(
    "-d",
    "--dbsnp",
    help="dbSNP vcf file Supplying this file will annotate variants with \
         their dbSNP refSNP ID numbers.",
    default=None,
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "-b",
    "--bed",
    help="Region BED file. Supplying this file will limit variant calling \
    to the intervals inside the BED file.",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--interval_padding",
    help="Amount to pad all intervals.",
)
@arg(
    "-t",
    "--cores",
    help="Number of threads/processes to use. %(default)s",
)
@arg(
    "--pcr-free",
    help="Use arguments for PCR-free data processing",
    action="store_true",
)
@arg(
    "-g",
    "--gvcf",
    help="Generate a gVCF output file along with the VCF."
    " (default generates only the VCF)",
    action="store_true",
)
@arg(
    "--duplicate-marking",
    help="Options for duplicate marking.",
    choices=["markdup", "rmdup", "none"],
)
@arg(
    "--assay",
    help="The type of assay, WGS or WES.",
    choices=["WGS", "WES"],
)
@arg(
    "--consensus",
    help="Generate consensus reads during dedup",
    action="store_true",
)
@arg(
    "--dry-run",
    help="Print the commands without running them.",
)
@arg(
    "--skip-small-variants",
    help="Skip small variant (SNV/indel) calling",
)
@arg(
    "--skip-svs",
    help="Skip SV calling",
)
@arg(
    "--skip-multiqc",
    help="Skip multiQC report generation",
)
@arg(
    "--align",
    help="Align the reads in the input uBAM/uCRAM file to the reference "
    "genome. Assumes paired reads are collated in the input file.",
    action="store_true",
)
@arg(
    "--collate-align",
    help="Collate and align the reads in the input BAM/CRAM file to the "
    "reference genome. Suitable for coordinate-sorted BAM/CRAM input.",
    action="store_true",
)
@arg(
    "--input_ref",
    help="Used to decode the input alignment file. Required if the input file"
    " is in the CRAM/uCRAM formats",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--bam_format",
    help="Use the BAM format instead of CRAM for output aligned files",
    action="store_true",
)
@arg(
    "--bwa_args",
    help="Extra arguments for sentieon bwa",
)
@arg(
    "--bwa_k",
    help="The '-K' argument in bwa",
)
@arg(
    "--util_sort_args",
    help="Extra arguments for sentieon util sort",
)
@arg(
    "--skip-version-check",
    help=argparse.SUPPRESS,
    action="store_true",
)
# Manually set `bwt_max_mem`
@arg(
    "--bwt-max-mem",
    help=argparse.SUPPRESS,
)
# Do not use /dev/shm, even on high memory machines
@arg(
    "--no-ramdisk",
    help=argparse.SUPPRESS,
    action="store_true",
)
# Do not split alignment into multiple jobs
@arg(
    "--no-split-alignment",
    help=argparse.SUPPRESS,
    action="store_true",
)
def dnascope(
    output_vcf: pathlib.Path,
    reference: Optional[pathlib.Path] = None,
    sample_input: Optional[List[pathlib.Path]] = None,
    r1_fastq: Optional[List[pathlib.Path]] = None,
    r2_fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    model_bundle: Optional[pathlib.Path] = None,
    dbsnp: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    bed: Optional[pathlib.Path] = None,
    interval_padding: int = 0,  # pylint: disable=W0613
    cores: int = mp.cpu_count(),
    pcr_free: bool = False,  # pylint: disable=W0613
    gvcf: bool = False,  # pylint: disable=W0613
    duplicate_marking: str = "markdup",
    assay: str = "WGS",
    consensus: bool = False,  # pylint: disable=W0613
    dry_run: bool = False,
    skip_small_variants: bool = False,
    skip_svs: bool = False,
    skip_multiqc: bool = False,
    align: bool = False,
    collate_align: bool = False,
    input_ref: Optional[pathlib.Path] = None,
    bam_format: bool = False,  # pylint: disable=W0613
    bwa_args: str = "",  # pylint: disable=W0613
    bwa_k: int = 100000000,
    util_sort_args: str = (
        "--cram_write_options version=3.0,compressor=rans"
    ),  # pylint: disable=W0613
    skip_version_check: bool = False,  # pylint: disable=W0613
    bwt_max_mem: Optional[str] = None,
    no_ramdisk: bool = False,
    no_split_alignment: bool = False,
    **kwargs: str,
):
    """
    Run the DNAscope pipeline
    """
    assert reference
    assert sample_input or (r1_fastq and readgroups)
    assert model_bundle
    assert str(output_vcf).endswith(".vcf.gz")
    bwa_k_arg = str(bwa_k)

    assert logger.parent
    logger.parent.setLevel(kwargs["loglevel"])
    logger.info("Starting sentieon-cli version: %s", __version__)

    if not library_preloaded("libjemalloc.so"):
        logger.warning(
            "jemalloc is recommended, but is not preloaded. See "
            "https://support.sentieon.com/appnotes/jemalloc/"
        )

    if bed is None:
        if assay == "WES":
            logger.warning(
                "A BED file is recommended with WES assays to restrict "
                "variant calling to target regions."
            )
        else:
            logger.info(
                "A BED file is recommended to avoid variant calling across "
                "decoy and unplaced contigs."
            )

    # split large alignment tasks into smaller tasks on large machines
    n_alignment_jobs = 1
    numa_nodes: List[str] = []
    if not no_split_alignment:
        numa_nodes = find_numa_nodes()
        if len(numa_nodes) > 0 and cores / len(numa_nodes) > 48:
            numa_nodes = split_numa_nodes(numa_nodes)
        if cores > 32 and numa_nodes:
            n_alignment_jobs = len(numa_nodes)
        else:
            numa_nodes = []

    total_input_size = None
    if not no_ramdisk:
        total_input_size = check_shm(
            sample_input, r1_fastq, r2_fastq, n_alignment_jobs
        )
    if total_input_size:
        os.environ["SENTIEON_TMPDIR"] = "/dev/shm"

    set_bwt_max_mem(bwt_max_mem, total_input_size, n_alignment_jobs)

    tmp_dir_str = tmp()
    tmp_dir = pathlib.Path(tmp_dir_str)  # type: ignore  # pylint: disable=W0641  # noqa: E501

    logger.info("Building the DAG")
    dag = DAG()

    align_jobs: Set[Job] = set()
    sample_input = sample_input if sample_input else []
    bam_rm_job = None
    if align or collate_align:
        sample_input, align_jobs, bam_rm_job = align_inputs(**locals())
        for job in align_jobs:
            dag.add_job(job)
    aligned_fastq, align_fastq_jobs, fq_rm_job = align_fastq(**locals())
    for job in align_fastq_jobs:
        dag.add_job(job)
    sample_input.extend(aligned_fastq)

    deduped, lc_job, dedup_job, metrics_job, rehead_job = dedup_and_metrics(
        **locals()
    )  # pylint: disable=W0641
    dag.add_job(lc_job, align_jobs.union(align_fastq_jobs))
    if dedup_job:
        dag.add_job(dedup_job, {lc_job})
        if metrics_job:
            dag.add_job(metrics_job, {dedup_job})
            if rehead_job:
                dag.add_job(rehead_job, {metrics_job})
        if bam_rm_job:
            dag.add_job(bam_rm_job, {dedup_job})
        if fq_rm_job:
            dag.add_job(fq_rm_job, {dedup_job})

    if not skip_small_variants:
        (
            call_job,
            apply_job,
            rm_job,
            gvcftyper_job,
            svsolver_job,
            sv_rm_job,
        ) = call_variants(**locals())
        call_dependencies: Set[Job] = set()
        if dedup_job:
            call_dependencies.add(dedup_job)
        else:
            call_dependencies.update(align_jobs)
            call_dependencies.update(align_fastq_jobs)
        dag.add_job(call_job, call_dependencies)
        dag.add_job(apply_job, {call_job})
        dag.add_job(rm_job, {apply_job})
        if gvcftyper_job:
            dag.add_job(gvcftyper_job, {apply_job})
        if svsolver_job and sv_rm_job:
            dag.add_job(svsolver_job, {call_job})
            dag.add_job(sv_rm_job, {svsolver_job})

    if not skip_multiqc:
        multiqc_job = multiqc(**locals())
        multiqc_dependencies: Set[Job] = set()
        multiqc_dependencies.add(lc_job)
        if metrics_job:
            multiqc_dependencies.add(metrics_job)
        if rehead_job:
            multiqc_dependencies.add(rehead_job)

        if multiqc_job:
            dag.add_job(multiqc_job, multiqc_dependencies)

    logger.debug("Creating the scheduler")
    resources: Dict[str, int] = {}
    for i in range(len(numa_nodes)):
        resources[f"node{i}"] = 1
    scheduler = ThreadScheduler(
        dag,
        cores,
        resources,
    )

    logger.debug("Creating the executor")
    Executor = DryRunExecutor if dry_run else LocalExecutor
    executor = Executor(scheduler)

    logger.info("Starting execution")
    executor.execute()
    shutil.rmtree(tmp_dir_str)

    if executor.jobs_with_errors:
        raise CommandError("Execution failed")

    if len(dag.waiting_jobs) > 0 or len(dag.ready_jobs) > 0:
        raise CommandError(
            "The DAG has some unexecuted jobs\n"
            f"Waiting jobs: {dag.waiting_jobs}\n"
            f"Ready jobs: {dag.ready_jobs}\n"
        )
