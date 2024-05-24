"""
DNAscope alignment and variant calling
"""

import argparse
import itertools
import multiprocessing as mp
import pathlib
import shlex
import shutil
import sys
from typing import Any, Callable, List, Optional

import packaging.version

from argh import arg

from . import command_strings as cmds
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
from .util import (
    __version__,
    check_version,
    library_preloaded,
    logger,
    path_arg,
    tmp,
)

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


def align_inputs(
    run: Callable[[str], None],
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    cores: int = mp.cpu_count(),
    dry_run: bool = False,
    collate_align: bool = False,
    skip_version_check: bool = False,
    bam_format: bool = False,
    bwa_args: str = "-K 100000000",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    input_ref: Optional[pathlib.Path] = None,
    **_kwargs: Any,
) -> List[pathlib.Path]:
    """Align input BAM/CRAM/uBAM/uCRAM files with bwa"""
    if not skip_version_check:
        for cmd, min_version in ALN_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    res: List[pathlib.Path] = []
    suffix = "bam" if bam_format else "cram"
    for i, input_aln in enumerate(sample_input):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_bwa_sorted_{i}.{suffix}")
        )
        rg_lines = cmds.get_rg_lines(input_aln, dry_run)
        rg_header = tmp_dir.joinpath(f"input_{i}.hdr")
        with open(rg_header, "w", encoding="utf-8") as rg_fh:
            for line in rg_lines:
                print(line, file=rg_fh)

        run(
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
                util_sort_args=util_sort_args,
            )
        )
        res.append(out_aln)
    return res


def align_fastq(
    run: Callable[[str], None],
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int = mp.cpu_count(),
    r1_fastq: Optional[List[pathlib.Path]] = None,
    r2_fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    skip_version_check: bool = False,
    bam_format: bool = False,
    bwa_args: str = "-K 100000000",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    **_kwargs: Any,
) -> List[pathlib.Path]:
    """Align fastq files to the reference genome using bwa"""
    res: List[pathlib.Path] = []
    if r1_fastq is None and readgroups is None:
        return res
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

    suffix = "bam" if bam_format else "cram"
    r2_fastq = [] if r2_fastq is None else r2_fastq
    for i, (r1, r2, rg) in enumerate(
        itertools.zip_longest(r1_fastq, r2_fastq, readgroups)
    ):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_bwa_sorted_fq_{i}.{suffix}")
        )
        run(
            cmds.cmd_fastq_bwa(
                out_aln,
                r1,
                r2,
                rg,
                reference,
                model_bundle,
                cores,
                unzip,
                bwa_args,
                util_sort_args,
            )
        )
        res.append(out_aln)
    return res


def dedup_and_metrics(
    run: Callable[[str], None],
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
) -> List[pathlib.Path]:
    """Perform dedup and metrics collection"""
    suffix = "bam" if bam_format else "cram"

    # Create the metrics directory
    metrics_dir = pathlib.Path(str(output_vcf).replace(".vcf.gz", "_metrics"))
    if not dry_run:
        metrics_dir.mkdir(exist_ok=True)

    # LocusCollector and Metrics
    out_score = metrics_dir.joinpath("score.txt.gz")
    is_metrics = metrics_dir.joinpath("insert_size.txt")
    mqbc_metrics = metrics_dir.joinpath("mean_qual_by_cycle.txt")
    bdbc_metrics = metrics_dir.joinpath("base_distribution_by_cycle.txt")
    qualdist_metrics = metrics_dir.joinpath("qual_distribution.txt")
    as_metrics = metrics_dir.joinpath("alignment_stat.txt")
    coverage_metrics = metrics_dir.joinpath("coverage")

    # WES metrics
    hs_metrics = metrics_dir.joinpath("hybrid-selection.txt")

    # WGS metrics
    wgs_metrics = metrics_dir.joinpath("wgs.txt")
    gc_metrics = metrics_dir.joinpath("gc_bias.txt")
    gc_summary = metrics_dir.joinpath("gc_bias_summary.txt")

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
    if assay == "WES" and not bed:
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))

    driver.add_algo(MeanQualityByCycle(mqbc_metrics))
    driver.add_algo(BaseDistributionByCycle(bdbc_metrics))
    driver.add_algo(QualDistribution(qualdist_metrics))
    driver.add_algo(AlignmentStat(as_metrics))
    if assay == "WGS":
        driver.add_algo(GCBias(gc_metrics, summary=gc_summary))

    run(shlex.join(driver.build_cmd()))

    if duplicate_marking == "none":
        return sample_input

    # Dedup
    deduped = pathlib.Path(
        str(output_vcf).replace(".vcf.gz", f"_deduped.{suffix}")
    )
    dedup_metrics = pathlib.Path(
        str(output_vcf).replace(".vcf.gz", "_dedup_metrics.txt")
    )
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
    run(shlex.join(driver.build_cmd()))

    # Run HsMetricAlgo after duplicate marking to account for duplicate reads
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=[deduped] if deduped else sample_input,
        interval=bed,
    )
    if assay == "WES" and bed:
        driver.add_algo(HsMetricAlgo(hs_metrics, bed, bed))
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        run(shlex.join(driver.build_cmd()))

    # Run WgsMetricsAlgo after duplicate marking to account for duplicate reads
    if assay == "WGS":
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        driver.add_algo(WgsMetricsAlgo(wgs_metrics, include_unpaired="true"))
        driver.add_algo(CoverageMetrics(coverage_metrics))
        run(shlex.join(driver.build_cmd()))
    return [deduped]


def multiqc(
    run: Callable[[str], None],
    output_vcf: pathlib.Path,
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> int:
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
            return 1

    metrics_dir = pathlib.Path(str(output_vcf).replace(".vcf.gz", "_metrics"))
    run(
        cmds.cmd_multiqc(
            metrics_dir,
            metrics_dir,
            f"Generated by the Sentieon-CLI version {__version__}",
        )
    )
    return 0


def call_variants(
    run: Callable[[str], None],
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
) -> int:
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
    run(shlex.join(driver.build_cmd()))

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
    run(shlex.join(driver.build_cmd()))

    # Remove the tmp_vcf
    tmp_vcf_idx = pathlib.Path(str(tmp_vcf) + ".tbi")
    tmp_vcf_idx.unlink(missing_ok=True)
    tmp_vcf.unlink()

    # Genotype gVCFs
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
        run(shlex.join(driver.build_cmd()))

    # Call SVs
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
        run(shlex.join(driver.build_cmd()))
    return 0


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
    "--util_sort_args",
    help="Extra arguments for sentieon util sort",
)
@arg(
    "--skip-version-check",
    help=argparse.SUPPRESS,
    action="store_true",
)
def dnascope(
    output_vcf: pathlib.Path,  # pylint: disable=W0613
    reference: Optional[pathlib.Path] = None,
    sample_input: Optional[List[pathlib.Path]] = None,
    r1_fastq: Optional[List[pathlib.Path]] = None,
    r2_fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    model_bundle: Optional[pathlib.Path] = None,
    dbsnp: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    bed: Optional[pathlib.Path] = None,
    interval_padding: int = 0,  # pylint: disable=W0613
    cores: int = mp.cpu_count(),  # pylint: disable=W0613
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
    bwa_args: str = "-K 100000000",  # pylint: disable=W0613
    util_sort_args: str = (
        "--cram_write_options version=3.0,compressor=rans"
    ),  # pylint: disable=W0613
    skip_version_check: bool = False,  # pylint: disable=W0613
    **kwargs: str,
):
    """
    Run the DNAscope pipeline
    """
    assert reference
    assert sample_input or (r1_fastq and readgroups)
    assert model_bundle

    logger.setLevel(kwargs["loglevel"])
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

    tmp_dir_str = tmp()
    tmp_dir = pathlib.Path(tmp_dir_str)  # type: ignore  # pylint: disable=W0641  # noqa: E501

    if dry_run:
        run = print  # type: ignore  # pylint: disable=W0641
    else:
        from .runner import run  # type: ignore[assignment]  # noqa: F401

    sample_input = sample_input if sample_input else []
    if align or collate_align:
        sample_input = align_inputs(**locals())
    sample_input.extend(align_fastq(**locals()))

    deduped = dedup_and_metrics(**locals())  # pylint: disable=W0641

    # Remove the bwa output before duplicate marking
    if duplicate_marking != "none" and not dry_run:
        for aln in sample_input:
            for idx_suffix in (".bai", ".crai"):
                idx = pathlib.Path(str(aln) + idx_suffix)
                idx.unlink(missing_ok=True)
            aln.unlink()

    if not skip_small_variants:
        res = call_variants(**locals())
        if res != 0:
            logger.error("Variant calling failed")
            return

    if not skip_multiqc:
        _res = multiqc(**locals())

    shutil.rmtree(tmp_dir_str)
