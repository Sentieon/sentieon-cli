"""
Functionality for the DNAscope LongRead pipeline
"""

import argparse
import multiprocessing as mp
import pathlib
import shlex
import shutil
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import packaging.version

from argh import arg, CommandError
from importlib_resources import files

from . import command_strings as cmds
from .dag import DAG
from .driver import (
    Driver,
    DNAscope,
    DNAscopeHP,
    DNAModelApply,
    LongReadSV,
    ReadWriter,
    RepeatModel,
    VariantPhaser,
)
from .executor import DryRunExecutor, LocalExecutor
from .job import Job
from .logging import get_logger
from .scheduler import ThreadScheduler
from .util import (
    __version__,
    check_version,
    path_arg,
    tmp,
    library_preloaded,
)

logger = get_logger(__name__)


TOOL_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308.01"),
    "bcftools": packaging.version.Version("1.10"),
    "bedtools": None,
}

SV_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}

ALN_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
    "samtools": packaging.version.Version("1.16"),
}

FQ_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}

MOSDEPTH_MIN_VERSIONS = {
    "mosdepth": packaging.version.Version("0.2.6"),
}

MERGE_MIN_VERSIONS = {
    "sentieon driver": None,
}

PBSV_MIN_VERSIONS = {
    "pbsv": packaging.version.Version("2.0.0"),
}

HIFICNV_MIN_VERSIONS = {
    "hificnv": packaging.version.Version("1.0.0"),
}


def align_inputs(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    cores: int = mp.cpu_count(),
    dry_run: bool = False,
    skip_version_check: bool = False,
    bam_format: bool = False,
    fastq_taglist: str = "*",
    minimap2_args: str = "-Y",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    input_ref: Optional[pathlib.Path] = None,
    **_kwargs: Any,
) -> Tuple[List[pathlib.Path], Set[Job]]:
    """
    Align reads to the reference genome using minimap2
    """
    if not skip_version_check:
        for cmd, min_version in ALN_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    res: List[pathlib.Path] = []
    suffix = "bam" if bam_format else "cram"
    sample_name = output_vcf.name.replace(".vcf.gz", "")
    realign_jobs = set()
    for i, input_aln in enumerate(sample_input):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_mm2_sorted_{i}.{suffix}")
        )
        rg_lines = cmds.get_rg_lines(
            input_aln,
            dry_run,
        )

        realign_jobs.add(
            Job(
                cmds.cmd_samtools_fastq_minimap2(
                    out_aln,
                    input_aln,
                    reference,
                    model_bundle,
                    cores,
                    rg_lines,
                    sample_name,
                    input_ref,
                    fastq_taglist,
                    minimap2_args,
                    util_sort_args,
                ),
                f"bam-realign-{i}",
                cores,
            )
        )
        res.append(out_aln)
    return (res, realign_jobs)


def align_fastq(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int = mp.cpu_count(),
    fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    skip_version_check: bool = False,
    bam_format: bool = False,
    minimap2_args: str = "-Y",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    **_kwargs: Any,
) -> Tuple[List[pathlib.Path], Set[Job]]:
    """
    Align fastq to the reference genome using minimap2
    """
    res: List[pathlib.Path] = []
    if fastq is None and readgroups is None:
        return (res, set())
    if (not fastq or not readgroups) or (len(fastq) != len(readgroups)):
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
    align_jobs = set()
    for i, (fq, rg) in enumerate(zip(fastq, readgroups)):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_mm2_sorted_fq_{i}.{suffix}")
        )
        align_jobs.add(
            Job(
                cmds.cmd_fastq_minimap2(
                    out_aln,
                    fq,
                    rg,
                    reference,
                    model_bundle,
                    cores,
                    unzip,
                    minimap2_args,
                    util_sort_args,
                ),
                "align-{i}",
                cores,
            )
        )
        res.append(out_aln)
    return (res, align_jobs)


def call_variants(
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    dbsnp: Optional[pathlib.Path] = None,
    bed: Optional[pathlib.Path] = None,
    haploid_bed: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),
    gvcf: bool = False,
    tech: str = "HiFi",
    dry_run: bool = False,
    repeat_model: Optional[pathlib.Path] = None,
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> Tuple[
    Job,
    Job,
    Job,
    Optional[Job],
    Optional[Job],
    Job,
    Optional[Job],
    Job,
    Set[Job],
    Job,
    Set[Job],
    Job,
    Job,
    Job,
    Job,
    Optional[Job],
    Optional[Job],
    Optional[Job],
    Optional[Job],
    Optional[Job],
    Optional[Job],
]:
    """
    Call SNVs and indels using the DNAscope LongRead pipeline
    """
    if haploid_bed and not bed:
        logger.error(
            "Please supply a BED file of diploid regions to distinguish "
            "haploid and diploid regions of the genome."
        )
        sys.exit(1)

    if not bed:
        logger.warning(
            "A BED file is recommended to restrict variant calling to diploid "
            "regions of the genome."
        )

    if not skip_version_check:
        for cmd, min_version in TOOL_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    # First pass - diploid calling
    diploid_gvcf_fn = tmp_dir.joinpath("out_diploid.g.vcf.gz")
    diploid_tmp_vcf = tmp_dir.joinpath("out_diploid_tmp.vcf.gz")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    if gvcf:
        driver.add_algo(
            DNAscope(
                diploid_gvcf_fn,
                dbsnp=dbsnp,
                emit_mode="gvcf",
                model=model_bundle.joinpath("gvcf_model"),
            )
        )
    driver.add_algo(
        DNAscope(
            diploid_tmp_vcf,
            dbsnp=dbsnp,
            model=model_bundle.joinpath("diploid_model"),
        )
    )
    first_calling_job = Job(
        shlex.join(driver.build_cmd()), "first-pass", cores
    )

    diploid_vcf = tmp_dir.joinpath("out_diploid.vcf.gz")
    driver = Driver(
        reference=reference,
        thread_count=cores,
    )
    driver.add_algo(
        DNAModelApply(
            model_bundle.joinpath("diploid_model"),
            diploid_tmp_vcf,
            diploid_vcf,
        )
    )
    first_modelapply_job = Job(
        shlex.join(driver.build_cmd()), "first-modelapply", cores
    )

    # Phasing and RepeatModel
    phased_bed = tmp_dir.joinpath("out_diploid_phased.bed")
    unphased_bed = tmp_dir.joinpath("out_diploid_unphased.bed")
    phased_vcf = tmp_dir.joinpath("out_diploid_phased.vcf.gz")
    phased_ext = tmp_dir.joinpath("out_diploid_phased.ext.vcf.gz")
    phased_unphased = tmp_dir.joinpath("out_diploid_phased_unphased.vcf.gz")
    phased_phased = tmp_dir.joinpath("out_diploid_phased_phased.vcf.gz")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    driver.add_algo(
        VariantPhaser(
            diploid_vcf,
            phased_vcf,
            out_bed=phased_bed,
            out_ext=phased_ext,
        )
    )
    phaser_job = Job(shlex.join(driver.build_cmd()), "variantphaser", cores)

    bcftools_subset_phased_job = None
    if tech.upper() == "ONT":
        bcftools_subset_phased_job = Job(
            f"bcftools view -T {phased_bed} {phased_vcf} \
            | sentieon util vcfconvert - {phased_phased}",
            "bcftools-subset-phased",
            0,
        )

    fai_to_bed_job = None
    if not bed:
        bed = tmp_dir.joinpath("reference.bed")
        fai_to_bed_job = Job(
            cmds.cmd_fai_to_bed(
                pathlib.Path(str(reference) + ".fai"),
                bed,
            ),
            "fai-to-bed",
            0,
        )

    bcftools_subtract_job = Job(
        cmds.cmd_bedtools_subtract(bed, phased_bed, unphased_bed),
        "bedtools-subtract",
        0,
    )

    repeatmodel_job = None
    if not repeat_model:
        repeat_model = tmp_dir.joinpath("out_repeat.model")
        driver = Driver(
            reference=reference,
            thread_count=cores,
            input=sample_input,
            interval=phased_bed,
            read_filter=[
                f"PhasedReadFilter,phased_vcf={phased_ext},phase_select=tag"
            ],  # noqa: E501
        )
        driver.add_algo(
            RepeatModel(
                repeat_model,
                phased=True,
                read_flag_mask="drop=supplementary",
            )
        )
        repeatmodel_job = Job(
            shlex.join(driver.build_cmd()), "repeatmodel", cores
        )

    bcftools_subset_unphased_job = Job(
        f"bcftools view -T {unphased_bed} {phased_vcf} \
        | sentieon util vcfconvert - {phased_unphased}",
        "bcftools-subset-unphased",
        0,
    )

    # Second pass - phased variants
    second_calling_job = set()
    for phase in (1, 2):
        hp_std_vcf = tmp_dir.joinpath(f"out_hap{phase}_nohp_tmp.vcf.gz")
        hp_vcf = tmp_dir.joinpath(f"out_hap{phase}_tmp.vcf.gz")
        driver = Driver(
            reference=reference,
            thread_count=cores,
            input=sample_input,
            interval=phased_bed,
            read_filter=[
                f"PhasedReadFilter,phased_vcf={phased_ext}"
                f",phase_select={phase}"
            ],
        )

        if tech.upper() == "HIFI":
            # ONT doesn't do DNAscope in 2nd pass.
            driver.add_algo(
                DNAscope(
                    hp_std_vcf,
                    dbsnp=dbsnp,
                    model=model_bundle.joinpath("haploid_model"),
                )
            )
        driver.add_algo(
            DNAscopeHP(
                hp_vcf,
                dbsnp=dbsnp,
                model=model_bundle.joinpath("haploid_hp_model"),
                pcr_indel_model=repeat_model,
            )
        )
        second_calling_job.add(
            Job(shlex.join(driver.build_cmd()), "second-pass", cores)
        )

    kwargs: Dict[str, str] = dict()
    kwargs["gvcf_combine_py"] = str(
        files("sentieon_cli.scripts").joinpath("gvcf_combine.py")
    )
    kwargs["vcf_mod_py"] = str(
        files("sentieon_cli.scripts").joinpath("vcf_mod.py")
    )

    patch_vcfs = [tmp_dir.joinpath(f"out_hap{i}_patch.vcf.gz") for i in (1, 2)]
    haploid_patch_job = Job(
        cmds.cmd_pyexec_vcf_mod_haploid_patch(
            str(patch_vcfs[0]),
            str(patch_vcfs[1]),
            f"{tmp_dir}/out_hap%d_%stmp.vcf.gz",
            tech,
            str(phased_phased),
            cores,
            kwargs,
        ),
        "patch",
        cores,
    )

    # apply trained model to the patched vcfs.
    hap_vcfs = [tmp_dir.joinpath(f"out_hap{i}.vcf.gz") for i in (1, 2)]
    second_modelapply_job = set()
    for patch_vcf, hap_vcf in zip(patch_vcfs, hap_vcfs):
        driver = Driver(
            reference=reference,
            thread_count=cores,
        )
        driver.add_algo(
            DNAModelApply(
                model_bundle.joinpath("haploid_model"),
                patch_vcf,
                hap_vcf,
            )
        )
        second_modelapply_job.add(
            Job(shlex.join(driver.build_cmd()), "second-modelapply", cores)
        )

    # Second pass - unphased regions
    diploid_unphased_hp = tmp_dir.joinpath(
        "out_diploid_phased_unphased_hp.vcf.gz"
    )
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=unphased_bed,
    )
    driver.add_algo(
        DNAscopeHP(
            diploid_unphased_hp,
            dbsnp=dbsnp,
            model=model_bundle.joinpath("diploid_hp_model"),
            pcr_indel_model=repeat_model,
        )
    )
    calling_unphased_job = Job(
        shlex.join(driver.build_cmd()), "calling-unphased", cores
    )

    # Patch DNA and DNAHP variants
    diploid_unphased_patch = tmp_dir.joinpath(
        "out_diploid_unphased_patch.vcf.gz"
    )
    diploid_unphased = tmp_dir.joinpath("out_diploid_unphased.vcf.gz")
    cmd = cmds.cmd_pyexec_vcf_mod_patch(
        str(diploid_unphased_patch),
        str(phased_unphased),
        str(diploid_unphased_hp),
        cores,
        kwargs,
    )
    diploid_patch_job = Job(cmd, "diploid-patch", cores)
    driver = Driver(
        reference=reference,
        thread_count=cores,
    )
    driver.add_algo(
        DNAModelApply(
            model_bundle.joinpath("diploid_model_unphased"),
            diploid_unphased_patch,
            diploid_unphased,
        )
    )
    modelapply_unphased_job = Job(
        shlex.join(driver.build_cmd()), "modelapply-unphased", cores
    )

    # merge calls to create the output
    diploid_merged_vcf = tmp_dir.joinpath("out_diploid_merged.vcf.gz")
    merge_out_vcf = diploid_merged_vcf if haploid_bed else output_vcf
    merge_job = Job(
        cmds.cmd_pyexec_vcf_mod_merge(
            str(hap_vcfs[0]),
            str(hap_vcfs[1]),
            str(diploid_unphased),
            str(phased_vcf),
            str(phased_bed),
            str(merge_out_vcf),
            cores,
            kwargs,
        ),
        "merge",
        cores,
    )

    gvcf_combine_job = None
    if gvcf:
        gvcf_combine_job = Job(
            cmds.cmd_pyexec_gvcf_combine(
                reference,
                str(diploid_gvcf_fn),
                str(merge_out_vcf),
                cores,
                kwargs,
            ),
            "gvcf-combine",
            0,
        )

    haploid_calling_job = None
    haploid_patch2_job = None
    haploid_concat_job = None
    haploid_gvcf_combine_job = None
    haploid_gvcf_concat_job = None
    if haploid_bed:
        # Haploid variant calling
        haploid_fn = tmp_dir.joinpath("haploid.vcf.gz")
        haploid_gvcf_fn = tmp_dir.joinpath("haploid.g.vcf.gz")
        haploid_hp_fn = tmp_dir.joinpath("haploid_hp.vcf.gz")
        haploid_out_fn = tmp_dir.joinpath("haploid_patched.vcf.gz")
        driver = Driver(
            reference=reference,
            thread_count=cores,
            input=sample_input,
            interval=haploid_bed,
        )
        driver.add_algo(
            DNAscope(
                haploid_fn,
                dbsnp=dbsnp,
                model=model_bundle.joinpath("haploid_model"),
            )
        )
        driver.add_algo(
            DNAscopeHP(
                haploid_hp_fn,
                dbsnp=dbsnp,
                model=model_bundle.joinpath("haploid_hp_model"),
                pcr_indel_model=repeat_model,
            )
        )
        if gvcf:
            driver.add_algo(
                DNAscope(
                    haploid_gvcf_fn,
                    dbsnp=dbsnp,
                    emit_mode="gvcf",
                    ploidy=1,
                    model=model_bundle.joinpath("gvcf_model"),
                )
            )
        haploid_calling_job = Job(
            shlex.join(driver.build_cmd()), "haploid-calling", cores
        )

        haploid_patch2_job = Job(
            cmds.cmd_pyexec_vcf_mod_haploid_patch2(
                str(haploid_out_fn),
                str(haploid_fn),
                str(haploid_hp_fn),
                cores,
                kwargs,
            ),
            "haploid-patch2",
            cores,
        )
        haploid_concat_job = Job(
            cmds.bcftools_concat(
                output_vcf,
                [diploid_merged_vcf, haploid_out_fn],
            ),
            "haploid-diploid-concat",
            0,
        )

        if gvcf:
            output_gvcf = pathlib.Path(
                str(output_vcf).replace(".vcf.gz", ".g.vcf.gz")
            )
            diploid_gvcf = pathlib.Path(
                str(diploid_merged_vcf).replace(".vcf.gz", ".g.vcf.gz")
            )
            haploid_gvcf = pathlib.Path(
                str(haploid_out_fn).replace(".vcf.gz", ".g.vcf.gz")
            )
            haploid_gvcf_combine_job = Job(
                cmds.cmd_pyexec_gvcf_combine(
                    reference,
                    str(haploid_gvcf_fn),
                    str(haploid_out_fn),
                    cores,
                    kwargs,
                ),
                "haploid-gvcf-combine",
                0,
            )
            haploid_gvcf_concat_job = Job(
                cmds.bcftools_concat(
                    output_gvcf,
                    [diploid_gvcf, haploid_gvcf],
                ),
                "haploid-gvcf-concat",
                0,
            )
    return (
        first_calling_job,
        first_modelapply_job,
        phaser_job,
        bcftools_subset_phased_job,
        fai_to_bed_job,
        bcftools_subtract_job,
        repeatmodel_job,
        bcftools_subset_unphased_job,
        second_calling_job,
        haploid_patch_job,
        second_modelapply_job,
        calling_unphased_job,
        diploid_patch_job,
        modelapply_unphased_job,
        merge_job,
        gvcf_combine_job,
        haploid_calling_job,
        haploid_patch2_job,
        haploid_concat_job,
        haploid_gvcf_combine_job,
        haploid_gvcf_concat_job,
    )


def call_svs(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    bed: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> Job:
    """
    Call SVs using Sentieon LongReadSV
    """
    if not skip_version_check:
        for cmd, min_version in SV_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    sv_vcf = pathlib.Path(str(output_vcf).replace(".vcf.gz", ".sv.vcf.gz"))
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    driver.add_algo(
        LongReadSV(
            sv_vcf,
            model_bundle.joinpath("longreadsv.model"),
        )
    )
    longreadsv_job = Job(shlex.join(driver.build_cmd()), "LongReadSV", cores)
    return longreadsv_job


def mosdepth(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> Set[Job]:
    """Run mosdepth for QC"""

    if not skip_version_check:
        if not all(
            [
                check_version(cmd, min_version)
                for (cmd, min_version) in MOSDEPTH_MIN_VERSIONS.items()
            ]
        ):
            logger.warning(
                "Skipping mosdepth. mosdepth version %s or later not found",
                MOSDEPTH_MIN_VERSIONS["mosdepth"],
            )
            return set()

    mosdpeth_jobs = set()
    for i, input_file in enumerate(sample_input):
        mosdepth_dir = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_mosdepth_{i}")
        )
        mosdpeth_jobs.add(
            Job(
                cmds.cmd_mosdepth(
                    input_file,
                    mosdepth_dir,
                    fasta=reference,
                    threads=cores,
                ),
                f"mosdepth-{i}",
                0,  # Run in background
            )
        )
    return mosdpeth_jobs


def merge_input_files(
    tmp_dir: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **kwargs: Any,
) -> Tuple[pathlib.Path, Job]:
    """
    Merge the aligned reads into a single BAM
    """
    if not skip_version_check:
        for cmd, min_version in MERGE_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    # Merge the sample_input into a single BAM
    merged_bam = tmp_dir.joinpath("longread_merged.bam")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
    )
    driver.add_algo(ReadWriter(merged_bam))
    merge_job = Job(
        shlex.join(driver.build_cmd()),
        "merge-bam",
        0,
    )
    return (merged_bam, merge_job)


def pbsv(
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    merged_bam: pathlib.Path,
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **kwargs: Any,
) -> Tuple[Job, Job]:
    """
    Call SVs with pbsv
    """
    if not skip_version_check:
        for cmd, min_version in PBSV_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    # pbsv discover
    pbsv_discover_fn = str(output_vcf).replace(
        ".vcf.gz",
        ".pbsv_discover.svsig.gz",
    )
    pbsv_discover = Job(
        f"pbsv discover --hifi {merged_bam} {pbsv_discover_fn}",
        "pbsv-discover",
        0,
    )

    # pbsv call
    pbsv_call_fn = str(output_vcf).replace(".vcf.gz", ".pbsv.vcf")
    pbsv_call = Job(
        f"pbsv call -j {cores} {reference} {pbsv_discover_fn} {pbsv_call_fn}",
        "pbsv-call",
        cores,
    )
    return (pbsv_discover, pbsv_call)


def hificnv(
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    cnv_excluded_regions: Optional[pathlib.Path],
    reference: pathlib.Path,
    merged_bam: pathlib.Path,
    sample_input: List[pathlib.Path],
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **kwargs: Any,
) -> Job:
    """
    Call CNVs with HiFiCNV
    """
    if not skip_version_check:
        for cmd, min_version in HIFICNV_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    hifi_cnv_fn = str(output_vcf).replace(".vcf.gz", ".hificnv")
    cmd = (
        f"hificnv --bam {merged_bam} --ref {reference} --threads {cores} "
        f"--output-prefix {hifi_cnv_fn}"
    )
    if cnv_excluded_regions:
        cmd += f" --exclude {cnv_excluded_regions}"
    hificnv_job = Job(cmd, "hificnv", cores)
    return hificnv_job


@arg(
    "-r",
    "--reference",
    required=True,
    help="fasta for reference genome",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "-i",
    "--sample_input",
    nargs="*",
    help="sample BAM or CRAM file",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--fastq",
    nargs="*",
    help="Sample fastq files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--readgroups",
    nargs="*",
    help="Readgroup information for the fastq files",
)
@arg(
    "-m",
    "--model_bundle",
    help="The model bundle file",
    required=True,
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "output_vcf",
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
    "--haploid_bed",
    help=(
        "A BED file of haploid regions. Supplying this file will perform "
        "haploid variant calling across these regions."
    ),
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--cnv_excluded_regions",
    help="Regions to exclude from CNV calling.",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "-t",
    "--cores",
    help="Number of threads/processes to use. %(default)s",
)
@arg(
    "-g",
    "--gvcf",
    help="Generate a gVCF output file along with the VCF."
    " (default generates only the VCF)",
    action="store_true",
)
@arg(
    "--tech",
    help="Sequencing technology used to generate the reads.",
    choices=["HiFi", "ONT"],
)
@arg(
    "--dry_run",
    help="Print the commands without running them.",
)
@arg(
    "--skip_small_variants",
    help="Skip small variant (SNV/indel) calling",
)
@arg(
    "--skip_svs",
    help="Skip SV calling",
)
@arg(
    "--skip_mosdepth",
    help="Skip QC with mosdepth",
)
@arg(
    "--skip_cnv",
    help="Skip CNV calling",
)
@arg(
    "--align",
    help="Align the input BAM/CRAM/uBAM file to the reference genome",
    action="store_true",
)
@arg(
    "--input_ref",
    help="Used to decode the input alignment file. Required if the input file"
    " is in the CRAM/uCRAM formats",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--fastq_taglist",
    help="A comma-separated list of tags to retain. Defaults to '%(default)s'"
    " and the 'RG' tag is required",
)
@arg(
    "--bam_format",
    help="Use the BAM format instead of CRAM for output aligned files",
    action="store_true",
)
@arg(
    "--minimap2_args",
    help="Extra arguments for sentieon minimap2",
)
@arg(
    "--util_sort_args",
    help="Extra arguments for sentieon util sort",
)
@arg(
    "--repeat_model",
    type=path_arg(exists=True, is_file=True),
    help=argparse.SUPPRESS,
)
@arg(
    "--skip_version_check",
    help=argparse.SUPPRESS,
    action="store_true",
)
@arg(
    "--retain_tmpdir",
    help=argparse.SUPPRESS,
    action="store_true",
)
@arg(
    "--use_pbsv",
    help=argparse.SUPPRESS,
    action="store_true",
)
def dnascope_longread(
    output_vcf: pathlib.Path,
    reference: Optional[pathlib.Path] = None,
    sample_input: Optional[List[pathlib.Path]] = None,
    fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    model_bundle: Optional[pathlib.Path] = None,
    dbsnp: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    bed: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    haploid_bed: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    cnv_excluded_regions: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),  # pylint: disable=W0613
    gvcf: bool = False,  # pylint: disable=W0613
    tech: str = "HiFi",
    dry_run: bool = False,
    skip_small_variants: bool = False,
    skip_svs: bool = False,
    skip_mosdepth: bool = False,
    skip_cnv: bool = False,
    align: bool = False,
    input_ref: Optional[pathlib.Path] = None,
    fastq_taglist: str = "*",  # pylint: disable=W0613
    bam_format: bool = False,  # pylint: disable=W0613
    minimap2_args: str = "-Y",  # pylint: disable=W0613
    util_sort_args: str = (
        "--cram_write_options version=3.0,compressor=rans"
    ),  # pylint: disable=W0613
    repeat_model: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    skip_version_check: bool = False,  # pylint: disable=W0613
    retain_tmpdir: bool = False,
    use_pbsv: bool = False,
    **kwargs: str,
):
    """
    Run sentieon cli with the algo DNAscope command.
    """
    assert reference
    assert sample_input or fastq
    assert model_bundle
    assert str(output_vcf).endswith(".vcf.gz")

    assert logger.parent
    logger.parent.setLevel(kwargs["loglevel"])
    logger.info("Starting sentieon-cli version: %s", __version__)

    if tech.upper() == "ONT":
        skip_cnv = True

    if not library_preloaded("libjemalloc.so"):
        logger.warning(
            "jemalloc is recommended, but is not preloaded. See "
            "https://support.sentieon.com/appnotes/jemalloc/"
        )

    if not cnv_excluded_regions and not skip_cnv:
        logger.warning(
            "Excluded regions are recommended for CNV calling. Please supply "
            "them with the `--cnv_excluded_regions` argument"
        )

    tmp_dir_str = tmp()
    tmp_dir = pathlib.Path(tmp_dir_str)  # type: ignore  # pylint: disable=W0641  # noqa: E501

    logger.info("Building the DAG")
    dag = DAG()

    sample_input = sample_input if sample_input else []
    realign_jobs: Set[Job] = set()
    if align:
        sample_input, realign_jobs = align_inputs(**locals())
        for job in realign_jobs:
            dag.add_job(job)
    aligned_fastq, align_jobs = align_fastq(**locals())
    sample_input.extend(aligned_fastq)
    for job in align_jobs:
        dag.add_job(job)

    if not skip_mosdepth:
        mosdpeth_jobs = mosdepth(**locals())
        for job in mosdpeth_jobs:
            dag.add_job(job, realign_jobs.union(align_jobs))

    if use_pbsv or not skip_cnv:
        merged_bam, merge_job = merge_input_files(**locals())
        dag.add_job(merge_job, realign_jobs.union(align_jobs))

        if use_pbsv:
            pbsv_discover, pbsv_call = pbsv(**locals())
            dag.add_job(pbsv_discover, {merge_job})
            dag.add_job(pbsv_call, {pbsv_discover})

        if not skip_cnv:
            hificnv_job = hificnv(**locals())
            dag.add_job(hificnv_job, {merge_job})

    if not skip_small_variants:
        (
            first_calling_job,
            first_modelapply_job,
            phaser_job,
            bcftools_subset_phased_job,
            fai_to_bed_job,
            bcftools_subtract_job,
            repeatmodel_job,
            bcftools_subset_unphased_job,
            second_calling_job,
            haploid_patch_job,
            second_modelapply_job,
            calling_unphased_job,
            diploid_patch_job,
            modelapply_unphased_job,
            merge_job,
            gvcf_combine_job,
            haploid_calling_job,
            haploid_patch2_job,
            haploid_concat_job,
            haploid_gvcf_combine_job,
            haploid_gvcf_concat_job,
        ) = call_variants(**locals())
        dag.add_job(first_calling_job, realign_jobs.union(align_jobs))
        dag.add_job(first_modelapply_job, {first_calling_job})
        dag.add_job(phaser_job, {first_modelapply_job})

        haploid_patch_deps = set()
        if bcftools_subset_phased_job:
            dag.add_job(bcftools_subset_phased_job, {phaser_job})
            haploid_patch_deps.add(bcftools_subset_phased_job)

        subtract_deps = {phaser_job}
        if fai_to_bed_job:
            dag.add_job(fai_to_bed_job)
            subtract_deps.add(fai_to_bed_job)
        dag.add_job(bcftools_subtract_job, subtract_deps)
        dag.add_job(bcftools_subset_unphased_job, {bcftools_subtract_job})

        second_pass_deps = {phaser_job}
        calling_unphased_deps = {bcftools_subtract_job}
        haploid_calling_deps = set()
        if repeatmodel_job:
            dag.add_job(repeatmodel_job, {phaser_job})
            second_pass_deps.add(repeatmodel_job)
            calling_unphased_deps.add(repeatmodel_job)
            haploid_calling_deps.add(repeatmodel_job)

        merge_deps = set()
        for job in second_calling_job:
            dag.add_job(job, second_pass_deps)
            haploid_patch_deps.add(job)
        dag.add_job(haploid_patch_job, haploid_patch_deps)
        for job in second_modelapply_job:
            dag.add_job(job, {haploid_patch_job})
            merge_deps.add(job)

        dag.add_job(calling_unphased_job, calling_unphased_deps)
        dag.add_job(
            diploid_patch_job,
            {bcftools_subset_unphased_job, calling_unphased_job},
        )
        dag.add_job(modelapply_unphased_job, {diploid_patch_job})
        merge_deps.add(modelapply_unphased_job)
        dag.add_job(merge_job, merge_deps)

        if gvcf_combine_job:  # if gvcf
            dag.add_job(gvcf_combine_job, {merge_job})

        if haploid_calling_job and haploid_patch2_job and haploid_concat_job:
            # if haploid bed
            dag.add_job(haploid_calling_job, haploid_calling_deps)
            dag.add_job(haploid_patch2_job, {haploid_calling_job})
            dag.add_job(haploid_concat_job, {haploid_patch2_job, merge_job})

            if (
                haploid_gvcf_combine_job
                and haploid_gvcf_concat_job
                and gvcf_combine_job
            ):
                # if gvcf
                dag.add_job(haploid_gvcf_combine_job, {haploid_patch2_job})
                dag.add_job(
                    haploid_gvcf_concat_job,
                    {haploid_gvcf_combine_job, gvcf_combine_job},
                )

    if not skip_svs:
        longreadsv_job = call_svs(**locals())
        dag.add_job(longreadsv_job, realign_jobs.union(align_jobs))

    logger.debug("Creating the scheduler")
    scheduler = ThreadScheduler(
        dag,
        cores,
    )

    logger.debug("Creating the executor")
    Executor = DryRunExecutor if dry_run else LocalExecutor
    executor = Executor(scheduler)

    logger.info("Starting execution")
    executor.execute()

    if not retain_tmpdir:
        shutil.rmtree(tmp_dir_str)

    if executor.jobs_with_errors:
        raise CommandError("Execution failed")

    if len(dag.waiting_jobs) > 0 or len(dag.ready_jobs) > 0:
        raise CommandError(
            "The DAG has some unexecuted jobs\n"
            f"Waiting jobs: {dag.waiting_jobs}\n"
            f"Ready jobs: {dag.ready_jobs}\n"
        )
    return
