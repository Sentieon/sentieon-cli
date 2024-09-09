"""
Functionality for the DNAscope LongRead pipeline
"""

import argparse
import multiprocessing as mp
import pathlib
import shlex
import shutil
import sys
from typing import Any, Callable, Dict, List, Optional

import packaging.version

from argh import arg
from importlib_resources import files

from . import command_strings as cmds
from .driver import (
    Driver,
    DNAscope,
    DNAscopeHP,
    DNAModelApply,
    LongReadSV,
    RepeatModel,
    VariantPhaser,
)
from .logging import get_logger
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


def align_inputs(
    run: Callable[[str], None],
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    cores: int = mp.cpu_count(),
    dry_run: bool = False,
    skip_version_check: bool = False,
    bam_format: bool = False,
    fastq_taglist: str = "*",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    input_ref: Optional[pathlib.Path] = None,
    **_kwargs: Any,
) -> List[pathlib.Path]:
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
    for i, input_aln in enumerate(sample_input):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_mm2_sorted_{i}.{suffix}")
        )
        rg_lines = cmds.get_rg_lines(
            input_aln,
            dry_run,
        )

        run(
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
                util_sort_args,
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
    fastq: Optional[List[pathlib.Path]] = None,
    readgroups: Optional[List[str]] = None,
    skip_version_check: bool = False,
    bam_format: bool = False,
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    **_kwargs: Any,
) -> List[pathlib.Path]:
    """
    Align fastq to the reference genome using minimap2
    """
    res: List[pathlib.Path] = []
    if fastq is None and readgroups is None:
        return res
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
    for i, (fq, rg) in enumerate(zip(fastq, readgroups)):
        out_aln = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_mm2_sorted_fq_{i}.{suffix}")
        )
        run(
            cmds.cmd_fastq_minimap2(
                out_aln,
                fq,
                rg,
                reference,
                model_bundle,
                cores,
                unzip,
                util_sort_args,
            )
        )
        res.append(out_aln)
    return res


def call_variants(
    run: Callable[[str], None],
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
) -> int:
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
    run(shlex.join(driver.build_cmd()))

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
    run(shlex.join(driver.build_cmd()))

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
    run(shlex.join(driver.build_cmd()))

    if tech.upper() == "ONT":
        run(
            f"bcftools view -T {phased_bed} {phased_vcf} \
            | sentieon util vcfconvert - {phased_phased}"
        )
    run(
        cmds.cmd_bedtools_subtract(
            bed, phased_bed, unphased_bed, tmp_dir, reference, dry_run
        )
    )

    if not repeat_model:
        repeat_model = tmp_dir.joinpath("out_repeat.model")
        driver = Driver(
            reference=reference,
            thread_count=cores,
            input=sample_input,
            interval=phased_bed,
            read_filter=f"PhasedReadFilter,phased_vcf={phased_ext},phase_select=tag",  # noqa: E501
        )
        driver.add_algo(
            RepeatModel(
                repeat_model,
                phased=True,
                read_flag_mask="drop=supplementary",
            )
        )
        run(shlex.join(driver.build_cmd()))

    run(
        f"bcftools view -T {unphased_bed} {phased_vcf} \
        | sentieon util vcfconvert - {phased_unphased}"
    )

    # Second pass - phased variants
    for phase in (1, 2):
        hp_std_vcf = tmp_dir.joinpath(f"out_hap{phase}_nohp_tmp.vcf.gz")
        hp_vcf = tmp_dir.joinpath(f"out_hap{phase}_tmp.vcf.gz")
        driver = Driver(
            reference=reference,
            thread_count=cores,
            input=sample_input,
            interval=phased_bed,
            read_filter=(
                f"PhasedReadFilter,phased_vcf={phased_ext}"
                f",phase_select={phase}"
            ),
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
        run(shlex.join(driver.build_cmd()))

    kwargs: Dict[str, str] = dict()
    kwargs["gvcf_combine_py"] = str(
        files("sentieon_cli.scripts").joinpath("gvcf_combine.py")
    )
    kwargs["vcf_mod_py"] = str(
        files("sentieon_cli.scripts").joinpath("vcf_mod.py")
    )

    patch_vcfs = [tmp_dir.joinpath(f"out_hap{i}_patch.vcf.gz") for i in (1, 2)]
    run(
        cmds.cmd_pyexec_vcf_mod_haploid_patch(
            str(patch_vcfs[0]),
            str(patch_vcfs[1]),
            f"{tmp_dir}/out_hap%d_%stmp.vcf.gz",
            tech,
            str(phased_phased),
            cores,
            kwargs,
        )
    )

    # apply trained model to the patched vcfs.
    hap_vcfs = [tmp_dir.joinpath(f"out_hap{i}.vcf.gz") for i in (1, 2)]
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
        run(shlex.join(driver.build_cmd()))

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
    run(shlex.join(driver.build_cmd()))

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
    run(cmd)
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
    run(shlex.join(driver.build_cmd()))

    # merge calls to create the output
    diploid_merged_vcf = tmp_dir.joinpath("out_diploid_merged.vcf.gz")
    merge_out_vcf = diploid_merged_vcf if haploid_bed else output_vcf
    run(
        cmds.cmd_pyexec_vcf_mod_merge(
            str(hap_vcfs[0]),
            str(hap_vcfs[1]),
            str(diploid_unphased),
            str(phased_vcf),
            str(phased_bed),
            str(merge_out_vcf),
            cores,
            kwargs,
        )
    )

    if gvcf:
        run(
            cmds.cmd_pyexec_gvcf_combine(
                reference,
                str(diploid_gvcf_fn),
                str(merge_out_vcf),
                cores,
                kwargs,
            )
        )

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
        run(shlex.join(driver.build_cmd()))

        run(
            cmds.cmd_pyexec_vcf_mod_haploid_patch2(
                str(haploid_out_fn),
                str(haploid_fn),
                str(haploid_hp_fn),
                cores,
                kwargs,
            )
        )
        run(
            cmds.bcftools_concat(
                output_vcf,
                [diploid_merged_vcf, haploid_out_fn],
            )
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
            run(
                cmds.cmd_pyexec_gvcf_combine(
                    reference,
                    str(haploid_gvcf_fn),
                    str(haploid_out_fn),
                    cores,
                    kwargs,
                )
            )
            run(
                cmds.bcftools_concat(
                    output_gvcf,
                    [diploid_gvcf, haploid_gvcf],
                )
            )
    return 0


def call_svs(
    run: Callable[[str], None],
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    bed: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> int:
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
    run(shlex.join(driver.build_cmd()))
    return 0


def mosdepth(
    run: Callable[[str], None],
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    cores: int = mp.cpu_count(),
    skip_version_check: bool = False,
    **_kwargs: Any,
) -> int:
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
            return 1

    for i, input_file in enumerate(sample_input):
        mosdepth_dir = pathlib.Path(
            str(output_vcf).replace(".vcf.gz", f"_mosdepth_{i}")
        )
        run(
            cmds.cmd_mosdepth(
                input_file,
                mosdepth_dir,
                fasta=reference,
                threads=cores,
            )
        )
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
    "--haploid-bed",
    help=(
        "A BED file of haploid regions. Supplying this file will perform "
        "haploid variant calling across these regions."
    ),
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
    "--skip-mosdepth",
    help="Skip QC with mosdepth",
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
    "--util_sort_args",
    help="Extra arguments for sentieon util sort",
)
@arg(
    "--repeat-model",
    type=path_arg(exists=True, is_file=True),
    help=argparse.SUPPRESS,
)
@arg(
    "--skip-version-check",
    help=argparse.SUPPRESS,
    action="store_true",
)
@arg(
    "--retain-tmpdir",
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
    cores: int = mp.cpu_count(),  # pylint: disable=W0613
    gvcf: bool = False,  # pylint: disable=W0613
    tech: str = "HiFi",  # pylint: disable=W0613
    dry_run: bool = False,
    skip_small_variants: bool = False,
    skip_svs: bool = False,
    skip_mosdepth: bool = False,
    align: bool = False,
    input_ref: Optional[pathlib.Path] = None,
    fastq_taglist: str = "*",  # pylint: disable=W0613
    bam_format: bool = False,  # pylint: disable=W0613
    util_sort_args: str = (
        "--cram_write_options version=3.0,compressor=rans"
    ),  # pylint: disable=W0613
    repeat_model: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    skip_version_check: bool = False,  # pylint: disable=W0613
    retain_tmpdir: bool = False,
    **kwargs: str,
):
    """
    Run sentieon cli with the algo DNAscope command.
    """
    assert reference
    assert sample_input or fastq
    assert model_bundle
    assert str(output_vcf).endswith(".vcf.gz")

    logger.parent.setLevel(kwargs["loglevel"])
    logger.info("Starting sentieon-cli version: %s", __version__)

    if not library_preloaded("libjemalloc.so"):
        logger.warning(
            "jemalloc is recommended, but is not preloaded. See "
            "https://support.sentieon.com/appnotes/jemalloc/"
        )

    tmp_dir_str = tmp()
    tmp_dir = pathlib.Path(tmp_dir_str)  # type: ignore  # pylint: disable=W0641  # noqa: E501

    if dry_run:
        run = print  # type: ignore  # pylint: disable=W0641
    else:
        from .runner import run  # type: ignore[assignment]  # noqa: F401

    sample_input = sample_input if sample_input else []
    if align:
        sample_input = align_inputs(**locals())
    sample_input.extend(align_fastq(**locals()))

    if not skip_mosdepth:
        _res = mosdepth(**locals())

    if not skip_small_variants:
        res = call_variants(**locals())
        if res != 0:
            logger.error("Small variant calling failed")
            return

    if not skip_svs:
        res = call_svs(**locals())
        if res != 0:
            logger.error("SV calling failed")
            return

    if not retain_tmpdir:
        shutil.rmtree(tmp_dir_str)
    return
