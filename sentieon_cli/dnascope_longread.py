"""
Functionality for the DNAscope LongRead pipeline
"""

import argparse
import multiprocessing as mp
import pathlib
import shlex
import shutil
from typing import Any, Callable, Dict, List, Optional

import packaging.version

from argh import arg
from importlib_resources import files

from . import command_strings as cmds
from .util import __version__, check_version, logger, path_arg, tmp

TOOL_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
    "bcftools": packaging.version.Version("1.10"),
    "bedtools": None,
}

SV_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
}


def call_variants(
    run: Callable[[str], None],
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sample_input: List[pathlib.Path],
    model_bundle: pathlib.Path,
    dbsnp: Optional[pathlib.Path] = None,
    bed: Optional[pathlib.Path] = None,
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
    if not skip_version_check:
        for cmd, min_version in TOOL_MIN_VERSIONS.items():
            check_version(cmd, min_version)

    # First pass - diploid calling
    diploid_gvcf_fn = tmp_dir.joinpath("out_diploid.g.vcf.gz")
    diploid_tmp_vcf = tmp_dir.joinpath("out_diploid_tmp.vcf.gz")
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    if gvcf:
        driver.add_algo(
            cmds.DNAscope(
                diploid_gvcf_fn,
                dbsnp=dbsnp,
                emit_mode="gvcf",
                model=model_bundle.joinpath("gvcf_model"),
            )
        )
    driver.add_algo(
        cmds.DNAscope(
            diploid_tmp_vcf,
            dbsnp=dbsnp,
            model=model_bundle.joinpath("diploid_model"),
        )
    )
    run(shlex.join(driver.build_cmd()))

    diploid_vcf = tmp_dir.joinpath("out_diploid.vcf.gz")
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
    )
    driver.add_algo(
        cmds.DNAModelApply(
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
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    driver.add_algo(
        cmds.VariantPhaser(
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
        driver = cmds.Driver(
            reference=reference,
            thread_count=cores,
            input=sample_input,
            interval=phased_bed,
            read_filter=f"PhasedReadFilter,phased_vcf={phased_ext},phase_select=tag",  # noqa: E501
        )
        driver.add_algo(
            cmds.RepeatModel(
                repeat_model,
                phased=True,
                read_flag_mask="drop=supplementary",
            )
        )
        shlex.join(driver.build_cmd())

    run(
        f"bcftools view -T {unphased_bed} {phased_vcf} \
        | sentieon util vcfconvert - {phased_unphased}"
    )

    # Second pass - phased variants
    for phase in (1, 2):
        hp_std_vcf = tmp_dir.joinpath(f"out_hap{phase}_nohp_tmp.vcf.gz")
        hp_vcf = tmp_dir.joinpath(f"out_hap{phase}_tmp.vcf.gz")
        driver = cmds.Driver(
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
                cmds.DNAscope(
                    hp_std_vcf,
                    dbsnp=dbsnp,
                    model=model_bundle.joinpath("haploid_model"),
                )
            )
        driver.add_algo(
            cmds.DNAscopeHP(
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
        driver = cmds.Driver(
            reference=reference,
            thread_count=cores,
        )
        driver.add_algo(
            cmds.DNAModelApply(
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
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=unphased_bed,
    )
    driver.add_algo(
        cmds.DNAscopeHP(
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
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
    )
    driver.add_algo(
        cmds.DNAModelApply(
            model_bundle.joinpath("diploid_model_unphased"),
            diploid_unphased_patch,
            diploid_unphased,
        )
    )
    run(shlex.join(driver.build_cmd()))

    # merge calls to create the output
    run(
        cmds.cmd_pyexec_vcf_mod_merge(
            str(hap_vcfs[0]),
            str(hap_vcfs[1]),
            str(diploid_unphased),
            str(phased_vcf),
            str(phased_bed),
            str(output_vcf),
            cores,
            kwargs,
        )
    )

    if gvcf:
        run(
            cmds.cmd_pyexec_gvcf_combine(
                reference,
                str(diploid_gvcf_fn),
                str(output_vcf),
                cores,
                kwargs,
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
            check_version(cmd, min_version)

    sv_vcf = pathlib.Path(str(output_vcf).replace(".vcf.gz", ".sv.vcf.gz"))
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    driver.add_algo(
        cmds.LongReadSV(
            sv_vcf,
            model_bundle.joinpath("longreadsv.model"),
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
    required=True,
    nargs="+",
    help="sample BAM or CRAM file",
    type=path_arg(exists=True, is_file=True),
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
    "--repeat-model",
    type=path_arg(exists=True, is_file=True),
    help=argparse.SUPPRESS,
)
@arg(
    "--skip-version-check",
    help=argparse.SUPPRESS,
)
@arg(
    "--retain-tmpdir",
    help=argparse.SUPPRESS,
)
def dnascope_longread(
    output_vcf: pathlib.Path,  # pylint: disable=W0613
    reference: Optional[pathlib.Path] = None,
    sample_input: Optional[List[pathlib.Path]] = None,
    model_bundle: Optional[pathlib.Path] = None,
    dbsnp: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    bed: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    cores: int = mp.cpu_count(),  # pylint: disable=W0613
    gvcf: bool = False,  # pylint: disable=W0613
    tech: str = "HiFi",  # pylint: disable=W0613
    dry_run: bool = False,
    skip_small_variants: bool = False,
    skip_svs: bool = False,
    repeat_model: Optional[pathlib.Path] = None,  # pylint: disable=W0613
    skip_version_check: bool = False,  # pylint: disable=W0613
    retain_tmpdir: bool = False,
    **kwargs: str,
):
    """
    Run sentieon cli with the algo DNAscope command.
    """
    assert reference
    assert sample_input
    assert model_bundle

    logger.setLevel(kwargs["loglevel"])
    logger.info("Starting sentieon-cli version: %s", __version__)

    tmp_dir_str = tmp()
    tmp_dir = pathlib.Path(tmp_dir_str)  # type: ignore  # pylint: disable=W0641  # noqa: E501

    if dry_run:
        run = print  # type: ignore  # pylint: disable=W0641
    else:
        from .runner import run  # type: ignore[assignment]  # noqa: F401

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
