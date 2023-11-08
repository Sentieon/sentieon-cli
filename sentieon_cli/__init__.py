import multiprocessing as mp
import argparse
import os
import sys
import subprocess as sp
import pathlib
import shutil
import tempfile
from typing import Any, Callable, Optional, List

import argh
import packaging.version

from argh import arg
from importlib_resources import files
from .logging import get_logger
from . import command_strings as cmds

__version__ = "0.1.0"

logger = get_logger(__name__)

TOOL_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202308"),
    "bcftools": packaging.version.Version("1.10"),
    "bedtools": None,
}


def tmp():
    """Create a temporary directory for the current process."""
    tmp_base = os.getenv("SENTIEON_TMPDIR")
    tmp_dir = tempfile.TemporaryDirectory(dir=tmp_base)
    return tmp_dir


def check_version(cmd: str, version: Optional[packaging.version.Version]):
    """Check the version of an executable"""
    cmd_list: List[str] = cmd.split()
    exec_file = shutil.which(cmd_list[0])
    if not exec_file:
        print(f"Error: no '{cmd}' found in the PATH")
        sys.exit(2)

    if version is None:
        return

    cmd_list.append("--version")
    cmd_version_str = sp.check_output(cmd_list).decode("utf-8").strip()
    if cmd_list[0] == "sentieon":
        cmd_version_str = cmd_version_str.split("-")[-1]
    else:
        # handle, e.g. bcftools which outputs multiple lines.
        cmd_version_str = (
            cmd_version_str.split("\n")[0].split()[-1].split("-")[0]
        )
    cmd_version = packaging.version.Version(cmd_version_str)
    if cmd_version < version:
        print(
            f"Error: the pipeline requires {cmd} version '{version}' or later "
            f"but {cmd} '{cmd_version}' was found in the PATH"
        )
        sys.exit(2)
    return


def path_arg(
    exists: Optional[bool] = None,
    is_dir: Optional[bool] = None,
    is_file: Optional[bool] = None,
    is_fifo: Optional[bool] = None,
) -> Callable[[str], pathlib.Path]:
    """pathlib checked types for argparse"""

    def _path_arg(arg: str) -> pathlib.Path:
        p = pathlib.Path(arg)

        attrs = [exists, is_dir, is_file, is_fifo]
        attr_names = ["exists", "is_dir", "is_file", "is_fifo"]

        for attr_val, attr_name in zip(attrs, attr_names):
            if attr_val is None:  # Skip attributes that are not defined
                continue

            m = getattr(p, attr_name)
            if m() != attr_val:
                raise argparse.ArgumentTypeError(
                    "The supplied path argument needs the attribute"
                    f" {attr_name}={attr_val}, but {attr_name}={m()}"
                )
        return p

    return _path_arg


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
    default=mp.cpu_count(),
    help="Number of threads/processes to use. %(default)s",
)
@arg(
    "-g",
    "--gvcf",
    help="Generate a gVCF output file along with the VCF."
    " (default generates only the VCF)",
    default=False,
    action="store_true",
)
@arg(
    "--tech",
    help="Sequencing technology used to generate the reads.",
    default="HiFi",
    choices=["HiFi", "ONT"],
)
@arg(
    "--dry-run",
    help="Print the commands without running them.",
    default=False,
    action="store_true",
)
def dnascope_longread(**kwargs: Any):
    """
    Run sentieon cli with the algo DNAscope command.
    """
    logger.info(kwargs)

    for cmd, min_version in TOOL_MIN_VERSIONS.items():
        check_version(cmd, min_version)

    tmp_dir_obj = tmp()
    tmp_dir = pathlib.Path(tmp_dir_obj.name)

    reference: pathlib.Path = kwargs["reference"]
    sample_input: List[pathlib.Path] = kwargs["sample_input"]
    model_bundle: pathlib.Path = kwargs["model_bundle"]
    output_vcf: pathlib.Path = kwargs["output-vcf"]
    dbsnp: Optional[pathlib.Path] = kwargs["dbsnp"]
    bed: Optional[pathlib.Path] = kwargs["bed"]
    cores: int = kwargs["cores"]
    use_gvcf: bool = kwargs["gvcf"]
    tech: str = kwargs["tech"]
    dry_run: bool = kwargs["dry_run"]

    if dry_run:
        run = print
    else:
        from .runner import run  # type: ignore[assignment]

    # First pass - diploid calling
    diploid_gvcf_fn = tmp_dir.joinpath("out_diploid.g.vcf.gz")
    diploid_tmp_vcf = tmp_dir.joinpath("out_diploid_tmp.vcf.gz")
    driver = cmds.Driver(
        reference=reference,
        thread_count=cores,
        input=sample_input,
        interval=bed,
    )
    if use_gvcf:
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
    run(" ".join(driver.build_cmd()))

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
    run(" ".join(driver.build_cmd()))

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
    run(" ".join(driver.build_cmd()))

    if tech == "ONT":
        run(
            f"bcftools view -T {phased_bed} {phased_vcf} \
            | sentieon util vcfconvert - {phased_phased}"
        )
    run(
        cmds.cmd_bedtools_subtract(
            bed, phased_bed, unphased_bed, tmp_dir, **kwargs
        )
    )

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
    run(" ".join(driver.build_cmd()))

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

        if tech == "HiFi":
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
        run(" ".join(driver.build_cmd()))

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
            kwargs["tech"],
            str(phased_phased),
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
        run(" ".join(driver.build_cmd()))

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
    run(" ".join(driver.build_cmd()))

    # Patch DNA and DNAHP variants
    diploid_unphased_patch = tmp_dir.joinpath(
        "out_diploid_unphased_patch.vcf.gz"
    )
    diploid_unphased = tmp_dir.joinpath("out_diploid_unphased.vcf.gz")
    cmd = cmds.cmd_pyexec_vcf_mod_patch(
        str(diploid_unphased_patch),
        str(phased_unphased),
        str(diploid_unphased_hp),
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
    run(" ".join(driver.build_cmd()))

    # merge calls to create the output
    run(
        cmds.cmd_pyexec_vcf_mod_merge(
            str(hap_vcfs[0]),
            str(hap_vcfs[1]),
            str(diploid_unphased),
            str(phased_vcf),
            str(phased_bed),
            str(output_vcf),
            kwargs,
        )
    )

    if use_gvcf:
        run(
            cmds.cmd_pyexec_gvcf_combine(
                str(diploid_gvcf_fn),
                str(output_vcf),
                kwargs,
            )
        )

    shutil.rmtree(tmp_dir)
    return


def main():
    """main entry point for this project"""
    logger.setLevel(os.environ.get("LOGLEVEL", "DEBUG").upper())
    logger.info("Starting sentieon-cli version: %s", __version__)

    argh.dispatch_commands([dnascope_longread])


if __name__ == "__main__":
    main()
