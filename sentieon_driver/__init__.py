import multiprocessing as mp
import argparse
import os
import sys
import subprocess as sp
import argh
import random
from argh import arg
from .logging import get_logger
from . import command_strings as cmds

if sys.version_info < (3, 9):
    from typing import Tuple as tuple

__version__ = "0.1.0"

logger = get_logger(__name__)


def tmp() -> tuple[str, str]:
    """Create a temporary directory for the current process."""
    tmp_base = os.getenv("SENTIEON_TMPDIR", os.getenv("TMPDIR", "/tmp"))
    tmp_base = os.path.join(tmp_base, str(os.getpid()))
    tmp_dir = tmp_base
    while os.path.exists(tmp_dir):
        tmp_dir = tmp_base + "_" + str(random.randint(0, 1000000))
    os.makedirs(tmp_dir)
    return tmp_dir, tmp_base


def check_version(cmd: str, version: str):
    cmd_version = sp.check_output([cmd, "--version"]).decode("utf-8").strip()
    # TODO: use regexp
    # https://github.com/Sentieon/sentieon-scripts/blob/master/dnascope_LongRead/dnascope_HiFi.sh#L146
    assert cmd_version


def cmd_algo_variantphaser():
    # https://github.com/Sentieon/sentieon-scripts/blob/8d33f29e442d5a1e782445f06bc1f11e278d8f87/dnascope_LongRead/dnascope_HiFi.sh#L364-L367
    pass


def common(fn):
    fn = arg("path", help="some path")(fn)
    fn = arg(
        "--extra", help="extra integer arguments: %(default)s", default=0
    )(fn)
    return fn


@common
def other(**kwargs):
    print("OTHER")


@arg(
    "-r",
    "--reference",
    help="fasta for reference genome",
    type=argparse.FileType("r"),
)
@arg(
    "-i",
    "--sample-input",
    help="sample BAM or CRAM file",
    type=argparse.FileType("r"),
)
@arg("-m", "--model-bundle", help="The model bundle directory")
@arg(
    "output-vcf",
    help="Output VCF File. The file name must end in .vcf.gz",
    type=argparse.FileType("w"),
)
@arg(
    "-d",
    "--dbsnp",
    help="dbSNP vcf file Supplying this file will annotate variants with \
         their dbSNP refSNP ID numbers.",
    default=None,
    type=argparse.FileType("r"),
)
@arg(
    "-b",
    "--bed",
    help="Region BED file. Supplying this file will limit variant calling \
    to the intervals inside the BED file.",
    type=argparse.FileType("r"),
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
def run_full_dnascope(**kwargs):
    """
    Run sentieon driver with the algo DNAscope command.
    """
    logger.info(kwargs)
    kwargs["tmp_base"], kwargs["tmp_dir"] = tmp()
    model = f"{kwargs['model_bundle']}/gvcf_model"
    commands = [cmds.cmd_algo_dnascope(model, **kwargs)]

    inp_vcf = f"{kwargs['tmp_base']}/out_diploid_tmp.vcf.gz"
    out_vcf = f"{kwargs['tmp_base']}/out_diploid.vcf.gz"
    commands.append(
        cmds.cmd_model_apply(
            f"{kwargs['model_bundle']}/diploid_model", inp_vcf, out_vcf, kwargs
        )
    )

    kwargs["inp_vcf"] = out_vcf
    commands.append(cmds.cmd_variant_phaser(kwargs))

    commands.append(
        cmds.cmd_bedtools_subtract(
            kwargs.get("bed"), kwargs["phased_bed"], kwargs["unphased_bed"]
        )
    )

    commands.append(cmds.cmd_repeat_model(kwargs))

    # TODO: bcftools view
    # https://github.com/Sentieon/sentieon-scripts/blob/master/dnascope_LongRead/dnascope_HiFi.sh#L384

    for phase in (1, 2):
        kwargs[
            "read-filter"
        ] = f"PhasedReadFilter,phased_vcf={kwargs['phased_ext']},phase_select={phase}"
        cmd = cmds.cmd_algo_dnascope(
            f"{kwargs['model_bundle']}/haploid_model",
            bed_key="unphased_bed",
            **kwargs,
        )
        cmd = ""
        commands.append("#PHASE %s" % phase)

        hp_vcf = f"{kwargs['tmp_base']}/out_hap{phase}_tmp.vcf.gz"
        hp_std_vcf = f"{kwargs['tmp_base']}/out_hap{phase}_nohp_tmp.vcf.gz"

        cmd += " " + cmds.cmd_dnascope_hp(
            f"{kwargs['model_bundle']}/haploid_hp_model",
            kwargs["repeat_model"],
            hp_vcf,
            hp_std_vcf,
            kwargs,
        )
        commands.append(cmd)

    return "\n".join(commands)


def main():
    """main entry point for this project"""
    logger.setLevel(os.environ.get("LOGLEVEL", "DEBUG").upper())
    logger.info("Starting sentieon_driver version: %s", __version__)
    argh.dispatch_commands([run_full_dnascope, other])


if __name__ == "__main__":
    main()