import multiprocessing as mp
import argparse
import os
import sys
import subprocess as sp
import argh
import random
import tempfile

from argh import arg
from typing import Any
from .logging import get_logger
from . import command_strings as cmds

if sys.version_info < (3, 9):
    from typing import Tuple as tuple

__version__ = "0.1.0"

logger = get_logger(__name__)


def tmp() -> tempfile.TemporaryDirectory[str]:
    """Create a temporary directory for the current process."""
    tmp_base = os.getenv("SENTIEON_TMPDIR")
    tmp_dir = tempfile.TemporaryDirectory(dir=tmp_base)
    return tmp_dir


def check_version(cmd: str, version: str):
    cmd_version = sp.check_output([cmd, "--version"]).decode("utf-8").strip()
    # TODO: use regexp
    # https://github.com/Sentieon/sentieon-scripts/blob/master/dnascope_LongRead/dnascope_HiFi.sh#L146
    assert cmd_version


@arg(
    "-r",
    "--reference",
    required=True,
    help="fasta for reference genome",
    type=argparse.FileType("r"),
)
@arg(
    "-i",
    "--sample-input",
    required=True,
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
@arg(
    "--tech",
    help="Sequencing technology used to generate the reads.",
    default="HiFi",
    choices=["HiFi", "ONT"],
)
def run_dnascope_longread(**kwargs: dict[str, Any]):
    """
    Run sentieon cli with the algo DNAscope command.
    """
    logger.info(kwargs)
    tmp_dir = tmp()
    gvcf = kwargs.pop("gvcf")
    if gvcf:
        gvcf = f"{tmp_dir.name}/out_diploid.g.vcf.gz"
    model = f"{kwargs['model_bundle']}/diploid_model"
    out_vcf = f"{tmp_dir.name}/out_diploid.vcf.tmp.gz"
    commands = [
        cmds.cmd_algo_dnascope(
            model, out_vcf, kwargs, bed_key="bed", gvcf=gvcf
        )
    ]

    inp_vcf = out_vcf
    out_vcf = f"{tmp_dir.name}/out_diploid.vcf.gz"
    commands.append(
        cmds.cmd_model_apply(
            f"{kwargs['model_bundle']}/diploid_model", inp_vcf, out_vcf, kwargs
        )
    )

    phased_bed = f"{tmp_dir.name}/out_diploid_phased.bed"
    unphased_bed = f"{tmp_dir.name}/out_diploid_unphased.bed"
    phased_vcf = f"{tmp_dir.name}/out_diploid_phased.vcf.gz"
    phased_ext = f"{tmp_dir.name}/out_diploid_phased.ext.vcf.gz"
    phased_unphased = (
        f"{tmp_dir.name}/out_diploid_phased_unphased.vcf.gz"
    )
    phased_phased = f"{tmp_dir.name}/out_diploid_phased_phased.vcf.gz"

    commands.append(
        cmds.cmd_variant_phaser(
            out_vcf, phased_bed, phased_vcf, phased_ext, kwargs
        )
    )

    if kwargs["tech"] == "ONT":
        commands.append(
            f"bcftools view -T {phased_bed} {phased_vcf} \
            | sentieon util vcfconvert - {phased_phased}"
        )
    commands.append(
        cmds.cmd_bedtools_subtract(
            kwargs.get("bed"), phased_bed, unphased_bed, tmp_dir.name, **kwargs
        )
    )

    commands.append(cmds.cmd_repeat_model(phased_bed, phased_ext, tmp_dir.name, kwargs))
    # TODO: difference with ONT here?
    commands.append(
        f"bcftools view -T {unphased_bed} {phased_vcf} \
        | sentieon util vcfconvert - {phased_unphased}"
    )

    kwargs["phased_bed"] = phased_bed
    kwargs["unphased_bed"] = unphased_bed
    for phase in (1, 2):
        kwargs["read-filter"] = f"PhasedReadFilter,phased_vcf={phased_ext}"
        kwargs["read-filter"] += f",phase_select={phase}"
        hp_std_vcf = f"{tmp_dir.name}/out_hap{phase}_nohp_tmp.vcf.gz"

        if kwargs["tech"] == "HiFi":
            cmd = cmds.cmd_algo_dnascope(
                f"{kwargs['model_bundle']}/haploid_model",
                hp_std_vcf,
                kwargs,
                bed_key="phased_bed",
                gvcf=None,
            )
        else:
            # ONT doesn't do DNAscope in 2nd pass.
            cmd = cmds.cmd_sentieon_driver(
                bed_key="phased_bed",
                **kwargs,
            )
        kwargs.pop("read-filter")
        commands.append("#PHASE %s" % phase)

        hp_vcf = f"{tmp_dir.name}/out_hap{phase}_tmp.vcf.gz"

        cmd += " " + cmds.cmd_dnascope_hp(
            f"{kwargs['model_bundle']}/haploid_hp_model",
            kwargs["repeat_model"],
            hp_vcf,
            kwargs,
        )
        commands.append(cmd)

    # TODO: set these properly.
    kwargs["vcf_mod_py"] = os.path.join(
        os.path.dirname(__file__), "vcf_mod.py"
    )
    kwargs["gvcf_combine_py"] = os.path.join(
        os.path.dirname(__file__), "gvcf_combine.py"
    )
    commands.append(
        cmds.cmd_pyexec_vcf_mod_haploid_patch(
            f"{tmp_dir.name}/out_hap1_patch.vcf.gz",
            f"{tmp_dir.name}/out_hap2_patch.vcf.gz",
            f"{tmp_dir.name}/out_hap%d_%stmp.vcf.gz",
            kwargs["tech"],
            phased_phased,
            kwargs,
        )
    )

    # apply trained model to the patched vcfs.
    for hap in (1, 2):
        hap_out = f"{tmp_dir.name}/out_hap{hap}.vcf.gz"
        commands.append(
            cmds.cmd_model_apply(
                f"{kwargs['model_bundle']}/haploid_model",
                f"{tmp_dir.name}/out_hap{hap}_patch.vcf.gz",
                hap_out,
                kwargs,
            )
        )

    # call variants on unphased regions
    # https://github.com/Sentieon/sentieon-scripts/blob/master/dnascope_LongRead/dnascope_HiFi.sh#L409
    cmd = cmds.cmd_sentieon_driver(
        bed_key="unphased_bed",
        skip_sample_input=False,
        **kwargs,
    )
    cmd += " " + cmds.cmd_dnascope_hp(
        f"{kwargs['model_bundle']}/diploid_hp_model",
        kwargs["repeat_model"],
        # f"{tmp_dir.name}/out_diploid_unphased.bed",
        f"{tmp_dir.name}/out_diploid_phased_unphased_hp.vcf.gz",
        kwargs,
    )

    commands.append(cmd)

    # Patch DNA and DNAHP variants
    cmd = cmds.cmd_pyexec_vcf_mod_patch(
        f"{tmp_dir.name}/out_diploid_unphased_patch.vcf.gz",
        phased_unphased,
        f"{tmp_dir.name}/out_diploid_phased_unphased_hp.vcf.gz",
        kwargs,
    )
    commands.append(cmd)
    cmd = cmds.cmd_model_apply(
        f"{kwargs['model_bundle']}/diploid_model_unphased",
        f"{tmp_dir.name}/out_diploid_unphased_patch.vcf.gz",
        f"{tmp_dir.name}/out_diploid_unphased.vcf.gz",
        kwargs,
    )
    commands.append(cmd)

    # merge calls to create the output
    # https://github.com/Sentieon/sentieon-scripts/blob/master/dnascope_LongRead/dnascope_HiFi.sh#L421

    commands.append(
        cmds.cmd_pyexec_vcf_mod_merge(
            f"{tmp_dir.name}/out_hap1.vcf.gz",
            f"{tmp_dir.name}/out_hap2.vcf.gz",
            f"{tmp_dir.name}/out_diploid_unphased.vcf.gz",
            phased_vcf,
            f"{tmp_dir.name}/out_diploid_phased.bed",
            cmds.name(kwargs["output-vcf"]),
            kwargs,
        )
    )

    if gvcf:
        commands.append(
            cmds.cmd_pyexec_gvcf_combine(
                gvcf,
                cmds.name(kwargs["output-vcf"]),
                kwargs,
            )
        )

    return "\n".join(commands)


def main():
    """main entry point for this project"""
    logger.setLevel(os.environ.get("LOGLEVEL", "DEBUG").upper())
    logger.info("Starting sentieon-cli version: %s", __version__)
    argh.dispatch_commands([run_dnascope_longread])


if __name__ == "__main__":
    main()
