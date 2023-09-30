import multiprocessing as mp
import argparse
import os
import subprocess as sp
import argh
import random
from argh import arg


def tmp():
    """Create a temporary directory for the current process."""
    tmp_base = os.getenv("SENTIEON_TMPDIR", os.getenv("TMPDIR", "/tmp"))
    tmp_base = os.path.join(tmp_base, str(os.getpid()))
    tmp_dir = tmp_base
    while os.path.exists(tmp_dir):
        tmp_dir = tmp_base + "_" + str(random.randint(0, 1000000))
    os.makedirs(tmp_dir)
    return tmp_dir


def check_version(cmd: str, version: str):
    cmd_version = sp.check_output([cmd, "--version"]).decode("utf-8").strip()
    # TODO: use regexp
    # https://github.com/Sentieon/sentieon-scripts/blob/master/dnascope_LongRead/dnascope_HiFi.sh#L146
    assert cmd_version


def run_algo_dnascope():
    # https://github.com/Sentieon/sentieon-scripts/blob/8d33f29e442d5a1e782445f06bc1f11e278d8f87/dnascope_LongRead/dnascope_HiFi.sh#L355-L359
    pass


def run_algo_variantphaser():
    # https://github.com/Sentieon/sentieon-scripts/blob/8d33f29e442d5a1e782445f06bc1f11e278d8f87/dnascope_LongRead/dnascope_HiFi.sh#L364-L367
    pass


def common(fn):
    fn = arg("path", help="some path")(fn)
    fn = arg("--extra", help="extra integer arguments: %(default)s", default=0)(fn)
    return fn


@common
def other(**kwargs):
    print("OTHER")


@arg(
    "output-vcf",
    help="Outpu VCF File. The file name must end in .vcf.gz",
    type=argparse.FileType("r"),
)
@arg("reference", help="fasta for reference genome", type=argparse.FileType("r"))
@arg("sample_input", help="sample BAM or CRAM file", type=argparse.FileType("r"))
@arg("model_file", help="The model bundle file", type=argparse.FileType("r"))
@arg(
    "-d",
    "--dbsnp",
    help="dbSNP vcf file Supplying this file will annotate variants with \
         their dbSNP refSNP ID numbers.",
    default="",
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
def c_main(**kwargs):
    assert kwargs["output_vcf"].endswith(".vcf.gz"), "Output file must end in .vcf.gz"


def main():
    argh.dispatch_commands([c_main, other])


if __name__ == "__main__":
    main()
