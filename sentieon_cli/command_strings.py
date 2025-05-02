"""
This module contains the functions that accept arguments and return the
command strings.
"""

import io
import pathlib
import shlex
import subprocess as sp
import typing
from typing import Any, Optional, List, Dict

from .driver import BaseDriver
from .logging import get_logger

logger = get_logger(__name__)


def cmd_fai_to_bed(
    fai: pathlib.Path,
    bed: pathlib.Path,
) -> str:
    cmd = ["awk", "-v", "OFS=\t", "{print $1,0,$2}", str(fai)]
    return shlex.join(cmd) + " > " + shlex.quote(str(bed))


def cmd_bedtools_subtract(
    regions_bed: pathlib.Path,
    phased_bed: pathlib.Path,
    unphased_bed: pathlib.Path,
) -> str:
    cmd = [
        "bedtools",
        "subtract",
        "-a",
        str(regions_bed),
        "-b",
        str(phased_bed),
    ]
    return shlex.join(cmd) + " > " + shlex.quote(str(unphased_bed))


def cmd_bedtools_merge(
    in_bed: pathlib.Path,
    out_bed: pathlib.Path,
    distance: int = 0,
) -> str:
    """Bedtools merge"""
    cmd = [
        "bedtools",
        "merge",
        "-d",
        str(distance),
        "-i",
        str(in_bed),
    ]
    return shlex.join(cmd) + " > " + shlex.quote(str(out_bed))


def cmd_bedtools_slop(
    in_bed: pathlib.Path,
    out_bed: pathlib.Path,
    ref_fai: pathlib.Path,
    slop_size: int = 0,
) -> str:
    """Bedtools slop"""
    cmd = [
        "bedtools",
        "slop",
        "-b",
        str(slop_size),
        "-g",
        str(ref_fai),
        "-i",
        str(in_bed),
    ]
    return shlex.join(cmd) + " > " + shlex.quote(str(out_bed))


def cmd_bedtools_cat_sort_merge(
    out_bed: pathlib.Path,
    in_bed: List[pathlib.Path],
    ref_fai: pathlib.Path,
) -> str:
    """Concat, sort and merge multiple BEDs"""
    cmds = []
    cmds.append(["cat"] + [shlex.quote(str(x)) for x in in_bed])
    cmds.append(
        [
            "bedtools",
            "sort",
            "-faidx",
            str(ref_fai),
            "-i",
            "-",
        ]
    )
    cmds.append(
        [
            "bedtools",
            "merge",
        ]
    )
    all_cmds = [shlex.join(x) for x in cmds]
    return " | ".join(all_cmds) + " > " + shlex.quote(str(out_bed))


def name(path: typing.Union[str, io.TextIOWrapper, pathlib.Path]) -> str:
    """Return the name of a file that may also be a text-wrapper."""
    if isinstance(path, io.TextIOWrapper):
        return path.name
    elif isinstance(path, pathlib.Path):
        return str(path)
    return path


def cmd_pyexec_vcf_mod_haploid_patch(
    hap1_patch: str,
    hap2_patch: str,
    hap_patt: str,
    tech: str,
    phased_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """
    merge dnascope and dnascope-hp variants

    >>> assert '--hap1 ' in cmd_pyexec_vcf_mod_haploid_patch("h1.vcf.gz",
    ... "h2.vcf.gz", "out_hap%d_%stmp.vcf.gz", "HiFi", "ph.vcf", 2,
    ... {'vcf_mod_py': 'vcf_mod.py'})
    >>> assert not "--hap1 " in cmd_pyexec_vcf_mod_haploid_patch("h1.vcf.gz",
    ... "h2.vcf.gz", "out_hap%d_%stmp.vcf.gz", "ONT", "ph.vcf", 2,
    ... {'vcf_mod_py': 'vcf_mod.py'})

    """
    assert tech.upper() in ("HIFI", "ONT")

    cmd = f"sentieon pyexec {kwargs['vcf_mod_py']} -t {cores} "
    cmd += "haploid_patch "
    cmd += f"--patch1 {hap1_patch} --patch2 {hap2_patch}"
    if tech.upper() == "HIFI":
        cmd += " --hap1 " + hap_patt % (1, "nohp_")
        cmd += " --hap2 " + hap_patt % (2, "nohp_")
    else:
        cmd += f" --phased {phased_vcf}"
    cmd += " --hap1_hp " + hap_patt % (1, "")
    cmd += " --hap2_hp " + hap_patt % (2, "")
    return cmd


def cmd_pyexec_vcf_mod_patch(
    out_vcf: str,
    vcf: str,
    vcf_hp: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Patch DNAscope and DNAscopeHP VCF files"""

    cmd = [
        "sentieon",
        "pyexec",
        str(kwargs["vcf_mod_py"]),
        "-t",
        str(cores),
        "patch",
        "--vcf",
        str(vcf),
        "--vcf_hp",
        str(vcf_hp),
        str(out_vcf),
    ]
    return shlex.join(cmd)


def cmd_pyexec_gvcf_combine(
    reference: pathlib.Path,
    gvcf: str,
    out_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Combine gVCF files"""

    cmd1 = [
        "sentieon",
        "pyexec",
        str(kwargs["gvcf_combine_py"]),
        "-t",
        str(cores),
        str(reference),
        gvcf,
        out_vcf,
        "-",
    ]
    cmd2 = [
        "sentieon",
        "util",
        "vcfconvert",
        "-",
        out_vcf.replace(".vcf.gz", ".g.vcf.gz"),
    ]
    return shlex.join(cmd1) + "|" + shlex.join(cmd2)


def cmd_pyexec_vcf_mod_merge(
    hap1_vcf: str,
    hap2_vcf: str,
    unphased_vcf: str,
    phased_vcf: str,
    phased_bed: str,
    out_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Merge haploid VCF files"""

    cmd = [
        "sentieon",
        "pyexec",
        kwargs["vcf_mod_py"],
        "-t",
        str(cores),
        "merge",
        "--hap1",
        hap1_vcf,
        "--hap2",
        hap2_vcf,
        "--unphased",
        unphased_vcf,
        "--phased",
        phased_vcf,
        "--bed",
        phased_bed,
        out_vcf,
    ]
    return shlex.join(cmd)


def cmd_pyexec_vcf_mod_haploid_patch2(
    out_vcf: str,
    vcf: str,
    vcf_hp: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Patch a single pair of haploid DNAscope and DNAscopeHP VCFs"""

    cmd = [
        "sentieon",
        "pyexec",
        str(kwargs["vcf_mod_py"]),
        "-t",
        str(cores),
        "haploid_patch2",
        "--vcf",
        vcf,
        "--vcf_hp",
        vcf_hp,
        "--patch_vcf",
        out_vcf,
    ]
    return shlex.join(cmd)


def cmd_pyexec_hybrid_select(
    out_bed: pathlib.Path,
    vcf: pathlib.Path,
    ref_fai: pathlib.Path,
    hybrid_select: pathlib.Path,
    threads: int,
    slop_size: int = 1000,
) -> str:
    cmds = []
    cmds.append(
        [
            "sentieon",
            "pyexec",
            str(hybrid_select),
            "-v",
            str(vcf),
            "-t",
            str(threads),
            "-",
        ]
    )
    cmds.append(
        [
            "bcftools",
            "view",
            "-f",
            "PASS,.",
            "-",
        ]
    )
    cmds.append(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\\t%POS0\\t%END\\n",
            "-",
        ]
    )
    cmds.append(
        [
            "bedtools",
            "slop",
            "-b",
            str(slop_size),
            "-g",
            str(ref_fai),
            "-i",
            "-",
        ]
    )
    return " | ".join([shlex.join(x) for x in cmds]) + " > " + str(out_bed)


def cmd_pyexec_hybrid_anno(
    out_vcf: pathlib.Path,
    vcf: pathlib.Path,
    bed: pathlib.Path,
    hybrid_anno: pathlib.Path,
    threads: int,
) -> str:
    cmd = [
        "sentieon",
        "pyexec",
        str(hybrid_anno),
        "-v",
        str(vcf),
        "-b",
        str(bed),
        "-t",
        str(threads),
        str(out_vcf),
    ]
    return shlex.join(cmd)


def hybrid_stage1(
    out_aln: pathlib.Path,
    reference: pathlib.Path,
    cores: int,
    readgroup: str,
    ins_driver: BaseDriver,
    stage1_driver: BaseDriver,
    bwa_model: pathlib.Path,
) -> str:
    unset_cmd = ["unset", "bwt_max_mem"]
    fq_cmds = []
    fq_cmds.append(shlex.join(stage1_driver.build_cmd()))
    fq_cmds.append(shlex.join(ins_driver.build_cmd()))

    aln_cmds = []
    aln_cmds.append(
        [
            "sentieon",
            "bwa",
            "mem",
            "-R",
            readgroup,
            "-t",
            str(cores),
            "-x",
            str(bwa_model),
            str(reference),
            "-",
        ]
    )
    aln_cmds.append(
        [
            "sentieon",
            "util",
            "sort",
            "-i",
            "-",
            "-t",
            str(cores),
            "-o",
            str(out_aln),
            "--sam2bam",
        ]
    )

    return (
        shlex.join(unset_cmd)
        + "; ("
        + " || echo -n 'error'; ".join(fq_cmds)
        + " || echo -n 'error') | "
        + " | ".join(shlex.join(x) for x in aln_cmds)
    )


def hybrid_stage3(
    out_bam: pathlib.Path,
    driver: BaseDriver,
    cores: int,
) -> str:
    cmds = []
    cmds.append(driver.build_cmd())
    cmds.append(
        [
            "sentieon",
            "util",
            "sort",
            "-i",
            "-",
            "-t",
            str(cores),
            "-o",
            str(out_bam),
        ]
    )
    return " | ".join([shlex.join(x) for x in cmds])


def bcftools_subset(
    out_vcf: pathlib.Path,
    mix_vcf: pathlib.Path,
    regions_bed: pathlib.Path,
) -> str:
    """Subset a vcf with bcftools"""
    cmds = []
    cmds.append(
        [
            "bcftools",
            "view",
            "-T",
            "^" + str(regions_bed),
            str(mix_vcf),
        ]
    )
    cmds.append(
        [
            "sentieon",
            "util",
            "vcfconvert",
            "-",
            str(out_vcf),
        ]
    )
    return " | ".join([shlex.join(x) for x in cmds])


def bcftools_concat(
    out_vcf: pathlib.Path,
    in_vcfs: List[pathlib.Path],
) -> str:
    """VCF processing through bcftools concat"""
    cmds = []
    cmds.append(
        [
            "bcftools",
            "concat",
            "-aD",
        ]
        + [str(x) for x in in_vcfs]
    )
    cmds.append(
        [
            "sentieon",
            "util",
            "vcfconvert",
            "-",
            str(out_vcf),
        ]
    )
    return " | ".join([shlex.join(x) for x in cmds])


def filter_norm(
    out_vcf: pathlib.Path,
    in_vcf: pathlib.Path,
    reference: pathlib.Path,
    exclude_homref: bool = True,
) -> str:
    """Trim and normalize"""
    cmds = []
    cmd = [
        "bcftools",
        "view",
        "-a",
    ]
    if exclude_homref:
        cmd.extend(
            [
                "-e",
                'GT="0/0"',
            ]
        )
    cmd.append(str(in_vcf))
    cmds.append(cmd)
    cmds.append(
        [
            "bcftools",
            "norm",
            "-f",
            str(reference),
        ]
    )
    cmds.append(
        [
            "sentieon",
            "util",
            "vcfconvert",
            "-",
            str(out_vcf),
        ]
    )
    return " | ".join([shlex.join(x) for x in cmds])


def get_rg_lines(
    input_aln: pathlib.Path,
    dry_run: bool,
) -> List[str]:
    """Read the @RG lines from an alignment file"""
    cmd = shlex.join(
        [
            "samtools",
            "view",
            "-H",
            str(input_aln),
        ]
    )
    rg_lines: List[str] = []
    if not dry_run:
        res = sp.run(cmd, shell=True, check=True, stdout=sp.PIPE)
        for line in res.stdout.decode("utf-8").split("\n"):
            if line.startswith("@RG"):
                rg_lines.append(line)
    else:
        rg_lines.append("@RG\tID:mysample-1\tSM:mysample")
    return rg_lines


def rehead_wgsmetrics(
    orig_metrics: pathlib.Path, tmp_metrics: pathlib.Path
) -> str:
    """Rehead Sentieon WGS metrics so the file is recognized by MultiQC"""
    cmd1 = shlex.join(["mv", str(orig_metrics), str(tmp_metrics)])
    cmd2 = (
        shlex.join(["echo", "'## METRICS CLASS WgsMetrics'"])
        + ">"
        + shlex.quote(str(orig_metrics))
    )
    cmd3 = (
        shlex.join(["tail", "-n", "+2", str(tmp_metrics)])
        + ">>"
        + shlex.quote(str(orig_metrics))
    )
    cmd4 = shlex.join(["rm", str(tmp_metrics)])
    return "; ".join((cmd1, cmd2, cmd3, cmd4))


def cmd_samtools_fastq_minimap2(
    out_aln: pathlib.Path,
    input_aln: pathlib.Path,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int,
    rg_lines: List[str],
    sample_name: str,
    input_ref: Optional[pathlib.Path] = None,
    fastq_taglist: str = "*",
    minimap2_args: str = "-Y",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
) -> str:
    """Re-align an input BAM/CRAM/uBAM/uCRAM file with minimap2"""

    ref_cmd: List[str] = []
    if input_ref:
        ref_cmd = ["--reference", str(reference)]
    cmd1 = (
        [
            "samtools",
            "fastq",
        ]
        + ref_cmd
        + [
            "-@",
            str(cores),
            "-T",
            fastq_taglist,
            str(input_aln),
        ]
    )
    cmd2 = [
        "sentieon",
        "minimap2",
        "-y",
        "-t",
        str(cores),
        "-a",
        minimap2_args,
        "-x",
        f"{model_bundle}/minimap2.model",
        str(reference),
        "/dev/stdin",
    ]
    cmd3 = [
        "sentieon",
        "util",
        "sort",
        "-i",
        "-",
        "-t",
        str(cores),
        "--reference",
        str(reference),
        "-o",
        str(out_aln),
        "--sam2bam",
    ] + util_sort_args.split()

    # Commands to replace the @RG lines in the header
    rg_cmds: List[List[str]] = []
    for rg_line in rg_lines:
        # Add a SM value, if missing
        if "\tSM:" not in rg_line:
            rg_line += f"\tSM:{sample_name}"
        rg_cmds.append(
            [
                "samtools",
                "addreplacerg",
                "-r",
                str(rg_line),
                "-m",
                "orphan_only",
                "-",
            ]
        )

    return " | ".join([shlex.join(x) for x in (cmd1, cmd2, *rg_cmds, cmd3)])


def cmd_samtools_fastq_bwa(
    out_aln: pathlib.Path,
    input_aln: pathlib.Path,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int,
    rg_header: pathlib.Path,
    input_ref: Optional[pathlib.Path] = None,
    collate: bool = False,
    bwa_args: str = "",
    bwa_k: str = "100000000",
    fastq_taglist: str = "RG",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
) -> str:
    """Re-align an input BAM/CRAM/uBAM/uCRAM file with bwa"""
    ref_cmd: List[str] = []
    if input_ref:
        ref_cmd = ["--reference", str(reference)]

    cmd0 = ["perl", "-MFcntl", "-e", "fcntl(STDOUT, 1031, 268435456)"]
    if collate:
        cmd1 = (
            [
                "samtools",
                "collate",
            ]
            + ref_cmd
            + [
                "-@",
                str(cores),
                "-O",
                "-u",
                "-f",
                str(input_aln),
            ]
        )
        cmd2 = [
            "samtools",
            "fastq",
            "-@",
            str(cores),
            "-T",
            fastq_taglist,
            "-",
        ]
    else:
        cmd1 = []
        cmd2 = (
            [
                "samtools",
                "fastq",
            ]
            + ref_cmd
            + [
                "-@",
                str(cores),
                "-T",
                fastq_taglist,
                str(input_aln),
            ]
        )
    cmd3 = (
        [
            "sentieon",
            "bwa",
            "mem",
            "-p",
            "-C",
            "-H",
            str(rg_header),
            "-t",
            str(cores),
            "-x",
            f"{model_bundle}/bwa.model",
        ]
        + bwa_args.split()
        + ["-K", bwa_k]
        + [str(reference), "/dev/stdin"]
    )
    cmd4 = [
        "sentieon",
        "util",
        "sort",
        "-i",
        "-",
        "-t",
        str(cores),
        "--reference",
        str(reference),
        "-o",
        str(out_aln),
        "--sam2bam",
    ] + util_sort_args.split()

    if cmd1:
        return (
            shlex.join(cmd0)
            + ";"
            + " | ".join([shlex.join(x) for x in (cmd1, cmd2, cmd3, cmd4)])
        )
    else:
        return (
            shlex.join(cmd0)
            + ";"
            + " | ".join([shlex.join(x) for x in (cmd2, cmd3, cmd4)])
        )


def cmd_fastq_minimap2(
    out_aln: pathlib.Path,
    fastq: pathlib.Path,
    readgroup: str,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int,
    unzip: str = "gzip",
    minimap2_args: str = "-Y",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
) -> str:
    """Align an input fastq file with minimap2"""

    cmd1 = [
        unzip,
        "-dc",
        str(fastq),
    ]
    cmd2 = [
        "sentieon",
        "minimap2",
        "-t",
        str(cores),
        "-a",
        minimap2_args,
        "-x",
        f"{model_bundle}/minimap2.model",
        "-R",
        readgroup,
        str(reference),
        "/dev/stdin",
    ]
    cmd3 = [
        "sentieon",
        "util",
        "sort",
        "-i",
        "-",
        "-t",
        str(cores),
        "--reference",
        str(reference),
        "-o",
        str(out_aln),
        "--sam2bam",
    ] + util_sort_args.split()
    return " | ".join([shlex.join(x) for x in (cmd1, cmd2, cmd3)])


def cmd_fastq_bwa(
    out_aln: pathlib.Path,
    r1: pathlib.Path,
    r2: Optional[pathlib.Path],
    readgroup: str,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int,
    unzip: str = "gzip",
    bwa_args: str = "",
    bwa_k: str = "100000000",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
    numa: Optional[str] = None,
    split: Optional[str] = None,
) -> str:
    """Align an input fastq file with bwa"""
    pipebuf_cmd = shlex.join(
        ["perl", "-MFcntl", "-e", "fcntl(STDOUT, 1031, 268435456)"]
    )
    cmd1 = (
        (["taskset", "-c", numa] if numa else [])
        + [
            "sentieon",
            "bwa",
            "mem",
            "-R",
            readgroup,
            "-t",
            str(cores),
            "-x",
            f"{model_bundle}/bwa.model",
        ]
        + bwa_args.split()
        + ["-K", bwa_k]
        + (["-p"] if (r2 and split) else [])
        + [
            str(reference),
        ]
    )
    cmd2 = [
        str(unzip),
        "-dc",
        str(r1),
    ]
    cmd3 = []
    if r2:
        cmd3 = [
            str(unzip),
            "-dc",
            str(r2),
        ]
    cmd4 = (
        (["taskset", "-c", numa] if numa else [])
        + [
            "sentieon",
            "util",
            "sort",
            "-i",
            "-",
            "-t",
            str(cores),
            "--reference",
            str(reference),
            "-o",
            str(out_aln),
            "--sam2bam",
        ]
        + util_sort_args.split()
    )

    cmds = [shlex.join(x) for x in (cmd1, cmd2, cmd3, cmd4)]
    cmd_str = ""
    if split:
        extract_cmd = shlex.join(
            ["sentieon", "fqidx", "extract", "-F", str(split), "-K", bwa_k]
        )
        cmd_str = (
            pipebuf_cmd
            + "; "
            + cmds[0]
            + " <("
            + pipebuf_cmd
            + "; "
            + extract_cmd
            + " <("
            + pipebuf_cmd
            + "; "
            + cmds[1]
            + ") "
        )
        if cmds[2]:
            cmd_str += " <(" + pipebuf_cmd + "; " + cmds[2] + ") "
        cmd_str += ") "
    else:
        cmd_str = (
            pipebuf_cmd
            + "; "
            + cmds[0]
            + " <("
            + pipebuf_cmd
            + "; "
            + cmds[1]
            + ") "
        )
        if cmds[2]:
            cmd_str += " <(" + pipebuf_cmd + "; " + cmds[2] + ") "
    cmd_str += " | " + cmds[3]
    return cmd_str


def cmd_multiqc(
    input_directory: pathlib.Path,
    output_directory: pathlib.Path,
    comment: Optional[str],
) -> str:
    """Run MultiQC"""
    cmd = [
        "multiqc",
        "-o",
        str(output_directory),
    ]
    if comment:
        cmd.extend(
            [
                "-b",
                comment,
            ]
        )
    cmd.extend(
        [
            str(input_directory),
        ]
    )
    return shlex.join(cmd)


def cmd_mosdepth(
    sample_input: pathlib.Path,
    output_directory: pathlib.Path,
    fasta: Optional[pathlib.Path] = None,
    threads: int = 1,
    xargs: str = (
        "--by 500 --no-per-base --use-median -T 1,3,5,10,15,20,30,40,50"
    ),
) -> str:
    cmd = [
        "mosdepth",
        "--fasta",
        str(fasta),
        "--threads",
        str(threads),
    ]
    cmd.extend(shlex.split(xargs))
    cmd.append(str(output_directory))
    cmd.append(str(sample_input))
    return shlex.join(cmd)
