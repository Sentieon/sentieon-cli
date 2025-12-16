"""
This module contains the functions that accept arguments and return the
command strings.
"""

import io
import os
import pathlib
import shlex
import subprocess as sp
import typing
from typing import Any, Dict, List, Optional, Union

from .driver import BaseDriver
from .logging import get_logger
from .shell_pipeline import Command, InputProcSub, Pipeline

logger = get_logger(__name__)


def cmd_fai_to_bed(
    fai: pathlib.Path,
    bed: pathlib.Path,
) -> Pipeline:
    cmd = ["awk", "-v", "OFS=\t", "{print $1,0,$2}", str(fai)]
    return Pipeline(Command(*cmd), file_output=bed)


def cmd_bedtools_subtract(
    regions_bed: pathlib.Path,
    phased_bed: pathlib.Path,
    unphased_bed: pathlib.Path,
) -> Pipeline:
    cmd = [
        "bedtools",
        "subtract",
        "-a",
        str(regions_bed),
        "-b",
        str(phased_bed),
    ]
    return Pipeline(Command(*cmd), file_output=unphased_bed)


def cmd_bedtools_merge(
    in_bed: pathlib.Path,
    out_bed: pathlib.Path,
    distance: int = 0,
) -> Pipeline:
    """Bedtools merge"""
    cmd = [
        "bedtools",
        "merge",
        "-d",
        str(distance),
        "-i",
        str(in_bed),
    ]
    return Pipeline(Command(*cmd), file_output=out_bed)


def cmd_bedtools_slop(
    in_bed: pathlib.Path,
    out_bed: pathlib.Path,
    ref_fai: pathlib.Path,
    slop_size: int = 0,
) -> Pipeline:
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
    return Pipeline(Command(*cmd), file_output=out_bed)


def cmd_bedtools_cat_sort_merge(
    out_bed: pathlib.Path,
    in_bed: List[pathlib.Path],
    ref_fai: pathlib.Path,
) -> Pipeline:
    """Concat, sort and merge multiple BEDs"""
    cat_cmd = Command("cat", *[str(x) for x in in_bed])
    sort_cmd = Command("bedtools", "sort", "-faidx", str(ref_fai), "-i", "-")
    merge_cmd = Command("bedtools", "merge")
    return Pipeline(cat_cmd, sort_cmd, merge_cmd, file_output=out_bed)


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
) -> Pipeline:
    """
    merge dnascope and dnascope-hp variants
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
    return Pipeline(Command(*shlex.split(cmd)))


def cmd_pyexec_vcf_mod_patch(
    out_vcf: str,
    vcf: str,
    vcf_hp: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> Pipeline:
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
    return Pipeline(Command(*cmd))


def cmd_pyexec_gvcf_combine(
    reference: pathlib.Path,
    gvcf: str,
    out_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> Pipeline:
    """Combine gVCF files"""

    cmd1 = Command(
        "sentieon",
        "pyexec",
        str(kwargs["gvcf_combine_py"]),
        "-t",
        str(cores),
        str(reference),
        gvcf,
        out_vcf,
        "-",
    )
    cmd2 = Command(
        "sentieon",
        "util",
        "vcfconvert",
        "-",
        out_vcf.replace(".vcf.gz", ".g.vcf.gz"),
    )
    return Pipeline(cmd1, cmd2)


def cmd_pyexec_vcf_mod_merge(
    hap1_vcf: str,
    hap2_vcf: str,
    unphased_vcf: str,
    phased_vcf: str,
    phased_bed: str,
    out_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> Pipeline:
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
    return Pipeline(Command(*cmd))


def cmd_pyexec_vcf_mod_haploid_patch2(
    out_vcf: str,
    vcf: str,
    vcf_hp: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> Pipeline:
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
    return Pipeline(Command(*cmd))


def cmd_pyexec_hybrid_select(
    out_bed: pathlib.Path,
    vcf: pathlib.Path,
    ref_fai: pathlib.Path,
    hybrid_select: pathlib.Path,
    threads: int,
    slop_size: int = 1000,
) -> Pipeline:
    select_cmd = Command(
        "sentieon",
        "pyexec",
        str(hybrid_select),
        "-v",
        str(vcf),
        "-t",
        str(threads),
        "-",
    )
    view_cmd = Command(
        "bcftools",
        "view",
        "-f",
        "PASS,.",
        "-",
    )
    query_cmd = Command(
        "bcftools",
        "query",
        "-f",
        "%CHROM\\t%POS0\\t%END\\n",
        "-",
    )
    slop_cmd = Command(
        "bedtools",
        "slop",
        "-b",
        str(slop_size),
        "-g",
        str(ref_fai),
        "-i",
        "-",
    )
    return Pipeline(
        select_cmd,
        view_cmd,
        query_cmd,
        slop_cmd,
        file_output=out_bed,
    )


def cmd_pyexec_hybrid_anno(
    out_vcf: pathlib.Path,
    vcf: pathlib.Path,
    bed: pathlib.Path,
    hybrid_anno: pathlib.Path,
    threads: int,
) -> Pipeline:
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
    return Pipeline(Command(*cmd))


def hybrid_stage1(
    out_aln: pathlib.Path,
    reference: pathlib.Path,
    cores: int,
    readgroup: str,
    ins_driver: BaseDriver,
    stage1_driver: BaseDriver,
    bwa_model: pathlib.Path,
) -> Pipeline:
    bwa_env = dict(os.environ)
    _ = bwa_env.pop("bwt_max_mem", None)

    # Send the input of both fq commands to bwa with cat
    fq1_cmd = Command(*stage1_driver.build_cmd())
    fq2_cmd = Command(*ins_driver.build_cmd())
    cat_cmd = Command(
        "cat",
        InputProcSub(Pipeline(fq1_cmd)),
        InputProcSub(Pipeline(fq2_cmd)),
    )

    bwa_cmd = Command(
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
        exec_kwargs={"env": bwa_env},
    )

    sort_cmd = Command(
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
    )

    return Pipeline(
        cat_cmd,
        bwa_cmd,
        sort_cmd,
    )


def hybrid_stage3(
    out_bam: pathlib.Path,
    driver: BaseDriver,
    cores: int,
) -> Pipeline:
    sort_cmd = Command(
        "sentieon",
        "util",
        "sort",
        "-i",
        "-",
        "-t",
        str(cores),
        "-o",
        str(out_bam),
    )
    return Pipeline(Command(*driver.build_cmd()), sort_cmd)


def bcftools_subset(
    out_vcf: pathlib.Path,
    mix_vcf: pathlib.Path,
    regions_bed: pathlib.Path,
) -> Pipeline:
    """Subset a vcf with bcftools"""
    bcftools_subset_cmd = shlex.join(
        [
            "bcftools",
            "view",
            "-T",
            "^" + str(regions_bed),
            str(mix_vcf),
        ]
    )
    vcfconvert_cmd = shlex.join(
        [
            "sentieon",
            "util",
            "vcfconvert",
            "-",
            str(out_vcf),
        ]
    )
    cp_cmd = shlex.join(
        [
            "sentieon",
            "util",
            "vcfconvert",
            str(mix_vcf),
            str(out_vcf),
        ]
    )
    shell_cmd = (
        f'if [ -s "{regions_bed}" ]; then   '  # noqa
        + bcftools_subset_cmd
        + " | "
        + vcfconvert_cmd
        + "; else   "
        + cp_cmd
        + "; fi"
    )
    return Pipeline(Command("sh", "-c", shell_cmd))


def bcftools_concat(
    out_vcf: pathlib.Path,
    in_vcfs: List[pathlib.Path],
    xargs: List[str] = ["-aD"],
) -> Pipeline:
    """VCF processing through bcftools concat"""
    concat_cmd = Command(
        "bcftools",
        "concat",
        "-W=tbi",
        "--output",
        str(out_vcf),
        *xargs,
        *[str(x) for x in in_vcfs],
    )
    return Pipeline(concat_cmd)


def filter_norm(
    out_vcf: pathlib.Path,
    in_vcf: pathlib.Path,
    reference: pathlib.Path,
    exclude_homref: bool = True,
) -> Pipeline:
    """Trim and normalize"""
    view_args = ["-a"]
    if exclude_homref:
        view_args.extend(["-e", 'GT="0/0"'])
    view_args.append(str(in_vcf))
    view_cmd = Command("bcftools", "view", *view_args)
    norm_cmd = Command("bcftools", "norm", "-f", str(reference))
    convert_cmd = Command("sentieon", "util", "vcfconvert", "-", str(out_vcf))
    return Pipeline(view_cmd, norm_cmd, convert_cmd)


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
    minimap2_args: str = "-YL",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
) -> Pipeline:
    """Re-align an input BAM/CRAM/uBAM/uCRAM file with minimap2"""

    ref_cmd: List[str] = []
    if input_ref:
        ref_cmd = ["--reference", str(reference)]
    cmd_list = [
        Command(
            "samtools",
            "fastq",
            *ref_cmd,
            "-@",
            str(cores),
            "-T",
            fastq_taglist,
            str(input_aln),
        )
    ]
    cmd_list.append(
        Command(
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
        )
    )

    # Commands to replace the @RG lines in the header
    for rg_line in rg_lines:
        # Add a SM value, if missing
        if "\tSM:" not in rg_line:
            rg_line += f"\tSM:{sample_name}"
        cmd_list.append(
            Command(
                "samtools",
                "addreplacerg",
                "-r",
                str(rg_line),
                "-m",
                "orphan_only",
                "-",
            )
        )
    cmd_list.append(
        Command(
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
            *util_sort_args.split(),
        )
    )
    return Pipeline(*cmd_list)


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
) -> Pipeline:
    """Re-align an input BAM/CRAM/uBAM/uCRAM file with bwa"""
    ref_cmd: List[str] = []
    if input_ref:
        ref_cmd = ["--reference", str(reference)]

    pipebuf_cmd = Command(
        "perl", "-MFcntl", "-e", "fcntl(STDOUT, 1031, 268435456)"
    )
    if collate:
        collate_cmd = Command(
            "samtools",
            "collate",
            *ref_cmd,
            "-@",
            str(cores),
            "-O",
            "-u",
            "-f",
            str(input_aln),
        )
        fastq_cmd = Command(
            "samtools",
            "fastq",
            "-@",
            str(cores),
            "-T",
            fastq_taglist,
            "-",
        )
    else:
        collate_cmd = None
        fastq_cmd = Command(
            "samtools",
            "fastq",
            *ref_cmd,
            "-@",
            str(cores),
            "-T",
            fastq_taglist,
            str(input_aln),
        )
    bwa_cmd = Command(
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
        *bwa_args.split(),
        "-K",
        bwa_k,
        str(reference),
        "/dev/stdin",
    )
    sort_cmd = Command(
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
        *util_sort_args.split(),
    )
    if collate_cmd:
        return Pipeline(
            pipebuf_cmd,
            collate_cmd,
            fastq_cmd,
            bwa_cmd,
            sort_cmd,
            skip_pipe=(0,),
        )
    else:
        return Pipeline(
            pipebuf_cmd,
            fastq_cmd,
            bwa_cmd,
            sort_cmd,
            skip_pipe=(0,),
        )


def cmd_fastq_minimap2(
    out_aln: pathlib.Path,
    fastq: pathlib.Path,
    readgroup: str,
    reference: pathlib.Path,
    model_bundle: pathlib.Path,
    cores: int,
    unzip: str = "gzip",
    minimap2_args: str = "-YL",
    util_sort_args: str = "--cram_write_options version=3.0,compressor=rans",
) -> Pipeline:
    """Align an input fastq file with minimap2"""

    cmd1 = Command(unzip, "-dc", str(fastq))
    cmd2 = Command(
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
    )
    cmd3 = Command(
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
        *util_sort_args.split(),
    )
    return Pipeline(cmd1, cmd2, cmd3)


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
) -> Pipeline:
    """Align an input fastq file with bwa"""
    pipebuf_cmd = Command(
        "perl", "-MFcntl", "-e", "fcntl(STDOUT, 1031, 268435456)"
    )
    r1_unzip = Pipeline(
        pipebuf_cmd,
        Command(str(unzip), "-dc", str(r1)),
        skip_pipe=(0,),
    )
    r2_unzip = None
    if r2:
        r2_unzip = Pipeline(
            pipebuf_cmd,
            Command(str(unzip), "-dc", str(r2)),
            skip_pipe=(0,),
        )
    bwa_cmd_args: List[Union[str, InputProcSub]] = [
        *(["taskset", "-c", numa] if numa else []),
        "sentieon",
        "bwa",
        "mem",
        "-R",
        readgroup,
        "-t",
        str(cores),
        "-x",
        f"{model_bundle}/bwa.model",
        *bwa_args.split(),
        "-K",
        bwa_k,
        *(["-p"] if (r2 and split) else []),
        str(reference),
    ]
    if split:
        extract_cmd = Pipeline(
            pipebuf_cmd,
            Command(
                "sentieon",
                "fqidx",
                "extract",
                "-F",
                str(split),
                "-K",
                bwa_k,
                InputProcSub(r1_unzip),
                *([InputProcSub(r2_unzip)] if r2_unzip else []),
            ),
            skip_pipe=(0,),
        )
        bwa_cmd_args.append(InputProcSub(extract_cmd))
    else:
        bwa_cmd_args.append(InputProcSub(r1_unzip))
        if r2_unzip:
            bwa_cmd_args.append(InputProcSub(r2_unzip))
    sort_cmd = Command(
        *(["taskset", "-c", numa] if numa else []),
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
        *util_sort_args.split(),
    )
    return Pipeline(
        pipebuf_cmd,
        Command(str(bwa_cmd_args[0]), *bwa_cmd_args[1:]),
        sort_cmd,
        skip_pipe=(0,),
    )


def cmd_multiqc(
    input_directory: pathlib.Path,
    output_directory: pathlib.Path,
    comment: Optional[str],
) -> Pipeline:
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
    return Pipeline(Command(*cmd))


def cmd_mosdepth(
    sample_input: pathlib.Path,
    output_directory: pathlib.Path,
    fasta: Optional[pathlib.Path] = None,
    threads: int = 1,
    xargs: str = (
        "--by 500 --no-per-base --use-median -T 1,3,5,10,15,20,30,40,50"
    ),
) -> Pipeline:
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
    return Pipeline(Command(*cmd))


def cmd_kmc(
    output_prefix: pathlib.Path,
    file_list: pathlib.Path,
    tmp_dir: pathlib.Path,
    k=29,
    memory=128,
    threads=1,
) -> Pipeline:
    cmd = [
        "kmc",
        f"-k{k}",
        f"-m{memory}",
        "-okff",
        f"-t{threads}",
        "-hp",
        f"@{file_list}",
        str(output_prefix),
        str(tmp_dir),
    ]
    return Pipeline(Command(*cmd))


def cmd_extract_kmc(
    output_prefix: pathlib.Path,
    out_fastq: pathlib.Path,
    input_aln: List[pathlib.Path],
    reference: pathlib.Path,
    extract_model: pathlib.Path,
    tmp_dir: pathlib.Path,
    threads=1,
) -> Pipeline:
    merge_cmd = Command(
        "samtools",
        "merge",
        "--reference",
        str(reference),
        "-@",
        str(threads),
        "-u",
        "-o",
        "/dev/stdout",
        *[str(x) for x in input_aln],
    )
    view_cmd = Command(
        "samtools",
        "view",
        "-@",
        str(threads),
        "-ub",
        "-F",
        "0xf00",
        "-",
    )
    extract_cmd = Command(
        "sentieon",
        "pgutil",
        "extract",
        "-t",
        str(threads),
        "-T",
        "tp,t0",
        "-S",
        "-PG",
        "-m",
        str(extract_model),
        "-O",
        "fastq_compression=1",
        "-o",
        str(out_fastq),
        "-a",
        "-",
    )
    kmc_cmd = Command(
        "kmc",
        "-k29",
        "-m32",
        "-okff",
        f"-t{threads}",
        "-fa",
        "/dev/stdin",
        str(output_prefix),
        str(tmp_dir),
    )
    return Pipeline(merge_cmd, view_cmd, extract_cmd, kmc_cmd)


def cmd_vg_haplotypes(
    output_gbz: pathlib.Path,
    kmer_file: pathlib.Path,
    hapl: pathlib.Path,
    gbz: pathlib.Path,
    threads=1,
    verbosity=1,
    xargs: List[str] = ["--include-reference", "--diploid-sampling"],
) -> Pipeline:
    cmd = [
        "vg",
        "haplotypes",
        "-g",
        str(output_gbz),
        "-t",
        str(threads),
        "-v",
        str(verbosity),
        "-i",
        str(hapl),
        "-k",
        str(kmer_file),
        str(gbz),
    ]
    cmd.extend(xargs)
    return Pipeline(Command(*cmd))


def cmd_vg_giraffe(
    sample_gam: pathlib.Path,
    sample_pangenome: pathlib.Path,
    fastq1: pathlib.Path,
    fastq2: Optional[pathlib.Path],
    readgroup="rg-1",
    sample="sample",
    threads=1,
    max_fragment_length=3000,
) -> Pipeline:
    vg_cmd = [
        "vg",
        "giraffe",
        "-t",
        str(threads),
        "-Z",
        str(sample_pangenome),
        "--read-group",
        readgroup,
        "--sample",
        sample,
        "-L",
        str(max_fragment_length),
        "--progress",
        "-f",
        str(fastq1),
    ]
    if fastq2:
        vg_cmd.extend(["-f", str(fastq2)])

    return Pipeline(Command(*vg_cmd), file_output=sample_gam)


def cmd_vg_pack(
    sample_pack: pathlib.Path,
    sample_gams: List[pathlib.Path],
    gbz: pathlib.Path,
    min_mapq=5,
    threads=1,
) -> Pipeline:
    cat_cmd = Command("cat", *[str(x) for x in sample_gams])
    vg_cmd = Command(
        "vg",
        "pack",
        "-x",
        str(gbz),
        "-g",
        "-",
        "-o",
        str(sample_pack),
        "-Q",
        str(min_mapq),
        "--threads",
        str(threads),
    )

    return Pipeline(cat_cmd, vg_cmd)


def cmd_sv_call(
    out_svs: pathlib.Path,
    sample_pack: pathlib.Path,
    sv_header: pathlib.Path,
    gbz: pathlib.Path,
    snarls: pathlib.Path,
    sample_name="",
    ref_path="GRCh38",
    min_indel_size=34,
) -> Pipeline:
    # Call SVs with vg call
    vg_call_cmd = Command(
        "vg",
        "call",
        "-r",
        str(snarls),
        "-k",
        str(sample_pack),
        "-s",
        sample_name,
        "-z",
        "--ref-sample",
        ref_path,
        "--progress",
        str(gbz),
    )

    # Use bcftools to remove smaller variants
    view_cmd = Command(
        "bcftools",
        "view",
        "-i",
        f"ABS(ILEN) > {min_indel_size}",
        "-H",
        "-",
    )

    # Rehead the vg call output
    cat_cmd = Command(
        "cat",
        str(sv_header),
        InputProcSub(Pipeline(vg_call_cmd, view_cmd)),
    )

    # Sort and compress the output
    sort_cmd = Command(
        "bcftools",
        "sort",
        "-O",
        "v",
    )
    convert_cmd = Command(
        "sentieon",
        "util",
        "vcfconvert",
        "-",
        str(out_svs),
    )
    return Pipeline(cat_cmd, sort_cmd, convert_cmd)


def cmd_vg_surject(
    out_bam: pathlib.Path,
    sample_gam: pathlib.Path,
    xg: pathlib.Path,
    surject_paths_dict: pathlib.Path,
    threads=1,
    surject_xargs: List[str] = [
        "--interleaved",
        "--progress",
    ],
) -> Pipeline:
    surject_cmd_list = [
        "vg",
        "surject",
        "--sam-output",
        "-t",
        str(threads),
        "-x",
        str(xg),
        "--into-paths",
        shlex.quote(str(surject_paths_dict)),
    ]
    surject_cmd_list.extend(surject_xargs)
    surject_cmd_list.append(str(sample_gam))

    sort_cmd = Command(
        "sentieon",
        "util",
        "sort",
        "--sam2bam",
        "-i",
        "-",
        "-o",
        str(out_bam),
    )

    return Pipeline(Command(*surject_cmd_list), sort_cmd)


def strip_ctg_prefix(
    output_header: pathlib.Path,
    bam: pathlib.Path,
    prefix: str,
) -> Pipeline:
    # Also remove the M5 tag
    samtools_cmd = Command("samtools", "view", "-H", str(bam))
    sed_cmd = Command(
        "sed",
        "-e",
        f"s/\tSN:{prefix}/\tSN:/",
        "-e",
        "s/\tM5:[0-9a-zA-Z]*\t/\t/",
        "-e",
        "s/\tM5:[0-9a-zA-Z]*$//",
    )
    return Pipeline(samtools_cmd, sed_cmd, file_output=output_header)


def cmd_estimate_ploidy(
    output_json: pathlib.Path,
    aln_files: List[pathlib.Path],
    ploidy_script: pathlib.Path,
) -> Pipeline:
    cmd = (
        [
            "python3",
            str(ploidy_script),
            "-i",
        ]
        + [str(x) for x in aln_files]
        + [
            "--outfile",
            str(output_json),
        ]
    )
    return Pipeline(Command(*cmd))


def cmd_t1k(
    out_prefix: pathlib.Path,
    deduped_bam: pathlib.Path,
    gene_seq: pathlib.Path,
    gene_coord: pathlib.Path,
    preset: str,
    threads=1,
) -> Pipeline:
    cmd = [
        "run-t1k",
        "--abnormalUnmapFlag",
        "-t",
        str(threads),
        "--preset",
        preset,
        "-f",
        str(gene_seq),
        "-c",
        str(gene_coord),
        "--od",
        str(out_prefix),
        "-b",
        str(deduped_bam),
    ]
    return Pipeline(Command(*cmd))


def cmd_expansion_hunter(
    out_expansions: pathlib.Path,
    input_alignments: pathlib.Path,
    reference: pathlib.Path,
    variant_catalog: pathlib.Path,
    sex="female",
    threads=1,
) -> Pipeline:
    cmd = [
        "ExpansionHunter",
        "--reads",
        str(input_alignments),
        "--reference",
        str(reference),
        "--variant-catalog",
        str(variant_catalog),
        "--sex",
        sex,
        "--threads",
        str(threads),
        "--output-prefix",
        str(out_expansions),
    ]
    return Pipeline(Command(*cmd))


def cmd_segdup_caller(
    out_segdup: pathlib.Path,
    sr_alignments: pathlib.Path,
    reference: pathlib.Path,
    sr_bundle: pathlib.Path,
    genes="CFH,CFHR3,CYP11B1,CYP2D6,GBA,NCF1,PMS2,SMN1,STRC",
) -> Pipeline:
    cmd = [
        "segdup-caller",
        "--short",
        str(sr_alignments),
        "--reference",
        str(reference),
        "--sr_model",
        str(sr_bundle),
        "--genes",
        genes,
        "--outdir",
        str(out_segdup),
    ]
    return Pipeline(Command(*cmd))


def cmd_bwa_extract(
    out_bam: pathlib.Path,
    out_fastq: pathlib.Path,
    reference: pathlib.Path,
    fastq1: List[pathlib.Path],
    fastq2: List[pathlib.Path],
    readgroup: str,
    extract_model: pathlib.Path,
    bwa_model: pathlib.Path,
    threads=1,
    unzip: str = "gzip",
) -> Pipeline:
    """BWA alignment with pgutil extract"""
    fq1_cmd = Pipeline(
        Command(
            unzip,
            "-dc",
            *[str(x) for x in fastq1],
        )
    )
    fq2_cmd = Pipeline(
        Command(
            unzip,
            "-dc",
            *[str(x) for x in fastq2],
        )
    )

    bwa_cmd = Command(
        "sentieon",
        "bwa",
        "mem",
        "-R",
        readgroup,
        "-t",
        str(threads),
        "-x",
        str(bwa_model),
        "-K",
        "10000000",
        str(reference),
        InputProcSub(fq1_cmd),
        InputProcSub(fq2_cmd),
    )

    extract_cmd = Command(
        "sentieon",
        "pgutil",
        "extract",
        "-t",
        str(threads),
        "-S",
        "-PG",
        "-m",
        str(extract_model),
        "-O",
        "fastq_compression=1",
        "-b",
        str(out_bam),
        "-o",
        str(out_fastq),
    )

    return Pipeline(bwa_cmd, extract_cmd)


def cmd_vg_convert_gfa(
    output_gfa: pathlib.Path,
    input_gbz: pathlib.Path,
    reference_name="GRCh38",
    threads=1,
) -> Pipeline:
    """Convert GBZ to GFA"""
    cmd = [
        "vg",
        "convert",
        "-t",
        str(threads),
        "-f",
        "-Q",
        reference_name,
        str(input_gbz),
    ]
    return Pipeline(Command(*cmd), file_output=output_gfa)


def cmd_vg_paths_fasta(
    output_fasta: pathlib.Path,
    input_gbz: pathlib.Path,
) -> Pipeline:
    """Extract haplotype sequences as FASTA"""
    cmd = [
        "vg",
        "paths",
        "-x",
        str(input_gbz),
        "-H",
        "-F",
    ]
    return Pipeline(Command(*cmd), file_output=output_fasta)


def cmd_minimap2_lift(
    out_bam: pathlib.Path,
    sample_fasta: pathlib.Path,
    input_fastq: pathlib.Path,
    gfa_file: pathlib.Path,
    reference_fasta: pathlib.Path,
    readgroup: str,
    minimap2_model: Union[pathlib.Path, str],
    minimap2_i: str = "16G",
    threads=1,
    mm2_xargs: List[str] = [],
) -> Pipeline:
    """Minimap2 alignment with pgutil lift"""
    mm2_cmd = Command(
        "sentieon",
        "minimap2",
        "-t",
        str(threads),
        "-R",
        readgroup,
        "-a",
        "-y",
        "-I",
        minimap2_i,
        "-x",
        str(minimap2_model),
        *mm2_xargs,
        str(sample_fasta),
        str(input_fastq),
    )

    lift_cmd = Command(
        "sentieon",
        "pgutil",
        "lift",
        "-t",
        str(threads),
        "-F",
        str(reference_fasta) + ".fai",
        "-g",
        str(gfa_file),
    )

    sort_cmd = Command(
        "sentieon",
        "util",
        "sort",
        "-i",
        "-",
        "-r",
        str(reference_fasta),
        "-t",
        str(threads),
        "-o",
        str(out_bam),
    )

    return Pipeline(mm2_cmd, lift_cmd, sort_cmd)


def cmd_bcftools_merge_trim(
    shard_vcf: pathlib.Path,
    raw_vcf: pathlib.Path,
    pop_vcf: pathlib.Path,
    trim_script: pathlib.Path,
    shard: str,
    merge_rules: str,
    merge_xargs: List[str] = [],
    view_xargs: List[str] = [],
):
    """Merge and trim variants"""
    merge_cmd = Command(
        "bcftools",
        "merge",
        "-r",
        shard,
        *merge_xargs,
        "-i",
        merge_rules,
        str(raw_vcf),
        str(pop_vcf),
    )

    trim_cmd = Command(
        "python3",
        str(trim_script),
    )

    view_cmd = Command(
        "bcftools",
        "view",
        *view_xargs,
        "-W=tbi",
        "-o",
        str(shard_vcf),
    )

    return Pipeline(merge_cmd, trim_cmd, view_cmd)


def cmd_bcftools_view_regions(
    out_vcf: pathlib.Path,
    in_vcf: pathlib.Path,
    regions: Optional[str] = None,
    regions_file: Optional[pathlib.Path] = None,
):
    xargs = []
    if regions:
        xargs.extend(["--regions", regions])
    if regions_file:
        xargs.extend(["--regions-file", str(regions_file)])
    return Pipeline(
        Command(
            "bcftools",
            "view",
            "--no-version",
            "-W=tbi",
            "-O",
            "z",
            "-o",
            str(out_vcf),
            *xargs,
            str(in_vcf),
        )
    )
