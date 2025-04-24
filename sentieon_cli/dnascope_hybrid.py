"""
Functionality for the DNAscope hybrid pipeline
"""

import argparse
from contextlib import contextmanager
import multiprocessing as mp
import os
import json
import pathlib
import shlex
import shutil
import sys
from typing import Any, List, Optional, Set, Tuple, Union

import packaging.version

from argh import arg, CommandError
from importlib_resources import files

from .archive import ar_load
from . import command_strings as cmds
from .dag import DAG
from .driver import (
    BaseDriver,
    Driver,
    DNAscope,
    DNAModelApply,
    HybridStage1,
    HybridStage2,
    HybridStage3,
)
from .dnascope import (
    align_fastq,
    call_cnvs,
    dedup_and_metrics,
    multiqc,
    set_bwt_max_mem,
)
from .dnascope_longread import align_inputs, call_svs, mosdepth
from .executor import DryRunExecutor, LocalExecutor
from .job import Job
from .logging import get_logger
from .scheduler import ThreadScheduler
from .util import (
    __version__,
    check_version,
    library_preloaded,
    parse_rg_line,
    path_arg,
    spit_alignment,
    tmp,
)

logger = get_logger(__name__)


CALLING_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202503"),
    "bedtools": None,
    "bcftools": packaging.version.Version("1.10"),
    "samtools": packaging.version.Version("1.16"),
}


@contextmanager
def cd(newdir: Union[str, pathlib.Path]):
    """from: https://stackoverflow.com/a/24176022/3203719"""
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def call_variants(
    tmp_dir: pathlib.Path,
    output_vcf: pathlib.Path,
    reference: pathlib.Path,
    sr_aln: List[pathlib.Path],
    lr_aln: List[pathlib.Path],
    model_bundle: pathlib.Path,
    hybrid_rg_sm: str = "",
    dbsnp: Optional[pathlib.Path] = None,
    bed: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),
    gvcf: bool = False,
    longread_tech: str = "HiFi",
    shortread_tech: str = "Illumina",
    dry_run: bool = False,
    skip_version_check: bool = False,
    stage3_ext: int = 1000,
    **_kwargs: Any,
) -> Tuple[
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
    Job,
]:
    """
    Call SNVs and indels using the DNAscope hybrid pipeline
    """

    ref_fai = pathlib.Path(str(reference) + ".fai")

    if not skip_version_check:
        for cmd, min_version in CALLING_MIN_VERSIONS.items():
            if not check_version(cmd, min_version):
                sys.exit(2)

    # readgroup handling for long-reads
    tr_read_filter: List[str] = []
    replace_rg_args: List[List[str]] = []
    for aln in lr_aln:
        aln_rgs = cmds.get_rg_lines(aln, dry_run)
        replace_rg_args.append([])
        for rg_line in aln_rgs:
            id = parse_rg_line(rg_line).get("ID")
            if not id:
                logger.error(
                    "Input file '%s' has a readgroup without an ID tag: %s",
                    aln,
                    rg_line,
                )
                sys.exit(2)
            replace_rg_args[-1].append(
                f"{id}={rg_line}\tLR:1".replace("\t", "\\t")
            )
            if longread_tech.upper() == "ONT":
                tr_read_filter.append(
                    f"TrimRepeat,max_repeat_unit=2,min_repeat_span=6,rgid={id}"
                )

    # readgroup handling for short-reads
    ultima_read_filter: List[str] = []
    if shortread_tech.upper() == "ULTIMA":
        for aln in sr_aln:
            aln_rgs = cmds.get_rg_lines(aln, dry_run)
            for rg_line in aln_rgs:
                id = parse_rg_line(rg_line).get("ID")
                if not id:
                    logger.error(
                        "Input file '%s' has a readgroup without an ID tag:%s",
                        aln,
                        rg_line,
                    )
                    sys.exit(2)
                ultima_read_filter.append(f"UltimaReadFilter,rgid={id}")

    # First pass - combined variant calling
    vcf_suffix = ".g.vcf.gz" if gvcf else ".vcf.gz"
    combined_vcf = tmp_dir.joinpath("initial" + vcf_suffix)
    emit_mode = "gvcf" if gvcf else None
    driver: BaseDriver = Driver(
        reference=reference,
        thread_count=cores,
        replace_rg=replace_rg_args,
        input=lr_aln + sr_aln,
        interval=bed,
        read_filter=tr_read_filter + ultima_read_filter,
    )
    driver.add_algo(
        DNAscope(
            combined_vcf,
            dbsnp=dbsnp,
            model=model_bundle.joinpath("hybrid.model"),
            pcr_indel_model="none",
            emit_mode=emit_mode,
        )
    )
    call_job = Job(shlex.join(driver.build_cmd()), "calling-1", cores)

    # Region selection
    hybrid_select = pathlib.Path(
        str(files("sentieon_cli.scripts").joinpath("hybrid_select.py"))
    ).resolve()
    selected_bed = tmp_dir.joinpath("selected.bed")
    select_job = Job(
        cmds.cmd_pyexec_hybrid_select(
            out_bed=selected_bed,
            vcf=combined_vcf,
            ref_fai=ref_fai,
            hybrid_select=hybrid_select,
            threads=cores,
        ),
        "hybrid-select",
        0,
    )

    mapq0_bed = tmp_dir.joinpath("hybrid_mapq0.bed")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        replace_rg=replace_rg_args,
        input=lr_aln + sr_aln,
    )
    driver.add_algo(
        HybridStage2(
            mode="4",
            all_bed=mapq0_bed,
        )
    )
    mapq0_job = Job(
        shlex.join(driver.build_cmd()),
        "mapq0-bed",
        cores,
    )

    mapq0_slop_bed = tmp_dir.joinpath("hybrid_mapq0.ex1000.bed")
    mapq0_slop_job = Job(
        cmds.cmd_bedtools_slop(
            mapq0_bed,
            mapq0_slop_bed,
            ref_fai,
            1000,
        ),
        "mapq0-bed-slop",
        0,
    )

    diff_bed = tmp_dir.joinpath("merged_diff.bed")
    cat_merge_job = Job(
        cmds.cmd_bedtools_cat_sort_merge(
            out_bed=diff_bed,
            in_bed=[selected_bed, mapq0_slop_bed],
            ref_fai=ref_fai,
        ),
        "concat-merge-bed",
        0,
    )
    rm_cmd = ["rm", str(selected_bed), str(mapq0_slop_bed)]
    rm_job1 = Job(shlex.join(rm_cmd), "rm-tmp1", 0, True)

    stage1_ins_fa = tmp_dir.joinpath("stage1_ins.fa")
    stage1_ins_bed = tmp_dir.joinpath("stage1_ins.bed")
    ins_driver = Driver(
        reference=reference,
        thread_count=cores,
        replace_rg=replace_rg_args,
        input=lr_aln,
    )
    ins_driver.add_algo(
        HybridStage1(
            "-",
            model=model_bundle.joinpath("HybirdStage1_ins.model"),
            fa_file=stage1_ins_fa,
            bed_file=stage1_ins_bed,
        )
    )

    stage1_hap_bam = tmp_dir.joinpath("stage1_hap.bam")
    stage1_hap_bed = tmp_dir.joinpath("stage1_hap.bed")
    stage1_hap_vcf = tmp_dir.joinpath("stage1_hap.vcf")
    stage1_driver = Driver(
        reference=reference,
        thread_count=cores,
        replace_rg=replace_rg_args,
        input=lr_aln,
        interval=diff_bed,
    )
    stage1_driver.add_algo(
        HybridStage1(
            "-",
            model=model_bundle.joinpath("HybridStage1.model"),
            hap_bam=stage1_hap_bam,
            hap_bed=stage1_hap_bed,
            hap_vcf=stage1_hap_vcf,
        )
    )

    stage1_bam = tmp_dir.joinpath("hybrid_stage1.bam")
    stage1_job = Job(
        cmds.hybrid_stage1(
            out_aln=stage1_bam,
            reference=reference,
            cores=cores,
            readgroup=f"@RG\\tID:hybrid-18893\\tSM:{hybrid_rg_sm}",
            ins_driver=ins_driver,
            stage1_driver=stage1_driver,
        ),
        "first-stage",
        cores,
    )
    rm_cmd = [
        "rm",
        str(stage1_ins_fa),
        str(stage1_ins_bed),
        str(stage1_hap_vcf),
    ]
    rm_job2 = Job(shlex.join(rm_cmd), "rm-tmp2", 0, True)

    stage2_bed = tmp_dir.joinpath("hybrid_stage2.bed")
    stage2_reg_bed = tmp_dir.joinpath("hybrid_stage2_reg.bed")
    stage2_unmap_bam = tmp_dir.joinpath("hybrid_stage2_unmap.bam")
    stage2_alt_bam = tmp_dir.joinpath("hybrid_stage2_alt.bam")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        input=[stage1_bam, stage1_hap_bam],
    )
    driver.add_algo(
        HybridStage2(
            model=model_bundle.joinpath("HybridStage2.model"),
            hap_bed=stage1_hap_bed,
            unmap_bam=stage2_unmap_bam,
            alt_bam=stage2_alt_bam,
            all_bed=stage2_bed,
        )
    )
    second_stage_job = Job(
        shlex.join(driver.build_cmd()),
        "second-stage",
        cores,
    )
    merge_job = Job(
        cmds.cmd_bedtools_merge(
            stage2_bed,
            stage2_reg_bed,
            distance=2 * stage3_ext,
        ),
        "merge-bed",
        0,
    )
    rm_cmd = ["rm", str(stage1_bam), str(stage1_hap_bam), str(stage2_bed)]
    rm_job3 = Job(shlex.join(rm_cmd), "rm-tmp3", 0, True)

    stage3_bam = tmp_dir.joinpath("hybrid_stage3.bam")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        replace_rg=replace_rg_args,
        input=lr_aln + [stage2_unmap_bam, stage2_alt_bam] + sr_aln,
        interval=stage2_reg_bed,
    )
    driver.add_algo(
        HybridStage3(
            "-",
            model=model_bundle.joinpath("HybridStage3.model"),
            region_ext=stage3_ext,
        )
    )
    third_stage_job = Job(
        cmds.hybrid_stage3(
            stage3_bam,
            driver=driver,
            cores=cores,
        ),
        "third-stage",
        cores,
    )
    rm_cmd = ["rm", str(stage2_unmap_bam), str(stage2_alt_bam)]
    rm_job4 = Job(shlex.join(rm_cmd), "rm-tmp4", 0, True)

    # pass 2 of variant calling
    pass2_vcf = tmp_dir.joinpath("hybrid_pass2.vcf.gz")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        replace_rg=replace_rg_args,
        input=lr_aln + [stage3_bam],
        interval=stage2_reg_bed,
        read_filter=tr_read_filter + ultima_read_filter,
    )
    driver.add_algo(
        DNAscope(
            pass2_vcf,
            dbsnp=dbsnp,
            model=model_bundle.joinpath("hybrid.model"),
            pcr_indel_model="none",
            emit_mode=emit_mode,
        )
    )
    call2_job = Job(shlex.join(driver.build_cmd()), "call2", cores)

    # Merge and normalize the VCFs
    subset_vcf = tmp_dir.joinpath("mix_subset.vcf.gz")
    combined_tmp_vcf = tmp_dir.joinpath("combined_tmp.vcf.gz")
    subset_job = Job(
        cmds.bcftools_subset(
            subset_vcf,
            combined_vcf,
            stage2_reg_bed,
        ),
        "subset-calls",
        0,
    )
    concat_job = Job(
        cmds.bcftools_concat(
            combined_tmp_vcf,
            [subset_vcf, pass2_vcf],
        ),
        "concat-calls",
        0,
    )
    rm_cmd = ["rm", str(combined_vcf), str(subset_vcf), str(pass2_vcf)]
    rm_job5 = Job(shlex.join(rm_cmd), "rm-tmp5", 0, True)

    # Annotate the output VCF
    hybrid_anno = pathlib.Path(
        str(files("sentieon_cli.scripts").joinpath("hybrid_anno.py"))
    )
    combined_anno_vcf = tmp_dir.joinpath("combined_tmp_anno.vcf.gz")
    anno_job = Job(
        cmds.cmd_pyexec_hybrid_anno(
            combined_anno_vcf,
            combined_tmp_vcf,
            stage1_hap_bed,
            hybrid_anno,
            cores,
        ),
        "anno-calls",
        0,
    )

    # Model Apply
    apply_vcf = tmp_dir.joinpath("combined_apply.vcf.gz")
    driver = Driver(
        reference=reference,
        thread_count=cores,
        interval=bed,
    )
    driver.add_algo(
        DNAModelApply(
            model=model_bundle.joinpath("hybrid.model"),
            vcf=combined_anno_vcf,
            output=apply_vcf,
        )
    )
    apply_job = Job(shlex.join(driver.build_cmd()), "model-apply", cores)

    # Final normalize
    norm_job = Job(
        cmds.filter_norm(
            output_vcf,
            apply_vcf,
            reference,
            exclude_homref=not gvcf,
        ),
        "final-norm",
        0,
    )
    return (
        call_job,
        select_job,
        mapq0_job,
        mapq0_slop_job,
        cat_merge_job,
        rm_job1,
        stage1_job,
        rm_job2,
        second_stage_job,
        merge_job,
        rm_job3,
        third_stage_job,
        rm_job4,
        call2_job,
        subset_job,
        concat_job,
        rm_job5,
        anno_job,
        apply_job,
        norm_job,
    )


@arg(
    "-r",
    "--reference",
    required=True,
    help="fasta for reference genome",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--sr_r1_fastq",
    nargs="*",
    help="Short-read R1 fastq files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--sr_r2_fastq",
    nargs="*",
    help="Short-read R2 fastq files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--sr_readgroups",
    nargs="*",
    help="Readgroup information for the short-read fastq files",
)
@arg(
    "--sr_aln",
    nargs="*",
    help="Short-read BAM or CRAM files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--sr_duplicate_marking",
    help="Options for duplicate marking.",
    choices=["markdup", "rmdup", "none"],
)
@arg(
    "--lr_aln",
    nargs="*",
    help="Long-read BAM or CRAM files",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "-m",
    "--model_bundle",
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
    "--longread_tech",
    help="Sequencing technology used to generate the long reads.",
    choices=["HiFi", "ONT"],
)
@arg(
    "--dry_run",
    help="Print the commands without running them.",
)
@arg(
    "--skip_svs",
    help="Skip SV calling",
)
@arg(
    "--skip_metrics",
    help="Skip all metrics collection and multiQC",
)
@arg(
    "--skip_mosdepth",
    help="Skip QC with mosdepth",
)
@arg(
    "--skip_multiqc",
    help="Skip multiQC report generation",
)
@arg(
    "--skip_cnv",
    help="Skip CNV calling",
)
@arg(
    "--lr_align_input",
    help="Align the input long-read BAM/CRAM/uBAM file to the reference"
    " genome",
    action="store_true",
)
@arg(
    "--lr_input_ref",
    help="Used to decode the input long-read alignment file. Required if the "
    "input file is in the CRAM/uCRAM formats",
    type=path_arg(exists=True, is_file=True),
)
@arg(
    "--lr_fastq_taglist",
    help="A comma-separated list of tags to retain. Defaults to '%(default)s'"
    " and the 'RG' tag is required",
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
    "--bwa_k",
    help="The '-K' argument in bwa",
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
    "--skip_version_check",
    help=argparse.SUPPRESS,
    action="store_true",
)
@arg(
    "--retain_tmpdir",
    help=argparse.SUPPRESS,
    action="store_true",
)
# Manually set `bwt_max_mem`
@arg(
    "--bwt_max_mem",
    help=argparse.SUPPRESS,
)
# Do not split alignment into multiple jobs
@arg(
    "--no_split_alignment",
    help=argparse.SUPPRESS,
    action="store_true",
)
def dnascope_hybrid(
    output_vcf: pathlib.Path,
    sr_r1_fastq: Optional[List[pathlib.Path]] = None,
    sr_r2_fastq: Optional[List[pathlib.Path]] = None,
    sr_readgroups: Optional[List[str]] = None,
    sr_aln: Optional[List[pathlib.Path]] = None,
    lr_aln: Optional[List[pathlib.Path]] = None,
    reference: Optional[pathlib.Path] = None,
    model_bundle: Optional[pathlib.Path] = None,
    dbsnp: Optional[pathlib.Path] = None,
    bed: Optional[pathlib.Path] = None,
    cores: int = mp.cpu_count(),
    gvcf: bool = False,  # pylint: disable=W0613
    sr_duplicate_marking: str = "markdup",
    longread_tech: str = "HiFi",  # pylint: disable=W0613
    dry_run: bool = False,
    skip_svs: bool = False,
    skip_mosdepth: bool = False,
    skip_metrics: bool = False,
    skip_multiqc: bool = False,
    skip_cnv: bool = False,
    lr_align_input: bool = False,
    lr_input_ref: Optional[pathlib.Path] = None,
    lr_fastq_taglist: str = "*",  # pylint: disable=W0613
    bam_format: bool = False,  # pylint: disable=W0613
    bwa_args: str = "",  # pylint: disable=W0613
    bwa_k: int = 100000000,
    minimap2_args: str = "-Y",  # pylint: disable=W0613
    util_sort_args: str = (
        "--cram_write_options version=3.0,compressor=rans"
    ),  # pylint: disable=W0613
    skip_version_check: bool = False,  # pylint: disable=W0613
    retain_tmpdir: bool = False,
    bwt_max_mem: Optional[str] = None,
    no_split_alignment: bool = False,
    **kwargs: str,
):
    """
    Run the DNAscope hybrid pipeline
    """
    assert reference
    assert model_bundle
    assert str(output_vcf).endswith(".vcf.gz")
    assert sr_aln or (sr_r1_fastq and sr_readgroups)
    sr_r1_fastq = sr_r1_fastq if sr_r1_fastq else []
    sr_r2_fastq = sr_r2_fastq if sr_r2_fastq else []
    sr_readgroups = sr_readgroups if sr_readgroups else []
    sr_aln = sr_aln if sr_aln else []
    lr_aln = lr_aln if lr_aln else []
    assert lr_aln
    bwa_k_arg = str(bwa_k)
    skip_multiqc = True if skip_metrics else skip_multiqc
    skip_mosdepth = True if skip_metrics else skip_mosdepth

    assert logger.parent
    logger.parent.setLevel(kwargs["loglevel"])
    logger.info("Starting sentieon-cli version: %s", __version__)

    n_alignment_jobs = 1
    numa_nodes: List[str] = []
    if not no_split_alignment:
        numa_nodes = spit_alignment(cores)
    n_alignment_jobs = max(1, len(numa_nodes))
    set_bwt_max_mem(bwt_max_mem, None, n_alignment_jobs)

    # Parse the bundle_info.json file
    longread_tech_cli = longread_tech
    bundle_info_bytes = ar_load(str(model_bundle) + "/bundle_info.json")
    if isinstance(bundle_info_bytes, list):
        bundle_info_bytes = b"{}"

    bundle_info = json.loads(bundle_info_bytes.decode())
    longread_tech = bundle_info.get("longReadPlatform")
    shortread_tech = bundle_info.get("shortReadPlatform")
    if not longread_tech or not shortread_tech:
        # logger.warning(
        #     "The model bundle file does not have the expected attributes. "
        #     "Are you using the latest version?"
        # )
        pass
    if not longread_tech:
        longread_tech = longread_tech_cli
    if not shortread_tech:
        shortread_tech = "Illumina"
    req_version = packaging.version.Version(
        bundle_info.get("minScriptVersion", __version__)
    )
    if req_version > packaging.version.Version(__version__):
        logger.error(
            "The model bundle requires version %s or later of the "
            "sentieon-cli.",
            req_version,
        )
        sys.exit(2)
    if bundle_info.get("pipeline", "") != "DNAscope Hybrid":
        logger.error("The model bundle is for a different pipeline.")
        sys.exit(2)

    # Confirm that all readgroups have the same RGSM
    hybrid_rg_sm = ""
    rg_sm_tag = None
    for aln_files in sr_aln, lr_aln:
        for aln in aln_files:
            rg_lines = cmds.get_rg_lines(aln, dry_run)
            for rg_line in rg_lines:
                sm = parse_rg_line(rg_line).get("SM")
                if not sm:
                    logger.error(
                        "Input file '%s' has a readgroup line without a SM"
                        " tag: %s",
                        aln,
                        rg_line,
                    )
                    sys.exit(2)
                if rg_sm_tag is None:
                    rg_sm_tag = sm
                    hybrid_rg_sm = sm
                    continue
                if dry_run:
                    continue
                elif rg_sm_tag != sm:
                    logger.error(
                        "Input file '%s' has a different RG-SM tag from"
                        " previously seen alignment files.\nfound='%s'"
                        " expected='%s'",
                        aln,
                        sm,
                        rg_sm_tag,
                    )
                    sys.exit(2)

    for rg in sr_readgroups:
        sm = parse_rg_line(rg.replace(r"\t", "\t")).get("SM")
        if not sm:
            logger.error("Input readgroup does not have a SM tag: %s", rg)
            sys.exit(2)
        if rg_sm_tag is None:
            rg_sm_tag = sm
            hybrid_rg_sm = sm
            continue
        if dry_run:
            continue
        elif rg_sm_tag != sm:
            logger.error(
                "Input readgroup, '%s' has a different RG-SM tag from"
                " previously seen tags.\nfound='%s'\nexpected='%s'",
                rg,
                sm,
                rg_sm_tag,
            )
            sys.exit(2)

    if not library_preloaded("libjemalloc.so"):
        logger.warning(
            "jemalloc is recommended, but is not preloaded. See "
            "https://support.sentieon.com/appnotes/jemalloc/"
        )

    if bed is None:
        logger.info(
            "A BED file is recommended to avoid variant calling across "
            "decoy and unplaced contigs."
        )

    tmp_dir_str = tmp()
    tmp_dir = pathlib.Path(tmp_dir_str)

    logger.info("Building the DAG")
    dag = DAG()

    # Short-read alignment
    aligned_fastq, align_fastq_jobs, fq_rm_job = align_fastq(
        duplicate_marking=sr_duplicate_marking,
        r1_fastq=sr_r1_fastq,
        r2_fastq=sr_r2_fastq,
        readgroups=sr_readgroups,
        **locals(),
    )
    for job in align_fastq_jobs:
        dag.add_job(job)
    sr_aln.extend(aligned_fastq)

    # Short-read dedup
    deduped, lc_job, dedup_job, metrics_job, rehead_job = dedup_and_metrics(
        sample_input=sr_aln,
        duplicate_marking=sr_duplicate_marking,
        assay="WGS",
        **locals(),
    )
    sr_aln = deduped
    if lc_job:
        dag.add_job(lc_job, align_fastq_jobs)
        if dedup_job:
            dag.add_job(dedup_job, {lc_job})
            if metrics_job and not skip_metrics:
                dag.add_job(metrics_job, {dedup_job})
                if rehead_job:
                    dag.add_job(rehead_job, {metrics_job})
            if fq_rm_job:
                dag.add_job(fq_rm_job, {dedup_job})

    # Long-read alignment
    realign_jobs: Set[Job] = set()
    if lr_align_input:
        lr_aln, realign_jobs = align_inputs(
            sample_input=lr_aln,
            fastq_taglist=lr_fastq_taglist,
            input_ref=lr_input_ref,
            **locals(),
        )
        for job in realign_jobs:
            dag.add_job(job)

    if not skip_mosdepth:
        mosdpeth_jobs = mosdepth(sample_input=lr_aln, **locals())
        for job in mosdpeth_jobs:
            dag.add_job(job, realign_jobs)

    if not skip_multiqc:
        multiqc_job = multiqc(**locals())
        multiqc_dependencies: Set[Job] = set()
        if lc_job:
            multiqc_dependencies.add(lc_job)
        if metrics_job:
            multiqc_dependencies.add(metrics_job)
        if rehead_job:
            multiqc_dependencies.add(rehead_job)

        if multiqc_job:
            dag.add_job(multiqc_job, multiqc_dependencies)

    if not skip_svs:
        sample_input = lr_aln
        longreadsv_job = call_svs(**locals())
        dag.add_job(longreadsv_job, realign_jobs)

    sr_preprocessing_jobs: Set[Job] = set()
    sr_preprocessing_jobs.update(align_fastq_jobs)
    if dedup_job:
        sr_preprocessing_jobs.add(dedup_job)
    if not skip_cnv:
        sample_input = sr_aln
        cnvscope_job, cnvmodelapply_job = call_cnvs(**locals())
        dag.add_job(cnvscope_job, sr_preprocessing_jobs)
        dag.add_job(cnvmodelapply_job, {cnvscope_job})

    (
        call_job,
        select_job,
        mapq0_job,
        mapq0_slop_job,
        cat_merge_job,
        rm_job1,
        stage1_job,
        rm_job2,
        second_stage_job,
        merge_job,
        rm_job3,
        third_stage_job,
        rm_job4,
        call2_job,
        subset_job,
        concat_job,
        rm_job5,
        anno_job,
        apply_job,
        norm_job,
    ) = call_variants(**locals())
    dag.add_job(call_job, realign_jobs | sr_preprocessing_jobs)
    dag.add_job(select_job, {call_job})
    dag.add_job(mapq0_job, realign_jobs | sr_preprocessing_jobs)
    dag.add_job(mapq0_slop_job, {mapq0_job})
    dag.add_job(cat_merge_job, {mapq0_slop_job, select_job})
    dag.add_job(stage1_job, {cat_merge_job})
    dag.add_job(second_stage_job, {stage1_job})
    dag.add_job(merge_job, {second_stage_job})
    dag.add_job(third_stage_job, {merge_job})
    dag.add_job(call2_job, {third_stage_job})
    dag.add_job(subset_job, {merge_job})
    dag.add_job(concat_job, {subset_job, call2_job})
    dag.add_job(anno_job, {concat_job})
    dag.add_job(apply_job, {anno_job})
    dag.add_job(norm_job, {apply_job})

    # Remove intermediate files during processing
    if not retain_tmpdir:
        dag.add_job(rm_job1, {cat_merge_job})
        dag.add_job(rm_job2, {stage1_job})
        dag.add_job(rm_job3, {merge_job})
        dag.add_job(rm_job4, {third_stage_job})
        dag.add_job(rm_job5, {concat_job})

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
