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
from .job import Job
from .pipeline import BasePipeline
from .util import (
    check_version,
    path_arg,
    library_preloaded,
)

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


class DNAscopeLRPipeline(BasePipeline):
    """The DNAscope LongRead pipeline"""

    params: Dict[str, Dict[str, Any]] = {
        "fastq": {
            "nargs": "*",
            "help": "Sample fastq files.",
            "type": path_arg(exists=True, is_file=True),
        },
        "model_bundle": {
            "flags": ["-m", "--model_bundle"],
            "help": "The model bundle file.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "readgroups": {
            "nargs": "*",
            "help": "Readgroup information for the fastq files.",
        },
        "reference": {
            "flags": ["-r", "--reference"],
            "required": True,
            "help": "fasta for reference genome.",
            "type": path_arg(exists=True, is_file=True),
        },
        "sample_input": {
            "flags": ["-i", "--sample_input"],
            "nargs": "*",
            "help": "sample BAM or CRAM file.",
            "type": path_arg(exists=True, is_file=True),
        },
        "align": {
            "help": (
                "Align the input BAM/CRAM/uBAM file to the reference genome."
            ),
            "action": "store_true",
        },
        "bam_format": {
            "help": (
                "Use the BAM format instead of CRAM for output aligned files."
            ),
            "action": "store_true",
        },
        "bed": {
            "flags": ["-b", "--bed"],
            "help": (
                "Region BED file. Supplying this file will restrict diploid "
                "variant calling to the intervals inside the BED file."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "cnv_excluded_regions": {
            "help": "Regions to exclude from CNV calling.",
            "type": path_arg(exists=True, is_file=True),
        },
        "cores": {
            "flags": ["-t", "--cores"],
            "help": (
                "Number of threads/processes to use. Defaults to all "
                "available."
            ),
            "default": mp.cpu_count(),
        },
        "dbsnp": {
            "flags": ["-d", "--dbsnp"],
            "help": (
                "dbSNP vcf file Supplying this file will annotate variants "
                "with their dbSNP refSNP ID numbers."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "dry_run": {
            "help": "Print the commands without running them.",
            "action": "store_true",
        },
        "gvcf": {
            "flags": ["-g", "--gvcf"],
            "help": (
                "Generate a gVCF output file along with the VCF."
                " (default generates only the VCF)"
            ),
            "action": "store_true",
        },
        "haploid_bed": {
            "help": (
                "A BED file of haploid regions. Supplying this file will "
                "perform haploid variant calling across these regions."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "input_ref": {
            "help": (
                "Used to decode the input alignment file. Required if the "
                "input file is in the CRAM/uCRAM formats."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "skip_cnv": {
            "help": "Skip CNV calling.",
            "action": "store_true",
        },
        "skip_mosdepth": {
            "help": "Skip QC with mosdepth.",
            "action": "store_true",
        },
        "skip_small_variants": {
            "help": "Skip small variant (SNV/indel) calling.",
            "action": "store_true",
        },
        "skip_svs": {
            "help": "Skip SV calling.",
            "action": "store_true",
        },
        "tech": {
            "help": "Sequencing technology used to generate the reads.",
            "choices": ["HiFi", "ONT"],
            "default": "HiFi",
        },
        "fastq_taglist": {
            # "help": (
            #    "A comma-separated list of tags to retain. Defaults to "
            #    "'%(default)s' and the 'RG' tag is required."
            # ),
            "help": argparse.SUPPRESS,
            "default": "*",
        },
        "minimap2_args": {
            # "help": "Extra arguments for sentieon minimap2.",
            "help": argparse.SUPPRESS,
            "default": "-Y",
        },
        "repeat_model": {
            "help": argparse.SUPPRESS,
            "type": path_arg(exists=True, is_file=True),
        },
        "retain_tmpdir": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
        },
        "skip_version_check": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
        },
        "use_pbsv": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
        },
        "util_sort_args": {
            # "help": "Extra arguments for sentieon util sort.",
            "help": argparse.SUPPRESS,
            "default": "--cram_write_options version=3.0,compressor=rans",
        },
    }

    positionals: Dict[str, Dict[str, Any]] = {
        "output_vcf": {
            "help": "Output VCF File. The file name must end in .vcf.gz",
            "type": path_arg(),
        },
    }

    def __init__(self) -> None:
        super().__init__()
        self.output_vcf: Optional[pathlib.Path] = None
        self.reference: Optional[pathlib.Path] = None
        self.sample_input: List[pathlib.Path] = []
        self.fastq: List[pathlib.Path] = []
        self.readgroups: List[str] = []
        self.model_bundle: Optional[pathlib.Path] = None
        self.dbsnp: Optional[pathlib.Path] = None
        self.bed: Optional[pathlib.Path] = None
        self.haploid_bed: Optional[pathlib.Path] = None
        self.cnv_excluded_regions: Optional[pathlib.Path] = None
        self.gvcf = False
        self.tech = "HiFi"
        self.dry_run = False
        self.skip_small_variants = False
        self.skip_svs = False
        self.skip_mosdepth = False
        self.skip_cnv = False
        self.align = False
        self.input_ref: Optional[pathlib.Path] = None
        self.fastq_taglist = "*"
        self.bam_format = False
        self.minimap2_args = "-Y"
        self.util_sort_args = (
            "--cram_write_options version=3.0,compressor=rans"
        )
        self.repeat_model: Optional[pathlib.Path] = None
        self.skip_version_check = False
        self.retain_tmpdir = False
        self.use_pbsv = False

    def validate(self) -> None:
        # uniquify pipeline attributes
        self.lr_aln = self.sample_input
        del self.sample_input

        # validate
        if not (self.lr_aln or self.fastq):
            self.logger.error(
                "Please suppy either the `--sample_input` or `--fastq` and "
                "`--readgroups` arguments"
            )
            sys.exit(2)
        if len(self.fastq) != len(self.readgroups):
            self.logger.error(
                "The number of readgroups does not equal the number of fastq "
                "files"
            )
            sys.exit(2)
        if not str(self.output_vcf).endswith(".vcf.gz"):
            self.logger.error("The output file should end with '.vcf.gz'")
            sys.exit(2)
        if self.tech.upper() == "ONT":
            self.logger.info("Skipping CNV calling with ONT data")
            self.skip_cnv = True

        if not library_preloaded("libjemalloc.so"):
            self.logger.warning(
                "jemalloc is recommended, but is not preloaded. See "
                "https://support.sentieon.com/appnotes/jemalloc/"
            )

        if not self.cnv_excluded_regions and not self.skip_cnv:
            self.logger.warning(
                "Excluded regions are recommended for CNV calling. Please "
                "supply them with the `--cnv_excluded_regions` argument"
            )

        if not self.skip_small_variants:
            if self.haploid_bed and not self.bed:
                self.logger.error(
                    "Please supply a BED file of diploid regions to "
                    "distinguish haploid and diploid regions of the genome."
                )
                sys.exit(2)
            if not self.bed:
                self.logger.warning(
                    "A BED file is recommended to restrict variant calling to "
                    "diploid regions of the genome."
                )

    def configure(self) -> None:
        pass

    def build_dag(self) -> DAG:
        self.logger.info("Building the DAG")
        dag = DAG()

        sample_input = self.lr_aln
        realign_jobs: Set[Job] = set()
        if self.align:
            sample_input, realign_jobs = self.lr_align_inputs()
            for job in realign_jobs:
                dag.add_job(job)
        aligned_fastq, align_jobs = self.lr_align_fastq()
        sample_input.extend(aligned_fastq)
        for job in align_jobs:
            dag.add_job(job)

        if not self.skip_mosdepth:
            mosdpeth_jobs = self.mosdepth(sample_input)
            for job in mosdpeth_jobs:
                dag.add_job(job, realign_jobs.union(align_jobs))

        if self.use_pbsv or not self.skip_cnv:
            merged_bam, merge_job = self.merge_input_files(sample_input)
            dag.add_job(merge_job, realign_jobs.union(align_jobs))

            if self.use_pbsv:
                pbsv_discover, pbsv_call = self.pbsv(merged_bam)
                dag.add_job(pbsv_discover, {merge_job})
                dag.add_job(pbsv_call, {pbsv_discover})

            if not self.skip_cnv:
                hificnv_job = self.hificnv(merged_bam)
                if hificnv_job:
                    dag.add_job(hificnv_job, {merge_job})

        if not self.skip_small_variants:
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
            ) = self.lr_call_variants(sample_input)
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

            if (
                haploid_calling_job
                and haploid_patch2_job
                and haploid_concat_job
            ):
                # if haploid bed
                dag.add_job(haploid_calling_job, haploid_calling_deps)
                dag.add_job(haploid_patch2_job, {haploid_calling_job})
                dag.add_job(
                    haploid_concat_job, {haploid_patch2_job, merge_job}
                )

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

        if not self.skip_svs:
            longreadsv_job = self.call_svs(sample_input)
            dag.add_job(longreadsv_job, realign_jobs.union(align_jobs))

        return dag

    def lr_align_inputs(self) -> Tuple[List[pathlib.Path], Set[Job]]:
        """
        Align reads to the reference genome using minimap2
        """
        assert self.output_vcf
        assert self.reference
        assert self.model_bundle
        if not self.skip_version_check:
            for cmd, min_version in ALN_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        res: List[pathlib.Path] = []
        suffix = "bam" if self.bam_format else "cram"
        sample_name = self.output_vcf.name.replace(".vcf.gz", "")
        realign_jobs = set()
        for i, input_aln in enumerate(self.lr_aln):
            out_aln = pathlib.Path(
                str(self.output_vcf).replace(
                    ".vcf.gz", f"_mm2_sorted_{i}.{suffix}"
                )
            )
            rg_lines = cmds.get_rg_lines(
                input_aln,
                self.dry_run,
            )

            realign_jobs.add(
                Job(
                    cmds.cmd_samtools_fastq_minimap2(
                        out_aln,
                        input_aln,
                        self.reference,
                        self.model_bundle,
                        self.cores,
                        rg_lines,
                        sample_name,
                        self.input_ref,
                        self.fastq_taglist,
                        self.minimap2_args,
                        self.util_sort_args,
                    ),
                    f"bam-realign-{i}",
                    self.cores,
                )
            )
            res.append(out_aln)
        return (res, realign_jobs)

    def lr_align_fastq(self) -> Tuple[List[pathlib.Path], Set[Job]]:
        """
        Align fastq to the reference genome using minimap2
        """
        assert self.reference
        assert self.model_bundle
        res: List[pathlib.Path] = []
        if self.fastq is None and self.readgroups is None:
            return (res, set())

        if not self.skip_version_check:
            for cmd, min_version in FQ_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        unzip = "igzip"
        if not shutil.which(unzip):
            self.logger.warning(
                "igzip is recommended for decompression, but is not "
                "available. Falling back to gzip."
            )
            unzip = "gzip"

        suffix = "bam" if self.bam_format else "cram"
        align_jobs = set()
        for i, (fq, rg) in enumerate(zip(self.fastq, self.readgroups)):
            out_aln = pathlib.Path(
                str(self.output_vcf).replace(
                    ".vcf.gz", f"_mm2_sorted_fq_{i}.{suffix}"
                )
            )
            align_jobs.add(
                Job(
                    cmds.cmd_fastq_minimap2(
                        out_aln,
                        fq,
                        rg,
                        self.reference,
                        self.model_bundle,
                        self.cores,
                        unzip,
                        self.minimap2_args,
                        self.util_sort_args,
                    ),
                    "align-{i}",
                    self.cores,
                )
            )
            res.append(out_aln)
        return (res, align_jobs)

    def mosdepth(self, sample_input: List[pathlib.Path]) -> Set[Job]:
        """Run mosdepth for QC"""

        if not self.skip_version_check:
            if not all(
                [
                    check_version(cmd, min_version)
                    for (cmd, min_version) in MOSDEPTH_MIN_VERSIONS.items()
                ]
            ):
                self.logger.warning(
                    "Skipping mosdepth. mosdepth version %s or later "
                    "not found",
                    MOSDEPTH_MIN_VERSIONS["mosdepth"],
                )
                return set()

        mosdpeth_jobs = set()
        for i, input_file in enumerate(sample_input):
            mosdepth_dir = pathlib.Path(
                str(self.output_vcf).replace(".vcf.gz", f"_mosdepth_{i}")
            )
            mosdpeth_jobs.add(
                Job(
                    cmds.cmd_mosdepth(
                        input_file,
                        mosdepth_dir,
                        fasta=self.reference,
                        threads=self.cores,
                    ),
                    f"mosdepth-{i}",
                    0,  # Run in background
                )
            )
        return mosdpeth_jobs

    def merge_input_files(
        self, sample_input: List[pathlib.Path]
    ) -> Tuple[pathlib.Path, Job]:
        """
        Merge the aligned reads into a single BAM
        """
        if not self.skip_version_check:
            for cmd, min_version in MERGE_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        # Merge the sample_input into a single BAM
        merged_bam = self.tmp_dir.joinpath("longread_merged.bam")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
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
        self,
        merged_bam: pathlib.Path,
    ) -> Tuple[Job, Job]:
        """
        Call SVs with pbsv
        """
        if not self.skip_version_check:
            for cmd, min_version in PBSV_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        # pbsv discover
        pbsv_discover_fn = str(self.output_vcf).replace(
            ".vcf.gz",
            ".pbsv_discover.svsig.gz",
        )
        pbsv_discover = Job(
            f"pbsv discover --hifi {merged_bam} {pbsv_discover_fn}",
            "pbsv-discover",
            0,
        )

        # pbsv call
        pbsv_call_fn = str(self.output_vcf).replace(".vcf.gz", ".pbsv.vcf")
        cmd = (
            f"pbsv call -j {self.cores} {self.reference} {pbsv_discover_fn} "
            f"{pbsv_call_fn}"
        )
        pbsv_call = Job(
            cmd,
            "pbsv-call",
            self.cores,
        )
        return (pbsv_discover, pbsv_call)

    def hificnv(
        self,
        merged_bam: pathlib.Path,
    ) -> Optional[Job]:
        """
        Call CNVs with HiFiCNV
        """
        if not self.skip_version_check:
            for cmd, min_version in HIFICNV_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    self.logger.warning(
                        "Skipping hificnv. hificnv version %s or later not "
                        "found",
                        HIFICNV_MIN_VERSIONS["hificnv"],
                    )
                    return None

        hifi_cnv_fn = str(self.output_vcf).replace(".vcf.gz", ".hificnv")
        cmd = (
            f"hificnv --bam {merged_bam} --ref {self.reference} "
            f"--threads {self.cores} "
            f"--output-prefix {hifi_cnv_fn}"
        )
        if self.cnv_excluded_regions:
            cmd += f" --exclude {self.cnv_excluded_regions}"
        hificnv_job = Job(cmd, "hificnv", self.cores)
        return hificnv_job

    def lr_call_variants(
        self,
        sample_input: List[pathlib.Path],
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
        assert self.model_bundle
        assert self.reference
        assert self.output_vcf
        if not self.skip_version_check:
            for cmd, min_version in TOOL_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        # First pass - diploid calling
        diploid_gvcf_fn = self.tmp_dir.joinpath("out_diploid.g.vcf.gz")
        diploid_tmp_vcf = self.tmp_dir.joinpath("out_diploid_tmp.vcf.gz")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=sample_input,
            interval=self.bed,
        )
        if self.gvcf:
            driver.add_algo(
                DNAscope(
                    diploid_gvcf_fn,
                    dbsnp=self.dbsnp,
                    emit_mode="gvcf",
                    model=self.model_bundle.joinpath("gvcf_model"),
                )
            )
        driver.add_algo(
            DNAscope(
                diploid_tmp_vcf,
                dbsnp=self.dbsnp,
                model=self.model_bundle.joinpath("diploid_model"),
            )
        )
        first_calling_job = Job(
            shlex.join(driver.build_cmd()), "first-pass", self.cores
        )

        diploid_vcf = self.tmp_dir.joinpath("out_diploid.vcf.gz")
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            DNAModelApply(
                self.model_bundle.joinpath("diploid_model"),
                diploid_tmp_vcf,
                diploid_vcf,
            )
        )
        first_modelapply_job = Job(
            shlex.join(driver.build_cmd()), "first-modelapply", self.cores
        )

        # Phasing and RepeatModel
        phased_bed = self.tmp_dir.joinpath("out_diploid_phased.bed")
        unphased_bed = self.tmp_dir.joinpath("out_diploid_unphased.bed")
        phased_vcf = self.tmp_dir.joinpath("out_diploid_phased.vcf.gz")
        phased_ext = self.tmp_dir.joinpath("out_diploid_phased.ext.vcf.gz")
        phased_unphased = self.tmp_dir.joinpath(
            "out_diploid_phased_unphased.vcf.gz"
        )
        phased_phased = self.tmp_dir.joinpath(
            "out_diploid_phased_phased.vcf.gz"
        )
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=sample_input,
            interval=self.bed,
        )
        driver.add_algo(
            VariantPhaser(
                diploid_vcf,
                phased_vcf,
                out_bed=phased_bed,
                out_ext=phased_ext,
            )
        )
        phaser_job = Job(
            shlex.join(driver.build_cmd()), "variantphaser", self.cores
        )

        bcftools_subset_phased_job = None
        if self.tech.upper() == "ONT":
            bcftools_subset_phased_job = Job(
                f"bcftools view -T {phased_bed} {phased_vcf} \
                | sentieon util vcfconvert - {phased_phased}",
                "bcftools-subset-phased",
                0,
            )

        fai_to_bed_job = None
        if self.bed:
            bed = self.bed
        else:
            bed = self.tmp_dir.joinpath("reference.bed")
            fai_to_bed_job = Job(
                cmds.cmd_fai_to_bed(
                    pathlib.Path(str(self.reference) + ".fai"),
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
        if not self.repeat_model:
            repeat_model = self.tmp_dir.joinpath("out_repeat.model")
            driver = Driver(
                reference=self.reference,
                thread_count=self.cores,
                input=sample_input,
                interval=phased_bed,
                read_filter=[
                    f"PhasedReadFilter,phased_vcf={phased_ext},"
                    "phase_select=tag"
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
                shlex.join(driver.build_cmd()), "repeatmodel", self.cores
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
            hp_std_vcf = self.tmp_dir.joinpath(
                f"out_hap{phase}_nohp_tmp.vcf.gz"
            )
            hp_vcf = self.tmp_dir.joinpath(f"out_hap{phase}_tmp.vcf.gz")
            driver = Driver(
                reference=self.reference,
                thread_count=self.cores,
                input=sample_input,
                interval=phased_bed,
                read_filter=[
                    f"PhasedReadFilter,phased_vcf={phased_ext}"
                    f",phase_select={phase}"
                ],
            )

            if self.tech.upper() == "HIFI":
                # ONT doesn't do DNAscope in 2nd pass.
                driver.add_algo(
                    DNAscope(
                        hp_std_vcf,
                        dbsnp=self.dbsnp,
                        model=self.model_bundle.joinpath("haploid_model"),
                    )
                )
            driver.add_algo(
                DNAscopeHP(
                    hp_vcf,
                    dbsnp=self.dbsnp,
                    model=self.model_bundle.joinpath("haploid_hp_model"),
                    pcr_indel_model=self.repeat_model,
                )
            )
            second_calling_job.add(
                Job(shlex.join(driver.build_cmd()), "second-pass", self.cores)
            )

        kwargs: Dict[str, str] = dict()
        kwargs["gvcf_combine_py"] = str(
            files("sentieon_cli.scripts").joinpath("gvcf_combine.py")
        )
        kwargs["vcf_mod_py"] = str(
            files("sentieon_cli.scripts").joinpath("vcf_mod.py")
        )

        patch_vcfs = [
            self.tmp_dir.joinpath(f"out_hap{i}_patch.vcf.gz") for i in (1, 2)
        ]
        haploid_patch_job = Job(
            cmds.cmd_pyexec_vcf_mod_haploid_patch(
                str(patch_vcfs[0]),
                str(patch_vcfs[1]),
                f"{self.tmp_dir}/out_hap%d_%stmp.vcf.gz",
                self.tech,
                str(phased_phased),
                self.cores,
                kwargs,
            ),
            "patch",
            self.cores,
        )

        # apply trained model to the patched vcfs.
        hap_vcfs = [
            self.tmp_dir.joinpath(f"out_hap{i}.vcf.gz") for i in (1, 2)
        ]
        second_modelapply_job = set()
        for patch_vcf, hap_vcf in zip(patch_vcfs, hap_vcfs):
            driver = Driver(
                reference=self.reference,
                thread_count=self.cores,
            )
            driver.add_algo(
                DNAModelApply(
                    self.model_bundle.joinpath("haploid_model"),
                    patch_vcf,
                    hap_vcf,
                )
            )
            second_modelapply_job.add(
                Job(
                    shlex.join(driver.build_cmd()),
                    "second-modelapply",
                    self.cores,
                )
            )

        # Second pass - unphased regions
        diploid_unphased_hp = self.tmp_dir.joinpath(
            "out_diploid_phased_unphased_hp.vcf.gz"
        )
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=sample_input,
            interval=unphased_bed,
        )
        driver.add_algo(
            DNAscopeHP(
                diploid_unphased_hp,
                dbsnp=self.dbsnp,
                model=self.model_bundle.joinpath("diploid_hp_model"),
                pcr_indel_model=self.repeat_model,
            )
        )
        calling_unphased_job = Job(
            shlex.join(driver.build_cmd()), "calling-unphased", self.cores
        )

        # Patch DNA and DNAHP variants
        diploid_unphased_patch = self.tmp_dir.joinpath(
            "out_diploid_unphased_patch.vcf.gz"
        )
        diploid_unphased = self.tmp_dir.joinpath("out_diploid_unphased.vcf.gz")
        cmd = cmds.cmd_pyexec_vcf_mod_patch(
            str(diploid_unphased_patch),
            str(phased_unphased),
            str(diploid_unphased_hp),
            self.cores,
            kwargs,
        )
        diploid_patch_job = Job(cmd, "diploid-patch", self.cores)
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            DNAModelApply(
                self.model_bundle.joinpath("diploid_model_unphased"),
                diploid_unphased_patch,
                diploid_unphased,
            )
        )
        modelapply_unphased_job = Job(
            shlex.join(driver.build_cmd()), "modelapply-unphased", self.cores
        )

        # merge calls to create the output
        diploid_merged_vcf = self.tmp_dir.joinpath("out_diploid_merged.vcf.gz")
        merge_out_vcf = (
            diploid_merged_vcf if self.haploid_bed else self.output_vcf
        )
        merge_job = Job(
            cmds.cmd_pyexec_vcf_mod_merge(
                str(hap_vcfs[0]),
                str(hap_vcfs[1]),
                str(diploid_unphased),
                str(phased_vcf),
                str(phased_bed),
                str(merge_out_vcf),
                self.cores,
                kwargs,
            ),
            "merge",
            self.cores,
        )

        gvcf_combine_job = None
        if self.gvcf:
            gvcf_combine_job = Job(
                cmds.cmd_pyexec_gvcf_combine(
                    self.reference,
                    str(diploid_gvcf_fn),
                    str(merge_out_vcf),
                    self.cores,
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
        if self.haploid_bed:
            # Haploid variant calling
            haploid_fn = self.tmp_dir.joinpath("haploid.vcf.gz")
            haploid_gvcf_fn = self.tmp_dir.joinpath("haploid.g.vcf.gz")
            haploid_hp_fn = self.tmp_dir.joinpath("haploid_hp.vcf.gz")
            haploid_out_fn = self.tmp_dir.joinpath("haploid_patched.vcf.gz")
            driver = Driver(
                reference=self.reference,
                thread_count=self.cores,
                input=sample_input,
                interval=self.haploid_bed,
            )
            driver.add_algo(
                DNAscope(
                    haploid_fn,
                    dbsnp=self.dbsnp,
                    model=self.model_bundle.joinpath("haploid_model"),
                )
            )
            driver.add_algo(
                DNAscopeHP(
                    haploid_hp_fn,
                    dbsnp=self.dbsnp,
                    model=self.model_bundle.joinpath("haploid_hp_model"),
                    pcr_indel_model=self.repeat_model,
                )
            )
            if self.gvcf:
                driver.add_algo(
                    DNAscope(
                        haploid_gvcf_fn,
                        dbsnp=self.dbsnp,
                        emit_mode="gvcf",
                        ploidy=1,
                        model=self.model_bundle.joinpath("gvcf_model"),
                    )
                )
            haploid_calling_job = Job(
                shlex.join(driver.build_cmd()), "haploid-calling", self.cores
            )

            haploid_patch2_job = Job(
                cmds.cmd_pyexec_vcf_mod_haploid_patch2(
                    str(haploid_out_fn),
                    str(haploid_fn),
                    str(haploid_hp_fn),
                    self.cores,
                    kwargs,
                ),
                "haploid-patch2",
                self.cores,
            )
            haploid_concat_job = Job(
                cmds.bcftools_concat(
                    self.output_vcf,
                    [diploid_merged_vcf, haploid_out_fn],
                ),
                "haploid-diploid-concat",
                0,
            )

            if self.gvcf:
                output_gvcf = pathlib.Path(
                    str(self.output_vcf).replace(".vcf.gz", ".g.vcf.gz")
                )
                diploid_gvcf = pathlib.Path(
                    str(diploid_merged_vcf).replace(".vcf.gz", ".g.vcf.gz")
                )
                haploid_gvcf = pathlib.Path(
                    str(haploid_out_fn).replace(".vcf.gz", ".g.vcf.gz")
                )
                haploid_gvcf_combine_job = Job(
                    cmds.cmd_pyexec_gvcf_combine(
                        self.reference,
                        str(haploid_gvcf_fn),
                        str(haploid_out_fn),
                        self.cores,
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
        self,
        sample_input: List[pathlib.Path],
    ) -> Job:
        """
        Call SVs using Sentieon LongReadSV
        """
        assert self.model_bundle
        if not self.skip_version_check:
            for cmd, min_version in SV_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        sv_vcf = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", ".sv.vcf.gz")
        )
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=sample_input,
            interval=self.bed,
        )
        driver.add_algo(
            LongReadSV(
                sv_vcf,
                self.model_bundle.joinpath("longreadsv.model"),
            )
        )
        longreadsv_job = Job(
            shlex.join(driver.build_cmd()), "LongReadSV", self.cores
        )
        return longreadsv_job
