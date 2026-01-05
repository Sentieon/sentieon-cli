"""
Sentieon's pangenome alignment and variant calling pipeline
"""

import argparse
import copy
import json
import pathlib
import re
import shutil
import subprocess as sp
import sys
from typing import Dict, List, NamedTuple, Optional, Set, Tuple, Union

import packaging.version

from importlib_resources import files

from . import command_strings as cmds
from .archive import ar_load
from .base_pangenome import BasePangenome
from .dag import DAG
from .driver import (
    AlignmentStat,
    BaseDistributionByCycle,
    CoverageMetrics,
    Dedup,
    DNAModelApply,
    DNAscope,
    Driver,
    GCBias,
    InsertSizeMetricAlgo,
    LocusCollector,
    MeanQualityByCycle,
    QualDistribution,
    WgsMetricsAlgo,
)
from .job import Job
from .logging import get_logger
from .shell_pipeline import Command, Pipeline
from .util import (
    __version__,
    check_version,
    parse_rg_line,
    path_arg,
    tmp,
    total_memory,
)

SENT_PANGENOME_MIN_VERSIONS = {
    "kmc": None,
    "sentieon driver": packaging.version.Version("202503.02"),
    "vg": None,
    "bcftools": packaging.version.Version("1.22"),
    "samtools": packaging.version.Version("1.16"),
}

GRCH38_CONTIGS: Dict[str, int] = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}

logger = get_logger(__name__)


class Shard(NamedTuple):
    contig: str
    start: int
    stop: int

    def __str__(self) -> str:
        return f"{self.contig}:{self.start}-{self.stop}"

    def bcftools_str(self) -> str:
        return f"{{{self.contig}}}:{self.start}-{self.stop}"


def parse_fai(ref_fai: pathlib.Path) -> Dict[str, Dict[str, int]]:
    """Parse a faidx index"""
    contigs: Dict[str, Dict[str, int]] = {}
    with open(ref_fai) as fh:
        for line in fh:
            try:
                chrom, length, offset, lb, lw = line.rstrip().split()
            except ValueError as err:
                logger.error(
                    "Reference fasta index (.fai) does not have the expected "
                    "format"
                )
                raise err
            contigs[chrom] = {
                "length": int(length),
                "offset": int(offset),
                "linebases": int(lb),
                "linewidth": int(lw),
            }
    return contigs


def determine_shards_from_fai(
    fai_data: Dict[str, Dict[str, int]], step: int
) -> List[Shard]:
    """Generate shards of the genome from the fasta index"""
    shards: List[Shard] = []
    for ctg, d in fai_data.items():
        pos = 1
        length = d["length"]
        while pos <= length:
            end = pos + step - 1
            end = end if end < length else length
            shards.append(Shard(ctg, pos, end))
            pos = end + 1
    return shards


def vcf_contigs(
    in_vcf: pathlib.Path, dry_run=False
) -> Dict[str, Optional[int]]:
    """Report the contigs in the input VCF"""
    if dry_run:
        return {
            "chr1": 100,
            "chr2": 200,
            "chr3": 300,
        }
    kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
    cmd = ["bcftools", "view", "-h", str(in_vcf)]
    p = sp.run(cmd, capture_output=True, text=True)
    contigs: Dict[str, Optional[int]] = {}
    for line in p.stdout.split("\n"):
        if not line.startswith("##contig"):
            continue
        s = line.index("<")
        e = line.index(">")
        d = dict(kvpat.findall(line[s + 1 : e]))  # noqa: E203
        ctg: str = d["ID"]
        length: Optional[str] = d.get("length", None)
        contigs[ctg] = int(length) if length else None
    return contigs


def vcf_id(in_vcf: pathlib.Path) -> Optional[str]:
    """Collect the SentieonVcfID header"""
    cmd = ["bcftools", "view", "-h", str(in_vcf)]
    p = sp.run(cmd, capture_output=True, text=True)
    for line in p.stdout.split("\n"):
        if line.startswith("##SentieonVcfID="):
            i = line.index("=")
            return line[i + 1 :]  # noqa: E203
    return None


class SentieonPangenome(BasePangenome):
    """The Sentieon pangenome pipeline"""

    params = copy.deepcopy(BasePangenome.params)
    params.update(
        {
            # Required arguments
            "readgroup": {
                "help": "Readgroup information for the fastq files.",
            },
            "sample_input": {
                "flags": ["-i", "--sample_input"],
                "nargs": "*",
                "help": "sample BAM or CRAM file.",
                "type": path_arg(exists=True, is_file=True),
            },
            "pop_vcf": {
                "flags": ["--pop_vcf"],
                "help": (
                    "A VCF containing annotations for use with DNAModelApply."
                ),
                "type": path_arg(exists=True, is_file=True),
                "required": True,
            },
            # Additional arguments
            "bed": {
                "flags": ["-b", "--bed"],
                "help": (
                    "Region BED file. Supplying this file will limit variant "
                    "calling to the intervals inside the BED file."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "skip_metrics": {
                "help": "Skip metrics collection and multiQC",
                "action": "store_true",
            },
            "skip_multiqc": {
                "help": "Skip multiQC report generation",
                "action": "store_true",
            },
            # Hidden arguments
            "skip_contig_checks": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_pangenome_name_checks": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_pop_vcf_id_check": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
        }
    )

    positionals = BasePangenome.positionals

    def __init__(self) -> None:
        super().__init__()
        self.readgroup: Optional[str] = None
        self.sample_input: List[pathlib.Path] = []
        self.bed: Optional[pathlib.Path] = None
        self.pop_vcf: Optional[pathlib.Path] = None
        self.skip_metrics = False
        self.skip_multiqc = False
        self.skip_contig_checks: bool = False
        self.skip_pangenome_name_checks: bool = False
        self.skip_pop_vcf_id_check: bool = False

    def main(self, args: argparse.Namespace) -> None:
        """Run the pipeline"""
        self.handle_arguments(args)
        self.setup_logging(args)
        self.validate_ref()

        self.fai_data = parse_fai(pathlib.Path(str(self.reference) + ".fai"))
        self.pop_vcf_contigs: Dict[str, Optional[int]] = {}
        if self.pop_vcf:
            self.pop_vcf_contigs = vcf_contigs(self.pop_vcf, self.dry_run)
            self.logger.debug("VCF contigs are: %s", self.pop_vcf_contigs)

        self.validate()
        self.shards = determine_shards_from_fai(
            self.fai_data, 10 * 1000 * 1000
        )

        tmp_dir_str = tmp()
        self.tmp_dir = pathlib.Path(tmp_dir_str)

        dag = self.build_first_dag()
        executor = self.run(dag)
        self.check_execution(dag, executor)

        if not self.retain_tmpdir:
            shutil.rmtree(tmp_dir_str)

    def validate(self) -> None:
        """Validate pipeline inputs"""
        self.validate_bundle()
        self.validate_fastq_rg()
        self.validate_output_vcf()
        self.validate_ref()
        self.collect_readgroups()

        if not self.sample_input and not self.r1_fastq:
            self.logger.error(
                "Please supply either the `--sample_input` or `--r1_fastq` "
                "and `--readgroups` arguments"
            )
            sys.exit(2)

        if self.sample_input and self.r1_fastq:
            self.logger.error(
                "Supplying both `--r1_fastq` and `--sample_input` is not "
                "supported"
            )
            sys.exit(2)

        if not self.skip_version_check:
            for cmd, min_version in SENT_PANGENOME_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        if self.bed is None:
            self.logger.info(
                "A BED file is recommended to avoid variant calling "
                "across decoy and unplaced contigs."
            )

        if not self.skip_pangenome_name_checks:
            if not str(self.gbz).endswith("grch38.gbz"):
                self.logger.error(
                    "The `--gbz` file does not have the expected suffix. "
                    "Check that you are using a GRCh38 pangenome."
                )
                sys.exit(2)

            if self.gbz.name != "hprc-v2.0-mc-grch38.gbz":
                self.logger.warning(
                    "The `--gbz` file name is not 'hprc-v2.0-mc-grch38.gbz'. "
                    "This pipeline is optimized for the HPRC v2.0 pangenome."
                )

            if not str(self.hapl).endswith("grch38.hapl"):
                self.logger.error(
                    "The `--hapl` file does not have the expected suffix. "
                    "Check that you re using a GRCh38 pangenome."
                )
                sys.exit(2)

        if not self.skip_contig_checks:
            # Check the fai file contigs
            mismatch_contigs: Set[str] = set()
            for ctg, length in GRCH38_CONTIGS.items():
                d = self.fai_data.get(ctg, {})
                fai_length = d.get("length", -1)
                if length != fai_length:
                    mismatch_contigs.add(ctg)
            if mismatch_contigs:
                mismatch_contigs_s = ", ".join(mismatch_contigs)
                self.logger.error(
                    "Reference contigs with unexpected lengths: %s",
                    mismatch_contigs_s,
                )
                sys.exit(2)

            # Check the pop VCF file contigs
            if not self.dry_run:
                mismatch_contigs = set()
                for ctg, length in GRCH38_CONTIGS.items():
                    vcf_length = self.pop_vcf_contigs.get(ctg, -1)
                    if length != vcf_length:
                        mismatch_contigs.add(ctg)
                if mismatch_contigs:
                    mismatch_contigs_s = ", ".join(mismatch_contigs)
                    self.logger.error(
                        "Pop VCF contigs with unexpected lengths: %s",
                        mismatch_contigs_s,
                    )
                    sys.exit(2)

    def validate_bundle(self) -> None:
        assert self.pop_vcf
        bundle_info_bytes = ar_load(
            str(self.model_bundle) + "/bundle_info.json"
        )
        if isinstance(bundle_info_bytes, list):
            bundle_info_bytes = b"{}"

        bundle_info = json.loads(bundle_info_bytes.decode())
        try:
            req_version = packaging.version.Version(
                bundle_info["minScriptVersion"]
            )
            bundle_pipeline = bundle_info["pipeline"]
            self.tech: str = bundle_info["platform"].upper()
            bundle_vcf_id = bundle_info["SentieonVcfID"]
        except KeyError:
            self.logger.error(
                "The model bundle does not have the expected attributes"
            )
            sys.exit(2)
        if req_version > packaging.version.Version(__version__):
            self.logger.error(
                "The model bundle requires version %s or later of the "
                "sentieon-cli.",
                req_version,
            )
            sys.exit(2)
        if bundle_pipeline != "Sentieon pangenome":
            self.logger.error("The model bundle is for a different pipeline.")
            sys.exit(2)

        bundle_members = set(ar_load(str(self.model_bundle)))
        if (
            "dnascope.model" not in bundle_members
            or "extract.model" not in bundle_members
            or "minimap2.model" not in bundle_members
        ):
            self.logger.error(
                "Expected model files not found in the model bundle file"
            )
            sys.exit(2)

        if not self.skip_pop_vcf_id_check and not self.dry_run:
            pop_vcf_id = vcf_id(self.pop_vcf)
            if bundle_vcf_id != pop_vcf_id:
                self.logger.error(
                    "The ID of the `--pop_vcf` does not match the model bundle"
                )
                sys.exit(2)

    def validate_fastq_rg(self) -> None:
        if len(self.r1_fastq) != len(self.r2_fastq):
            self.logger.error(
                "The number of input `--r1_fastq` files does not equal the "
                "number of `--r2_fastq` files"
            )
            sys.exit(2)

        if (len(self.r1_fastq) > 0 and not self.readgroup) or (
            self.readgroup and len(self.r1_fastq) < 1
        ):
            self.logger.error(
                "`--r1_fastq`, `--r2_fastq`, and `--readgroup` are required "
                "with fastq input. `--readgroup` cannot be used with bam/cram "
                "input."
            )
            sys.exit(2)

    def collect_readgroups(self) -> None:
        """Collect readgroup tags"""
        self.bam_readgroups: List[Dict[str, str]] = []
        for aln in self.sample_input:
            aln_rgs = cmds.get_rg_lines(aln, self.dry_run)
            for rg_line in aln_rgs:
                self.bam_readgroups.append(parse_rg_line(rg_line))
                break  # currently just the first RG line for each input

        self.fastq_readgroup: Dict[str, str] = {}
        if not self.readgroup:
            return

        parsed_rg = parse_rg_line(self.readgroup.replace(r"\t", "\t"))
        if not parsed_rg.get("ID"):
            self.logger.error(
                "Readgroup '%s' does not have a RGID tag",
                self.readgroup,
            )
            sys.exit(2)
        if parsed_rg.get("SM", None) is None:
            self.logger.error(
                "Readgroup '%s' does not have a RGSM tag",
                self.readgroup,
            )
        self.fastq_readgroup = parsed_rg

    def configure(self) -> None:
        """Configure pipeline parameters"""
        pass

    def build_dag(self) -> DAG:
        return DAG()

    def build_first_dag(self) -> DAG:
        """Build the main DAG for the Sentieon pangenome pipeline"""
        assert self.reference
        assert self.model_bundle

        self.logger.info("Building the Sentieon pangenome DAG")
        dag = DAG()

        # Output files
        suffix = "bam" if self.bam_format else "cram"
        out_bwa_aln = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", f"_bwa_deduped.{suffix}")
        )
        out_mm2_aln = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", f"_mm2_deduped.{suffix}")
        )
        self.ploidy_json = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_ploidy.json")
        )

        # Intermediate file paths
        bwa_bam = self.tmp_dir.joinpath("sample-bwa.bam")
        ext_fastq = self.tmp_dir.joinpath("sample-bwa.fq.gz")
        kmer_prefix = self.tmp_dir.joinpath("sample.fq")
        kmer_file = pathlib.Path(str(kmer_prefix) + ".kff")
        sample_pangenome = self.tmp_dir.joinpath("sample_pangenome.gbz")
        sample_gfa = self.tmp_dir.joinpath("sample-hap.gfa")
        sample_fasta = self.tmp_dir.joinpath("sample-hap.fa")
        mm2_bam = self.tmp_dir.joinpath("sample-mm2.bam")
        raw_vcf = self.tmp_dir.joinpath("sample-dnascope.vcf.gz")
        transfer_vcf = self.tmp_dir.joinpath("sample-dnascope_transfer.vcf.gz")
        # - Ensure we have a bam file for t1k
        self.rw_bam = self.tmp_dir.joinpath("sample_deduped.bam")

        bwa_lc_dependencies: Set[Job] = set()
        haplotype_dependencies: Set[Job] = set()
        mm2_dependencies: Set[Job] = set()
        dnascope_bams: List[pathlib.Path] = []

        total_mem_gb = total_memory() / (1024.0**3)

        if self.r1_fastq:
            # KMC k-mer counting
            kmc_job = self.build_kmc_job(kmer_prefix, 0)  # run in background
            dag.add_job(kmc_job)
            haplotype_dependencies.add(kmc_job)

            # BWA alignment and extraction
            bwa_job = self.build_alignment_job(bwa_bam, ext_fastq)
            dag.add_job(bwa_job)
            mm2_dependencies.add(bwa_job)
            bwa_lc_dependencies.add(bwa_job)
            # Do not run vg-haplotypes with bwa in low-mem environments
            if total_mem_gb < 70:
                haplotype_dependencies.add(bwa_job)
        else:
            dnascope_bams = copy.deepcopy(self.sample_input)
            extract_kmc_job = Job(
                cmds.cmd_extract_kmc(
                    kmer_prefix,
                    ext_fastq,
                    self.sample_input,
                    self.reference,
                    self.model_bundle.joinpath("extract.model"),
                    self.tmp_dir,
                    threads=self.cores,
                ),
                "extract-kmc",
                self.cores,
            )
            dag.add_job(extract_kmc_job)
            haplotype_dependencies.add(extract_kmc_job)

        # vg haplotypes - create a sample-specific pangenome
        haplotypes_job = self.build_haplotypes_job(sample_pangenome, kmer_file)
        dag.add_job(haplotypes_job, haplotype_dependencies)

        # convert the sample pangenome
        gfa_job = self.build_gfa_job(sample_gfa, sample_pangenome)
        fasta_job = self.build_fasta_job(sample_fasta, sample_pangenome)
        dag.add_job(gfa_job, {haplotypes_job})
        dag.add_job(fasta_job, {haplotypes_job})

        # minimap2 alignment of the extracted fastq
        dnascope_dependencies = set()
        mm2_job = self.build_minimap2_lift_job(
            mm2_bam,
            ext_fastq,
            sample_fasta,
            sample_gfa,
        )
        dag.add_job(mm2_job, mm2_dependencies | {gfa_job, fasta_job})
        dnascope_dependencies.add(mm2_job)

        # With fastq input, perform dedup and metrics
        if self.r1_fastq:
            dnascope_bams.append(out_bwa_aln)
            dnascope_bams.append(out_mm2_aln)
            bwa_lc_job, bwa_dedup_job = self.build_dedup_job(
                out_bwa_aln, [bwa_bam], "bwa"
            )
            mm2_lc_job, mm2_dedup_job = self.build_dedup_job(
                out_mm2_aln, [mm2_bam], "mm2", left_align=True
            )
            dag.add_job(bwa_lc_job, bwa_lc_dependencies)
            dag.add_job(bwa_dedup_job, {bwa_lc_job})
            dnascope_dependencies.add(bwa_dedup_job)
            dag.add_job(mm2_lc_job, {mm2_job})
            dag.add_job(mm2_dedup_job, {mm2_lc_job})
            dnascope_dependencies.add(mm2_dedup_job)

            if not self.skip_metrics:
                metrics_job, rehead_job = self.build_metrics_job(
                    [out_bwa_aln, out_mm2_aln],
                )
                dag.add_job(metrics_job, {bwa_dedup_job, mm2_dedup_job})
                dag.add_job(rehead_job, {metrics_job})
                if not self.skip_multiqc:
                    multiqc_job = self.multiqc()
                    if multiqc_job:
                        dag.add_job(multiqc_job, {rehead_job})
        else:
            dnascope_bams.append(mm2_bam)

        # DNAscope calling with bwa and mm2 input
        apply_dependencies = set()
        dnascope_job = self.build_dnascope_job(raw_vcf, dnascope_bams)
        dag.add_job(dnascope_job, dnascope_dependencies)
        apply_dependencies.add(dnascope_job)

        # transfer annotations from the pop_vcf
        if self.pop_vcf:
            transfer_jobs, concat_job = self.build_transfer_jobs(
                transfer_vcf, raw_vcf
            )
            for job in transfer_jobs:
                dag.add_job(job, {dnascope_job})
            dag.add_job(concat_job, set(transfer_jobs))
            apply_dependencies.add(concat_job)

        # DNAModelApply
        apply_job = self.build_dnamodelapply_job(transfer_vcf)
        dag.add_job(apply_job, apply_dependencies)

        return dag

    def build_alignment_job(
        self,
        sample_bam: pathlib.Path,
        sample_fastq: pathlib.Path,
    ) -> Job:
        """Build the alignment and extract jobs"""
        assert self.reference
        assert self.model_bundle

        unzip = "igzip"
        if not shutil.which(unzip):
            self.logger.info(
                "igzip is recommended for decompression, but is not "
                "available. Falling back to gzip."
            )
            unzip = "gzip"

        rg = copy.deepcopy(self.fastq_readgroup)
        rg["ID"] = rg["ID"] + "-bwa"
        bwa_job = Job(
            cmds.cmd_bwa_extract(
                sample_bam,
                sample_fastq,
                self.reference,
                self.r1_fastq,
                self.r2_fastq,
                "@RG\\t" + "\\t".join([f"{x[0]}:{x[1]}" for x in rg.items()]),
                self.model_bundle.joinpath("extract.model"),
                self.model_bundle.joinpath("bwa.model"),
                self.cores,
                unzip=unzip,
            ),
            "bwa-extract",
            self.cores,
        )
        return bwa_job

    def build_haplotypes_job(
        self, output_gbz: pathlib.Path, kmer_file: pathlib.Path
    ) -> Job:
        """Build vg haplotypes job"""
        assert self.hapl
        assert self.gbz

        haplotypes_job = Job(
            cmds.cmd_vg_haplotypes(
                output_gbz,
                kmer_file,
                self.hapl,
                self.gbz,
                threads=self.cores,
                xargs=[
                    "--include-reference",
                    "--diploid-sampling",
                    "--set-reference",
                    "GRCh38",
                ],
            ),
            "vg-haplotypes",
            self.cores,
        )
        return haplotypes_job

    def build_gfa_job(
        self, output_gfa: pathlib.Path, input_gbz: pathlib.Path
    ) -> Job:
        """Build vg convert to GFA job"""
        gfa_job = Job(
            cmds.cmd_vg_convert_gfa(
                output_gfa,
                input_gbz,
                threads=self.cores,
            ),
            "vg-convert-gfa",
            0,
        )
        return gfa_job

    def build_fasta_job(
        self, output_fasta: pathlib.Path, input_gbz: pathlib.Path
    ) -> Job:
        """Build vg paths to FASTA job"""
        fasta_job = Job(
            cmds.cmd_vg_paths_fasta(
                output_fasta,
                input_gbz,
            ),
            "vg-paths-fasta",
            0,
        )
        return fasta_job

    def build_minimap2_lift_job(
        self,
        mm2_bam: pathlib.Path,
        ext_fastq: pathlib.Path,
        sample_fasta: pathlib.Path,
        sample_gfa: pathlib.Path,
    ) -> Job:
        """Build minimap2 alignment with pgutil lift job"""
        assert self.model_bundle
        assert self.reference

        rg = (
            self.fastq_readgroup
            if self.fastq_readgroup
            else self.bam_readgroups[0]
        )
        rg2 = copy.deepcopy(rg)
        rg2["ID"] = rg2["ID"] + "-mm2"
        rg2["LR"] = "1"

        mm2_model: Union[str, pathlib.Path] = self.model_bundle.joinpath(
            "minimap2.model"
        )
        mm2_job = Job(
            cmds.cmd_minimap2_lift(
                mm2_bam,
                sample_fasta,
                ext_fastq,
                sample_gfa,
                self.reference,
                "@RG\\t" + "\\t".join([f"{x[0]}:{x[1]}" for x in rg2.items()]),
                mm2_model,
                threads=self.cores,
            ),
            "mm2-lift",
            self.cores,
        )
        return mm2_job

    def build_dedup_job(
        self,
        output_bam,
        input_bam: List[pathlib.Path],
        tag: str,
        left_align=False,
    ) -> Tuple[Job, Job]:
        """Build deduplication job"""
        score_file = self.tmp_dir.joinpath(f"sample-{tag}-score.txt.gz")

        read_filters = []
        if left_align:
            read_filters.append(
                "IndelLeftAlignReadTransform,"
                f"rgid={self.fastq_readgroup['ID']}-mm2"
            )

        # LocusCollector + Dedup
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bam,
            read_filter=read_filters,
        )
        driver.add_algo(LocusCollector(score_file))

        lc_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            f"locuscollector-{tag}",
            self.cores,
        )

        driver2 = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bam,
            read_filter=read_filters,
        )
        driver2.add_algo(Dedup(output_bam, score_file))

        dedup_job = Job(
            Pipeline(Command(*driver2.build_cmd())),
            f"dedup-{tag}",
            self.cores,
        )

        return lc_job, dedup_job

    def build_metrics_job(
        self,
        sample_input: List[pathlib.Path],
    ) -> Tuple[Job, Job]:
        """Build a metrics job"""
        assert self.output_vcf

        # Create the metrics directory
        sample_name = self.output_vcf.name.replace(".vcf.gz", "")
        metric_base = sample_name + ".txt"
        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )
        if not self.dry_run:
            metrics_dir.mkdir(exist_ok=True)

        is_metrics = metrics_dir.joinpath(metric_base + ".insert_size.txt")
        mqbc_metrics = metrics_dir.joinpath(
            metric_base + ".mean_qual_by_cycle.txt"
        )
        bdbc_metrics = metrics_dir.joinpath(
            metric_base + ".base_distribution_by_cycle.txt"
        )
        qualdist_metrics = metrics_dir.joinpath(
            metric_base + ".qual_distribution.txt"
        )
        as_metrics = metrics_dir.joinpath(metric_base + ".alignment_stat.txt")
        coverage_metrics = metrics_dir.joinpath("coverage")

        # WGS metrics
        wgs_metrics = metrics_dir.joinpath(metric_base + ".wgs.txt")
        gc_metrics = metrics_dir.joinpath(metric_base + ".gc_bias.txt")
        gc_summary = metrics_dir.joinpath(metric_base + ".gc_bias_summary.txt")

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=sample_input,
        )

        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        driver.add_algo(MeanQualityByCycle(mqbc_metrics))
        driver.add_algo(BaseDistributionByCycle(bdbc_metrics))
        driver.add_algo(QualDistribution(qualdist_metrics))
        driver.add_algo(AlignmentStat(as_metrics))
        driver.add_algo(GCBias(gc_metrics, summary=gc_summary))
        driver.add_algo(WgsMetricsAlgo(wgs_metrics, include_unpaired="true"))
        driver.add_algo(CoverageMetrics(coverage_metrics))

        metrics_job = Job(Pipeline(Command(*driver.build_cmd())), "metrics", 0)

        rehead_script = pathlib.Path(
            str(
                files("sentieon_cli.scripts").joinpath("rehead_wgs_metrics.py")
            )
        )
        rehead_job = Job(
            Pipeline(
                Command(
                    "sentieon",
                    "pyexec",
                    str(rehead_script),
                    "--metrics_file",
                    str(wgs_metrics),
                )
            ),
            "Rehead metrics",
            0,
        )
        return (metrics_job, rehead_job)

    def build_dnascope_job(
        self,
        out_vcf: pathlib.Path,
        input_bams: List[pathlib.Path],
    ) -> Job:
        assert self.model_bundle

        read_filters = []
        if self.tech.upper() == "ULTIMA":
            read_filters.append("UltimaReadFilter")

        pcr_indel_model = "NONE" if self.pcr_free else "CONSERVATIVE"
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bams,
            interval=self.bed,
            read_filter=read_filters,
        )
        driver.add_algo(
            DNAscope(
                out_vcf,
                model=self.model_bundle.joinpath("dnascope.model"),
                pcr_indel_model=pcr_indel_model,
                dbsnp=self.dbsnp,
            )
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            "dnascope-raw",
            self.cores,
        )

    def build_transfer_jobs(
        self, out_vcf: pathlib.Path, raw_vcf: pathlib.Path
    ) -> Tuple[List[Job], Job]:
        """Transfer annotations from the pop_vcf to the raw_vcf"""
        assert self.pop_vcf

        # Generate merge rules from the population VCF
        merge_rules = "AC_v20:sum,AF_v20:sum,AC_genomes:sum,AF_genomes:sum"
        if not self.dry_run:
            kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
            cmd = ["bcftools", "view", "-h", str(self.pop_vcf)]
            p = sp.run(cmd, capture_output=True, text=True)
            id_fields: List[str] = []
            for line in p.stdout.split("\n"):
                if not line.startswith("##INFO"):
                    continue
                if ",Number=A" not in line:
                    continue
                s = line.index("<")
                e = line.index(">")
                d = dict(kvpat.findall(line[s + 1 : e]))  # noqa: E203
                id_fields.append(d["ID"])
            merge_rules = ",".join([x + ":sum" for x in id_fields])

        # Merge VCFs by shards
        sharded_vcfs: List[pathlib.Path] = []
        sharded_merge_jobs: List[Job] = []
        trim_script = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("trimalt.py"))
        ).resolve()
        seen_contigs: Set[str] = set()
        for i, shard in enumerate(self.shards):
            # Use a BED file for unusual contig names
            subset_bed = self.tmp_dir.joinpath(
                f"sample-dnascope_transfer-subset{i}.bed"
            )

            # Extract contigs not in the pop vcf as merge will fail
            if shard.contig not in self.pop_vcf_contigs:
                if shard.contig in seen_contigs:
                    continue
                self.logger.info(
                    "Skipping transfer for contig: %s", shard.contig
                )
                seen_contigs.add(shard.contig)
                subset_vcf = self.tmp_dir.joinpath(
                    f"sample-dnascope_transfer-subset{i}.vcf.gz"
                )

                ctg_len = self.fai_data[shard.contig]["length"]
                if not self.dry_run:
                    with open(subset_bed, "w") as fh:
                        print(f"{shard.contig}\t0\t{ctg_len}", file=fh)

                view_job = Job(
                    cmds.cmd_bcftools_view_regions(
                        subset_vcf,
                        raw_vcf,
                        regions_file=subset_bed,
                    ),
                    "merge-trim-extra",
                    1,
                )
                sharded_merge_jobs.append(view_job)
                sharded_vcfs.append(subset_vcf)
            else:
                if not self.dry_run:
                    with open(subset_bed, "w") as fh:
                        print(
                            f"{shard.contig}\t{shard.start}\t{shard.stop}",
                            file=fh,
                        )

                self.logger.debug("Transferring shard: %s", shard)
                shard_vcf = self.tmp_dir.joinpath(
                    f"sample-dnascope_transfer-shard{i}.vcf.gz"
                )
                merge_job = Job(
                    cmds.cmd_bcftools_merge_trim(
                        shard_vcf,
                        raw_vcf,
                        self.pop_vcf,
                        trim_script,
                        subset_bed,
                        merge_rules=merge_rules,
                        merge_xargs=[
                            "--no-version",
                            "--regions-overlap",
                            "pos",
                            "-m",
                            "all",
                        ],
                        view_xargs=["--no-version"],
                    ),
                    f"merge-trim-{i}",
                    1,
                )
                sharded_merge_jobs.append(merge_job)
                sharded_vcfs.append(shard_vcf)

        # Concat all shards
        concat_job = Job(
            cmds.bcftools_concat(
                out_vcf,
                sharded_vcfs,
                xargs=["--no-version", "--threads", str(self.cores)],
            ),
            "merge-trim-concat",
            self.cores,
        )
        return (sharded_merge_jobs, concat_job)

    def build_dnamodelapply_job(self, in_vcf) -> Job:
        assert self.output_vcf
        assert self.model_bundle

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            DNAModelApply(
                model=self.model_bundle.joinpath("dnascope.model"),
                vcf=in_vcf,
                output=self.output_vcf,
            )
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            "model-apply",
            self.cores,
        )
