"""
Sentieon's pangenome alignment and variant calling pipeline
"""

import argparse
import copy
import itertools
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
    Dedup,
    DNAModelApply,
    DNAscope,
    Driver,
    LocusCollector,
)
from .job import Job
from .logging import get_logger
from .shell_pipeline import Command, Pipeline
from .util import __version__, check_version, parse_rg_line, path_arg, tmp

SENT_PANGENOME_MIN_VERSIONS = {
    "kmc": None,
    "sentieon driver": packaging.version.Version("202503.02"),
    "vg": None,
    "bcftools": packaging.version.Version("1.20"),
}

# With multiple aligned input files, require kff-tools
MULTI_INPUT_MIN_VERSIONS = {
    "kff-tools": None,
}

logger = get_logger(__name__)


class Shard(NamedTuple):
    contig: str
    start: int
    stop: int

    def __str__(self) -> str:
        return f"{self.contig}:{self.start}-{self.stop}"


def determine_shards_from_fai(ref_fai: pathlib.Path, step: int) -> List[Shard]:
    """Generate shards of the genome from the fasta index"""
    shards: List[Shard] = []
    with open(ref_fai) as fh:
        for line in fh:
            pos = 1
            try:
                chrom, length, _offset, _ll, _lb = line.rstrip().split()
            except ValueError:
                logger.error(
                    "Reference fasta index (.fai) does not have the expected "
                    "format"
                )
                raise

            length_i = int(length)
            while pos <= length_i:
                end = pos + step - 1
                end = end if end < length_i else length_i
                shards.append(Shard(chrom, pos, end))
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


class SentieonPangenome(BasePangenome):
    """The Sentieon pangenome pipeline"""

    params = copy.deepcopy(BasePangenome.params)
    params.update(
        {
            # Required arguments
            "sample_input": {
                "flags": ["-i", "--sample_input"],
                "nargs": "*",
                "help": "sample BAM or CRAM file.",
                "type": path_arg(exists=True, is_file=True),
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
            "pop_vcf": {
                "flags": ["--pop_vcf"],
                "help": (
                    "A VCF containing annotations for use with DNAModelApply."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
        }
    )

    positionals = BasePangenome.positionals

    def __init__(self) -> None:
        super().__init__()
        self.sample_input: List[pathlib.Path] = []
        self.bed: Optional[pathlib.Path] = None
        self.pop_vcf: Optional[pathlib.Path] = None

    def main(self, args: argparse.Namespace) -> None:
        """Run the pipeline"""
        self.handle_arguments(args)
        self.setup_logging(args)
        self.validate()
        self.shards = determine_shards_from_fai(
            pathlib.Path(str(self.reference) + ".fai"), 10 * 1000 * 1000
        )
        if self.pop_vcf:
            self.pop_vcf_contigs = vcf_contigs(self.pop_vcf, self.dry_run)

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
            if self.sample_input and len(self.sample_input) > 1:
                for cmd, min_version in MULTI_INPUT_MIN_VERSIONS.items():
                    if not check_version(cmd, min_version):
                        sys.exit(2)

        if self.bed is None:
            self.logger.info(
                "A BED file is recommended to avoid variant calling "
                "across decoy and unplaced contigs."
            )

    def validate_bundle(self) -> None:
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
        ):
            self.logger.error(
                "Expected model files not found in the model bundle file"
            )
            sys.exit(2)

    def collect_readgroups(self) -> None:
        """Collect readgroup tags from all input files"""
        self.fastq_readgroups: List[Dict[str, str]] = []
        sample_sm = None
        for rg in self.readgroups:
            parsed_rg = parse_rg_line(rg.replace(r"\t", "\t"))
            if not parsed_rg.get("ID"):
                self.logger.error(
                    "Readgroup '%s' does not have a RGID tag",
                    rg,
                )
                sys.exit(2)
            if sample_sm:
                if parsed_rg.get("SM") != sample_sm:
                    self.logger.error(
                        "Readgroup '%s' has a RGSM tag inconsistent with "
                        "earlier readgroups",
                        rg,
                    )
                    sys.exit(2)
            self.fastq_readgroups.append(parsed_rg)

        self.bam_readgroups: List[Dict[str, str]] = []
        for aln in self.sample_input:
            aln_rgs = cmds.get_rg_lines(aln, self.dry_run)
            for rg_line in aln_rgs:
                self.bam_readgroups.append(parse_rg_line(rg_line))
                break  # currently just the first RG line for each input

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

        # Intermediate file paths
        kmer_prefix = self.tmp_dir.joinpath("sample.fq")
        kmer_file = pathlib.Path(str(kmer_prefix) + ".kff")
        sample_pangenome = self.tmp_dir.joinpath("sample_pangenome.gbz")
        sample_gfa = self.tmp_dir.joinpath("sample-hap.gfa")
        sample_fasta = self.tmp_dir.joinpath("sample-hap.fa")
        raw_vcf = self.tmp_dir.joinpath("sample-dnascope.vcf.gz")
        transfer_vcf = self.tmp_dir.joinpath("sample-dnascope_transfer.vcf.gz")

        haplotype_dependencies: Set[Job] = set()
        mm2_dependencies: Set[Job] = set()
        bwa_bams: List[pathlib.Path] = []
        lmr_fastqs: List[pathlib.Path] = []
        bwa_jobs: List[Job] = []
        dnascope_bams: List[pathlib.Path] = []
        if self.r1_fastq:
            # KMC k-mer counting
            kmc_job = self.build_kmc_job(kmer_prefix, 0)  # run in background
            dag.add_job(kmc_job)
            haplotype_dependencies.add(kmc_job)

            # BWA alignment and extraction
            bwa_bams, lmr_fastqs, bwa_jobs = self.build_alignment_jobs()
            for bwa_job in bwa_jobs:
                dag.add_job(bwa_job)
                mm2_dependencies.add(bwa_job)
        else:
            dnascope_bams = copy.deepcopy(self.sample_input)
            all_kmer_files: List[pathlib.Path] = []
            all_extract_jobs: Set = set()
            for i, aln in enumerate(self.sample_input):
                kmer_prefix_i = self.tmp_dir.joinpath(f"sample.fq.kmer-{i}")
                kmer_file_i = pathlib.Path(str(kmer_prefix_i) + ".kff")
                lmr_fq = self.tmp_dir.joinpath(f"sample_lmr-{i}.fq.gz")
                extract_kmc_job = Job(
                    cmds.cmd_extract_kmc(
                        kmer_prefix_i,
                        lmr_fq,
                        aln,
                        self.reference,
                        self.model_bundle.joinpath("extract.model"),
                        self.tmp_dir,
                        threads=self.cores,
                    ),
                    "extract-kmc-{i}",
                    self.cores,
                )
                dag.add_job(extract_kmc_job)
                haplotype_dependencies.add(extract_kmc_job)
                mm2_dependencies.add(extract_kmc_job)
                lmr_fastqs.append(lmr_fq)
                all_kmer_files.append(kmer_file_i)
                all_extract_jobs.add(extract_kmc_job)

            if len(all_kmer_files) > 1:
                # Merge kff files
                merge_kff_job = Job(
                    cmds.cmd_merge_kff(
                        kmer_file,
                        all_kmer_files,
                    ),
                    "merge-kff",
                    self.cores,
                )
                dag.add_job(merge_kff_job, all_extract_jobs)
                haplotype_dependencies.add(merge_kff_job)
            else:
                kmer_file = all_kmer_files[0]

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
        mm2_bams, mm2_jobs = self.build_minimap2_lift_job(
            lmr_fastqs,
            sample_fasta,
            sample_gfa,
        )
        for mm2_job in mm2_jobs:
            dag.add_job(mm2_job, mm2_dependencies | {gfa_job, fasta_job})
            dnascope_dependencies.add(mm2_job)

        # With fastq input, perform dedup
        if self.r1_fastq:
            dnascope_bams.append(out_bwa_aln)
            dnascope_bams.append(out_mm2_aln)
            bwa_lc_job, bwa_dedup_job = self.build_dedup_job(
                out_bwa_aln, bwa_bams, "bwa"
            )
            mm2_lc_job, mm2_dedup_job = self.build_dedup_job(
                out_mm2_aln, mm2_bams, "mm2", left_align=True
            )
            dag.add_job(bwa_lc_job, set(bwa_jobs))
            dag.add_job(bwa_dedup_job, {bwa_lc_job})
            dnascope_dependencies.add(bwa_dedup_job)
            dag.add_job(mm2_lc_job, set(mm2_jobs))
            dag.add_job(mm2_dedup_job, {mm2_lc_job})
            dnascope_dependencies.add(mm2_dedup_job)
        else:
            dnascope_bams.extend(mm2_bams)

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

    def build_alignment_jobs(
        self,
    ) -> Tuple[List[pathlib.Path], List[pathlib.Path], List[Job]]:
        """Build the alignment and extract jobs"""
        assert self.reference
        assert self.model_bundle

        # Process each pair of fastq files
        sample_bams: List[pathlib.Path] = []
        sample_fastqs: List[pathlib.Path] = []
        bwa_jobs: List[Job] = []
        for i, (r1, r2, rg) in enumerate(
            itertools.zip_longest(
                self.r1_fastq, self.r2_fastq, self.fastq_readgroups
            )
        ):
            sample_bam = self.tmp_dir.joinpath(f"sample_aligned_{i}.bam")
            sample_fastq = self.tmp_dir.joinpath(f"sample_aligned_{i}.fq.gz")
            rg2 = copy.deepcopy(rg)
            rg2["ID"] = rg2["ID"] + "-bwa"

            bwa_job = Job(
                cmds.cmd_bwa_extract(
                    sample_bam,
                    sample_fastq,
                    self.reference,
                    r1,
                    r2,
                    "@RG\\t"
                    + "\\t".join([f"{x[0]}:{x[1]}" for x in rg2.items()]),
                    self.model_bundle.joinpath("extract.model"),
                    self.model_bundle.joinpath("bwa.model"),
                    self.cores,
                ),
                f"bwa-extract-{i}",
                self.cores,
            )

            sample_bams.append(sample_bam)
            sample_fastqs.append(sample_fastq)
            bwa_jobs.append(bwa_job)
        return (sample_bams, sample_fastqs, bwa_jobs)

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
        lmr_fastqs: List[pathlib.Path],
        sample_fasta: pathlib.Path,
        sample_gfa: pathlib.Path,
    ) -> Tuple[List[pathlib.Path], List[Job]]:
        """Build minimap2 alignment with pgutil lift job"""
        assert self.model_bundle
        assert self.reference

        mm2_bams: List[pathlib.Path] = []
        mm2_jobs: List[Job] = []
        for i, (fq, rg) in enumerate(
            zip(lmr_fastqs, self.fastq_readgroups + self.bam_readgroups)
        ):
            mm2_bam = self.tmp_dir.joinpath(f"sample_mm2_{i}.bam")
            rg2 = copy.deepcopy(rg)
            rg2["ID"] = rg2["ID"] + "-mm2"
            rg2["LR"] = "1"

            mm2_model: Union[str, pathlib.Path] = self.model_bundle.joinpath(
                "minimap2.model"
            )
            mm2_xargs = []
            if self.tech.upper() == "ULTIMA":
                mm2_model = "sr"
                mm2_xargs = ["-2", "--secondary=yes"]
            mm2_job = Job(
                cmds.cmd_minimap2_lift(
                    mm2_bam,
                    sample_fasta,
                    fq,
                    sample_gfa,
                    self.reference,
                    "@RG\\t"
                    + "\\t".join([f"{x[0]}:{x[1]}" for x in rg2.items()]),
                    mm2_model,
                    threads=self.cores,
                    mm2_xargs=mm2_xargs,
                ),
                f"mm2-lift-{i}",
                self.cores,
            )

            mm2_bams.append(mm2_bam)
            mm2_jobs.append(mm2_job)
        return (mm2_bams, mm2_jobs)

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
            for rg in self.fastq_readgroups:
                read_filters.append(
                    f"IndelLeftAlignReadTransform,rgid={rg['ID']}-mm2"
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

    def build_dnascope_job(
        self,
        out_vcf: pathlib.Path,
        input_bams: List[pathlib.Path],
    ) -> Job:
        assert self.model_bundle

        read_filters = []
        if self.tech.upper() == "ULTIMA":
            read_filters.append("UltimaReadFilter")

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
                pcr_indel_model="NONE",
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
            # Extract contigs not in the pop vcf as merge will fail
            if shard.contig not in self.pop_vcf_contigs:
                if shard.contig in seen_contigs:
                    continue
                seen_contigs.add(shard.contig)
                subset_vcf = self.tmp_dir.joinpath(
                    f"sample-dnascope_transfer-subset{i}.vcf.gz"
                )
                view_job = Job(
                    cmds.cmd_bcftools_view_regions(
                        subset_vcf,
                        raw_vcf,
                        regions=shard.contig,
                    ),
                    "merge-trim-extra",
                    1,
                )
                sharded_merge_jobs.append(view_job)
                sharded_vcfs.append(subset_vcf)
            else:
                shard_vcf = self.tmp_dir.joinpath(
                    f"sample-dnascope_transfer-shard{i}.vcf.gz"
                )
                merge_job = Job(
                    cmds.cmd_bcftools_merge_trim(
                        shard_vcf,
                        raw_vcf,
                        self.pop_vcf,
                        trim_script,
                        str(shard),
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
