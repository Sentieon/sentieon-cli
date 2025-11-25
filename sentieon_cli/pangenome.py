"""
Pangenome alignment and variant calling pipeline
"""

import argparse
from enum import Enum
import itertools
import json
import multiprocessing as mp
import os
import pathlib
import shutil
import sys
from typing import Any, Dict, List, Optional, Tuple

import packaging.version

from importlib_resources import files

from . import command_strings as cmds
from .archive import ar_load
from .dag import DAG
from .driver import (
    AlignmentStat,
    BaseDistributionByCycle,
    CNVscope,
    CNVModelApply,
    Dedup,
    DNAModelApply,
    DNAscope,
    Driver,
    GCBias,
    InsertSizeMetricAlgo,
    MeanQualityByCycle,
    LocusCollector,
    ReadWriter,
    QualDistribution,
    WgsMetricsAlgo,
)
from .job import Job
from .pipeline import BasePipeline
from .shell_pipeline import Command, Pipeline
from .util import __version__, check_version, parse_rg_line, path_arg, tmp


PANGENOME_MIN_VERSIONS = {
    "kmc": None,
    "sentieon driver": packaging.version.Version("202503.01"),
    "vg": None,
    "bcftools": packaging.version.Version("1.10"),
    "samtools": packaging.version.Version("1.16"),
    "run-t1k": None,
    "ExpansionHunter": None,
    "segdup-caller": None,
}

MULTIQC_MIN_VERSION = {
    "multiqc": packaging.version.Version("1.18"),
}


SV_HDR_ATTR = {
    "INFO": {
        "AT": {
            "Number": "R",
            "Type": "String",
            "Description": '"Allele Traversal as path in graph"',
        },
        "DP": {
            "Number": "1",
            "Type": "Integer",
            "Description": '"Total Depth"',
        },
    },
    "FORMAT": {
        "AD": {
            "Number": ".",
            "Type": "Integer",
            "Description": (
                '"Allelic depths for the ref and alt alleles in the order '
                'listed"'
            ),
        },
        "DP": {
            "Number": "1",
            "Type": "Integer",
            "Description": '"Read Depth"',
        },
        "GL": {
            "Number": "G",
            "Type": "Float",
            "Description": (
                '"Genotype Likelihood, log10-scaled likelihoods of the data '
                "given the called genotype for each possible genotype "
                "generated from the reference and alternate alleles given the "
                'sample ploidy"'
            ),
        },
        "GP": {
            "Number": "1",
            "Type": "Float",
            "Description": (
                '"Genotype Probability, the log-scaled posterior probability '
                'of the called genotype"'
            ),
        },
        "GQ": {
            "Number": "1",
            "Type": "Integer",
            "Description": (
                '"Genotype Quality, the Phred-scaled probability estimate of '
                'the called genotype"'
            ),
        },
        "GT": {
            "Number": "1",
            "Type": "String",
            "Description": '"Genotype"',
        },
        "MAD": {
            "Number": "1",
            "Type": "Integer",
            "Description": '"Minimum site allele depth"',
        },
        "XD": {
            "Number": "1",
            "Type": "Float",
            "Description": (
                '"eXpected Depth, background coverage as used for the Poisson '
                'model"'
            ),
        },
    },
    "FILTER": {
        "lowad": {
            "Description": (
                '"Variant does not meet minimum allele read support threshold '
                'of 1"'
            ),
        },
        "lowdepth": {"Description": '"Variant has read depth less than 4"'},
        "PASS": {"Description": '"All filters passed"'},
    },
}


GRCh38_CONTIGS_NOT_IN_PANGENOME = {
    "chr14_KI270726v1_random",
    "chr15_KI270727v1_random",
    "chr1_KI270711v1_random",
    "chr22_KI270737v1_random",
    "chr22_KI270739v1_random",
    "chrUn_GL000213v1",
    "chrUn_KI270338v1",
    "chrUn_KI270364v1",
    "chrUn_KI270371v1",
    "chrUn_KI270372v1",
    "chrUn_KI270374v1",
    "chrUn_KI270375v1",
    "chrUn_KI270424v1",
    "chrUn_KI270528v1",
    "chrUn_KI270581v1",
    "chrUn_KI270587v1",
}


class SampleSex(Enum):
    FEMALE = 1
    MALE = 2
    UNKNOWN = 3


class PangenomePipeline(BasePipeline):
    """The Pangenome pipeline"""

    params: Dict[str, Dict[str, Any]] = {
        # Required arguments
        "reference": {
            "flags": ["-r", "--reference"],
            "required": True,
            "help": "fasta for reference genome.",
            "type": path_arg(exists=True, is_file=True),
        },
        "gbz": {
            "help": "The pangenome graph file in GBZ format.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "hapl": {
            "help": "The haplotype file.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "snarls": {
            "help": (
                "The vg snarls file for the .gbz file. Can be created "
                "using vg."
            ),
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "xg": {
            "help": "The xg file for the .gbz file. Can be created using vg",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "model_bundle": {
            "flags": ["-m", "--model_bundle"],
            "help": "The model bundle file.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "r1_fastq": {
            "nargs": "*",
            "help": "Sample R1 fastq files.",
            "type": path_arg(exists=True, is_file=True),
        },
        "r2_fastq": {
            "nargs": "*",
            "help": "Sample R2 fastq files.",
            "type": path_arg(exists=True, is_file=True),
        },
        "readgroups": {
            "nargs": "*",
            "help": (
                "Readgroup information for the fastq files. Only the ID and "
                "SM attributes are used."
            ),
        },
        # Additional arguments
        "bam_format": {
            "help": (
                "Use the BAM format instead of CRAM for output aligned files."
            ),
            "action": "store_true",
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
        "expansion_catalog": {
            "help": (
                "An ExpansionHunter variant catalog. Required for expansion "
                "calling."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "kmer_memory": {
            "help": "Memory limit for KMC in GB.",
            "default": 58,
            "type": int,
        },
        "segdup_caller_genes": {
            "type": str,
            "help": (
                "Genes for SegDup calling. Ex: "
                "'CFH,CFHR3,CYP11B1,CYP2D6,GBA,NCF1,PMS2,SMN1,STRC'. Required "
                "for SegDup calling."
            ),
        },
        "t1k_hla_seq": {
            "help": (
                "The DNA HLA seq FASTA file for T1K. Required for HLA calling."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "t1k_hla_coord": {
            "help": (
                "The DNA HLA coord FASTA file for T1K. Required for HLA "
                "calling."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "t1k_kir_seq": {
            "help": (
                "The DNA KIR seq FASTA file for T1K. Required for KIA calling."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "t1k_kir_coord": {
            "help": (
                "The DNA KIR coord FASTA file for T1K. Required for KIA "
                "calling."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        # Hidden arguments
        "retain_tmpdir": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
        },
        "skip_version_check": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
        },
    }

    positionals: Dict[str, Dict[str, Any]] = {
        "output_vcf": {
            "help": "Output VCF file. The file name must end in .vcf.gz",
            "type": path_arg(),
        },
    }

    def __init__(self) -> None:
        super().__init__()
        self.output_vcf: Optional[pathlib.Path] = None
        self.reference: Optional[pathlib.Path] = None
        self.gbz: Optional[pathlib.Path] = None
        self.hapl: Optional[pathlib.Path] = None
        self.snarls: Optional[pathlib.Path] = None
        self.xg: Optional[pathlib.Path] = None
        self.r1_fastq: List[pathlib.Path] = []
        self.r2_fastq: List[pathlib.Path] = []
        self.readgroups: List[str] = []
        self.model_bundle: Optional[pathlib.Path] = None
        self.bam_format = False
        self.dbsnp: Optional[pathlib.Path] = None
        self.cores = mp.cpu_count()
        self.expansion_catalog: Optional[pathlib.Path] = None
        self.kmer_memory = 128
        self.segdup_caller_genes: Optional[str] = None
        self.t1k_hla_seq: Optional[pathlib.Path] = None
        self.t1k_hla_coord: Optional[pathlib.Path] = None
        self.t1k_kir_seq: Optional[pathlib.Path] = None
        self.t1k_kir_coord: Optional[pathlib.Path] = None
        self.retain_tmpdir = False
        self.skip_version_check = False

    def main(self, args: argparse.Namespace) -> None:
        """Run the pipeline"""
        self.handle_arguments(args)
        self.setup_logging(args)
        self.validate()

        tmp_dir_str = tmp()
        self.tmp_dir = pathlib.Path(tmp_dir_str)

        dag = self.build_first_dag()
        executor = self.run(dag)
        self.check_execution(dag, executor)

        self.get_sex()
        dag = self.build_second_dag()
        executor = self.run(dag)
        self.check_execution(dag, executor)

        if not self.retain_tmpdir:
            shutil.rmtree(tmp_dir_str)

    def build_dag(self) -> DAG:
        return DAG()

    def validate(self) -> None:
        """Validate pipeline inputs"""
        assert self.output_vcf
        assert self.reference

        self.validate_bundle()

        if not self.r1_fastq:
            self.logger.error("Please supply --r1_fastq arguments")
            sys.exit(2)

        if not str(self.output_vcf).endswith(".vcf.gz"):
            self.logger.error("The output file should end with '.vcf.gz'")
            sys.exit(2)

        if len(self.r1_fastq) != len(self.readgroups):
            self.logger.error(
                "The number of readgroups does not equal the number of fastq "
                "files"
            )
            sys.exit(2)

        # Confirm the presence of the reference index file
        fai_file = str(self.reference) + ".fai"
        if not os.path.isfile(fai_file):
            self.logger.error(
                "Fasta index file %s does not exist. Please index the "
                "reference genome with 'samtools faidx'",
                fai_file,
            )
            sys.exit(2)

        # Validate readgroups
        rg_sample = None
        for rg in self.readgroups:
            rg_dict = parse_rg_line(rg.replace(r"\t", "\t"))
            rg_sm = rg_dict.get("SM")
            if not rg_sm:
                self.logger.error(
                    "Found a readgroup without a SM tag: %s",
                    str(rg),
                )
                sys.exit(2)
            if rg_sample and rg_sample != rg_sm:
                self.logger.error(
                    "Inconsistent readgroup sample information found in: %s",
                    str(rg),
                )
                sys.exit(2)
            rg_sample = rg_sm
            if "ID" not in rg_dict:
                self.logger.error(
                    "Found a readgroup without an ID tag: %s",
                    str(rg),
                )
                sys.exit(2)

        if not self.skip_version_check:
            for cmd, min_version in PANGENOME_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

        # Check the T1K files
        if (self.t1k_hla_seq and not self.t1k_hla_coord) or (
            self.t1k_hla_coord and not self.t1k_hla_seq
        ):
            self.logger.error(
                "For HLA calling, both the seq and coord fasta files need to "
                "be supplied. Exiting"
            )
            sys.exit(2)

        if (self.t1k_kir_seq and not self.t1k_kir_coord) or (
            self.t1k_kir_coord and not self.t1k_kir_seq
        ):
            self.logger.error(
                "For KIR calling, both the seq and coord fasta files need to "
                "be supplied. Exiting"
            )
            sys.exit(2)

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
        if bundle_pipeline != "DNAscope pangenome":
            self.logger.error("The model bundle is for a different pipeline.")
            sys.exit(2)

        bundle_members = set(ar_load(str(self.model_bundle)))
        if (
            "bwa.model" not in bundle_members
            or "cnv.model" not in bundle_members
            or "dnascope.model" not in bundle_members
        ):
            self.logger.error(
                "Expected model files not found in the model bundle file"
            )
            sys.exit(2)

    def configure(self) -> None:
        """Configure pipeline parameters"""
        pass

    def build_first_dag(self) -> DAG:
        """Build the first DAG for the pangenome pipeline"""
        self.logger.info("Building the first pangenome DAG")
        dag = DAG()

        # Find the sample name
        sample_name = parse_rg_line(self.readgroups[0].replace(r"\t", "\t"))[
            "SM"
        ]

        # Output files
        suffix = "bam" if self.bam_format else "cram"
        out_svs = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_svs.vcf.gz")
        )
        self.deduped_cram = pathlib.Path(
            str(self.output_vcf).replace(
                ".vcf.gz", f"_pangenome-aligned.{suffix}"
            )
        )
        self.ploidy_json = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_ploidy.json")
        )
        out_hla = pathlib.Path(str(self.output_vcf).replace(".vcf.gz", "_hla"))
        out_kir = pathlib.Path(str(self.output_vcf).replace(".vcf.gz", "_kir"))
        out_segdup = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_segdups")
        )

        # Intermediate file paths
        surject_paths_dict = self.tmp_dir.joinpath("GRCh38_pansn_paths.dict")
        sv_header = self.tmp_dir.joinpath("sample_sv.hdr")
        kmer_prefix = self.tmp_dir.joinpath("sample.fq")
        kmer_file = pathlib.Path(str(kmer_prefix) + ".kff")
        sample_pangenome = self.tmp_dir.joinpath("sample_pangenome.gbz")
        sample_pack = self.tmp_dir.joinpath("sample_coverage.pack")
        lc_score = self.tmp_dir.joinpath("sample_lc_score.txt.gz")
        # - Ensure we have a bam file for t1k
        rw_bam = self.tmp_dir.joinpath("sample_deduped.bam")
        dnascope_tmp = self.tmp_dir.joinpath("sample_dnascope_tmp.vcf.gz")

        # Create files used to generate headers
        self.create_paths_dict(surject_paths_dict)
        self.create_sv_header(sv_header, sample_name)

        # KMC k-mer counting
        kmc_job = self.build_kmc_job(kmer_prefix)
        dag.add_job(kmc_job)

        # vg haplotypes - create sample-specific pangenome
        haplotypes_job = self.build_haplotypes_job(sample_pangenome, kmer_file)
        dag.add_job(haplotypes_job, {kmc_job})

        # vg giraffe - align reads to pangenome
        sample_gams, giraffe_jobs = self.build_alignment_jobs(sample_pangenome)
        alignment_dependencies = {haplotypes_job}
        for job in giraffe_jobs:
            dag.add_job(job, alignment_dependencies)

        # vg pack - compute read support
        pack_job = self.build_pack_job(sample_pack, sample_gams)
        dag.add_job(pack_job, set(giraffe_jobs))

        # vg call - call SVs
        call_job = self.build_call_job(
            out_svs, sample_pack, sv_header, sample_name
        )
        dag.add_job(call_job, {pack_job})

        # surject and sort - move reads back to the linear reference
        aligned_bams, surject_jobs = self.build_surject_jobs(
            sample_gams, surject_paths_dict
        )
        for job in surject_jobs:
            dag.add_job(job, set(giraffe_jobs))

        # duplicate marking
        hdr_jobs, lc_job, dedup_job = self.build_dedup_job(
            lc_score, self.deduped_cram, aligned_bams
        )
        for job in hdr_jobs:
            dag.add_job(job, set(surject_jobs))
        dag.add_job(lc_job, set(hdr_jobs))
        dag.add_job(dedup_job, {lc_job})

        # metrics - collect some metrics after dedup
        metrics_job, rehead_job = self.build_metrics_job(
            rw_bam, self.deduped_cram
        )
        dag.add_job(metrics_job, {dedup_job})
        dag.add_job(rehead_job, {metrics_job})

        # multiqc
        multiqc_job = self.multiqc()
        if multiqc_job:
            dag.add_job(multiqc_job, {rehead_job})

        # small variant calling
        dnascope_job, dnamodelapply_job = self.build_smallvar_job(
            dnascope_tmp, self.deduped_cram
        )
        dag.add_job(dnascope_job, {dedup_job})
        dag.add_job(dnamodelapply_job, {dnascope_job})

        # estimate ploidy
        ploidy_job = self.build_ploidy_job(
            self.ploidy_json,
            [rw_bam],
        )
        dag.add_job(ploidy_job, {metrics_job})

        # HLA and KIR calling with T1K
        if self.t1k_hla_seq and self.t1k_hla_coord:
            hla_job = self.build_t1k_job(
                out_hla,
                rw_bam,
                self.t1k_hla_seq,
                self.t1k_hla_coord,
                "hla-wgs",
            )
            dag.add_job(hla_job, {metrics_job})
        if self.t1k_kir_seq and self.t1k_kir_coord:
            kir_job = self.build_t1k_job(
                out_kir,
                rw_bam,
                self.t1k_kir_seq,
                self.t1k_kir_coord,
                "kir-wgs",
            )
            dag.add_job(kir_job, {metrics_job})

        # special-caller
        if self.segdup_caller_genes:
            segdup_job = self.build_segdup_job(
                out_segdup, self.deduped_cram, genes=self.segdup_caller_genes
            )
            dag.add_job(segdup_job, {dedup_job})

        return dag

    def create_paths_dict(
        self,
        surject_paths_dict: pathlib.Path,
        ctg_prefix="GRCh38#0#",
        assembly="38",
        species="Homo sapiens",
    ):
        """Create a .dict file for surjection"""
        contigs: List[Tuple[str, str]] = []
        fai_file = str(self.reference) + ".fai"
        if not self.dry_run:
            with open(fai_file) as fh:
                for line in fh:
                    line_split = line.rstrip().split("\t")
                    contigs.append((line_split[0], line_split[1]))

            with open(surject_paths_dict, "w") as ofh:
                print("@HD\tVN:1.5", file=ofh)
                for ctg_name, ctg_len in contigs:
                    if ctg_name in GRCh38_CONTIGS_NOT_IN_PANGENOME:
                        continue
                    flds = {
                        "SN": f"{ctg_prefix}{ctg_name}",
                        "LN": ctg_len,
                        "AS": assembly,
                        "SP": species,
                    }
                    print(
                        "@SQ\t"
                        + "\t".join([f"{x}:{y}" for x, y in flds.items()]),
                        file=ofh,
                    )

    def create_sv_header(self, sv_header: pathlib.Path, sample_name: str):
        """Create a header file for the SV VCF"""
        contigs: List[Tuple[str, str]] = []
        fai_file = str(self.reference) + ".fai"
        if not self.dry_run:
            with open(fai_file) as fh:
                for line in fh:
                    line_split = line.rstrip().split("\t")
                    contigs.append((line_split[0], line_split[1]))

        if not self.dry_run:
            with open(sv_header, "w") as fh:
                print(r"##fileformat=VCFv4.2", file=fh)
                for contig, contig_l in contigs:
                    print(
                        f"##contig=<ID={contig},length={contig_l}>",
                        file=fh,
                    )

                for fld in ("FILTER", "INFO", "FORMAT"):
                    for hdr_id, hdr_d in SV_HDR_ATTR[fld].items():
                        hdr_kv = [f"ID={hdr_id}"]
                        for k, v in hdr_d.items():
                            hdr_kv.append(f"{k}={v}")
                        hdr_kv_str = ",".join(hdr_kv)
                        print(f"##{fld}=<{hdr_kv_str}>", file=fh)

                last_hdr_line = [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    sample_name,
                ]
                print("\t".join(last_hdr_line), file=fh)

    def build_kmc_job(self, kmer_prefix: pathlib.Path) -> Job:
        """Build KMC k-mer counting jobs"""
        # Create file list for KMC
        file_list = pathlib.Path(str(kmer_prefix) + ".paths")
        all_fastqs = []

        # Add R1 files
        all_fastqs.extend(self.r1_fastq)

        # Add R2 files if present
        if self.r2_fastq:
            all_fastqs.extend(self.r2_fastq)

        # Write file list
        if not self.dry_run:
            with open(file_list, "w") as f:
                for fq in all_fastqs:
                    f.write(f"{fq}\n")

        # Create KMC job
        kmc_job = Job(
            cmds.cmd_kmc(
                kmer_prefix,
                file_list,
                self.tmp_dir,
                memory=self.kmer_memory,
                threads=self.cores,
            ),
            "kmc",
            self.cores,
        )

        return kmc_job

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
            ),
            "vg-haplotypes",
            self.cores,
        )
        return haplotypes_job

    def build_alignment_jobs(
        self, sample_pangenome: pathlib.Path
    ) -> Tuple[List[pathlib.Path], List[Job]]:
        """Build vg giraffe alignment jobs"""
        # Process each pair of fastq files
        sample_gams: List[pathlib.Path] = []
        giraffe_jobs: List[Job] = []
        for i, (r1, r2, rg) in enumerate(
            itertools.zip_longest(
                self.r1_fastq, self.r2_fastq, self.readgroups
            )
        ):
            sample_gam = self.tmp_dir.joinpath(f"sample_aligned_{i}.gam")
            rg_dict = parse_rg_line(rg.replace(r"\t", "\t"))
            giraffe_job = Job(
                cmds.cmd_vg_giraffe(
                    sample_gam,
                    sample_pangenome,
                    fastq1=r1,
                    fastq2=r2,
                    readgroup=rg_dict["ID"],
                    sample=rg_dict["SM"],
                    threads=self.cores,
                ),
                "vg-giraffe",
                self.cores,
            )

            sample_gams.append(sample_gam)
            giraffe_jobs.append(giraffe_job)
        return (sample_gams, giraffe_jobs)

    def build_pack_job(
        self,
        sample_pack: pathlib.Path,
        sample_gams: List[pathlib.Path],
    ) -> Job:
        """Compute read support with vg pack"""
        assert self.gbz

        pack_job = Job(
            cmds.cmd_vg_pack(
                sample_pack,
                sample_gams,
                self.gbz,
                threads=self.cores,
            ),
            "vg-pack",
            self.cores,
        )
        return pack_job

    def build_call_job(
        self,
        out_svs: pathlib.Path,
        sample_pack: pathlib.Path,
        sv_header: pathlib.Path,
        sample_name="",
    ) -> Job:
        """Call SVs with vg call"""
        assert self.gbz
        assert self.snarls

        call_job = Job(
            cmds.cmd_sv_call(
                out_svs,
                sample_pack,
                sv_header,
                self.gbz,
                self.snarls,
                sample_name=sample_name,
            ),
            "vg-call",
            self.cores,
        )
        return call_job

    def build_surject_jobs(
        self,
        sample_gams: List[pathlib.Path],
        surject_paths_dict: pathlib.Path,
    ) -> Tuple[List[pathlib.Path], List[Job]]:
        """Surject gam files on the linear reference"""
        assert self.xg

        aligned_bams: List[pathlib.Path] = []
        surject_jobs: List[Job] = []

        for i, sample_gam in enumerate(sample_gams):
            sample_bam = self.tmp_dir.joinpath(f"sample_aligned_{i}.bam")
            surject_job = Job(
                cmds.cmd_vg_surject(
                    sample_bam,
                    sample_gam,
                    self.xg,
                    surject_paths_dict,
                    threads=self.cores,
                ),
                f"vg_surject_{i}",
                self.cores,
            )

            aligned_bams.append(sample_bam)
            surject_jobs.append(surject_job)
        return (aligned_bams, surject_jobs)

    def build_dedup_job(
        self,
        lc_score: pathlib.Path,
        deduped_bam: pathlib.Path,
        aligned_bams: List[pathlib.Path],
    ) -> Tuple[List[Job], Job, Job]:
        """Mark duplicates"""
        assert self.output_vcf

        # Create .hdr files for each input file
        hdr_jobs: List[Job] = []
        for i, bam in enumerate(aligned_bams):
            hdr = pathlib.Path(str(bam) + ".hdr")
            hdr_jobs.append(
                Job(
                    cmds.strip_ctg_prefix(
                        hdr,
                        bam,
                        "GRCh38#0#",
                    ),
                    f"create_hdr_{i}",
                    0,
                )
            )

        # Output metrics
        metric_base = self.output_vcf.name.replace(".vcf.gz", "") + ".txt"
        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )
        try:
            os.mkdir(metrics_dir)
        except FileExistsError:
            pass
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
        gc_metrics = metrics_dir.joinpath(metric_base + ".gc_bias.txt")
        gc_summary = metrics_dir.joinpath(metric_base + ".gc_bias_summary.txt")
        dedup_metrics = metrics_dir.joinpath(
            metric_base + ".dedup_metrics.txt"
        )

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=aligned_bams,
        )
        driver.add_algo(LocusCollector(lc_score))
        driver.add_algo(MeanQualityByCycle(mqbc_metrics))
        driver.add_algo(BaseDistributionByCycle(bdbc_metrics))
        driver.add_algo(QualDistribution(qualdist_metrics))
        driver.add_algo(AlignmentStat(as_metrics))
        driver.add_algo(GCBias(gc_metrics, summary=gc_summary))

        lc_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "locuscollector",
            self.cores,
        )

        # Dedup
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=aligned_bams,
            read_filter="IndelLeftAlignReadTransform",
        )
        driver.add_algo(
            Dedup(
                deduped_bam,
                lc_score,
                cram_write_options="version=3.0,compressor=rans",
                metrics=dedup_metrics,
            )
        )
        dedup_job = Job(
            Pipeline(Command(*driver.build_cmd())), "dedup", self.cores
        )

        return (hdr_jobs, lc_job, dedup_job)

    def build_metrics_job(
        self,
        deduped_bam: pathlib.Path,
        deduped_cram: pathlib.Path,
    ) -> Tuple[Job, Job]:
        assert self.output_vcf

        # Output metrics
        metric_base = self.output_vcf.name.replace(".vcf.gz", "") + ".txt"
        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )

        # Output metrics
        # - Run some metrics after duplicate marking
        is_metrics = metrics_dir.joinpath(metric_base + ".insert_size.txt")
        wgs_metrics = metrics_dir.joinpath(metric_base + ".wgs.txt")

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[deduped_cram],
        )
        driver.add_algo(ReadWriter(deduped_bam))
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        driver.add_algo(WgsMetricsAlgo(wgs_metrics, include_unpaired="true"))
        metrics_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "metrics",
            self.cores,
        )

        # Rehead WGS metrics so they are recognized by MultiQC
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

    def multiqc(self) -> Optional[Job]:
        """Run MultiQC on the metrics files"""
        assert self.output_vcf
        if not self.skip_version_check:
            if not all(
                [
                    check_version(cmd, min_version)
                    for (cmd, min_version) in MULTIQC_MIN_VERSION.items()
                ]
            ):
                self.logger.warning(
                    "Skipping MultiQC. MultiQC version %s or later not found",
                    MULTIQC_MIN_VERSION["multiqc"],
                )
                return None

        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )
        multiqc_job = Job(
            cmds.cmd_multiqc(
                metrics_dir,
                metrics_dir,
                f"Generated by the Sentieon-CLI version {__version__}",
            ),
            "multiqc",
            0,
        )
        return multiqc_job

    def build_smallvar_job(
        self,
        dnascope_tmp: pathlib.Path,
        realigned_cram: pathlib.Path,
    ) -> Tuple[Job, Job]:
        """Small variant calling"""
        assert self.model_bundle
        assert self.output_vcf

        model = self.model_bundle.joinpath("dnascope.model")
        interval_ctgs = list(["chr" + str(x) for x in range(1, 23)])
        interval_ctgs.extend(["chrX", "chrY"])
        interval = ",".join(interval_ctgs)

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[realigned_cram],
            interval=interval,
        )
        driver.add_algo(
            DNAscope(
                dnascope_tmp,
                dbsnp=self.dbsnp,
                pcr_indel_model="NONE",
                model=model,
            )
        )
        dnascope_job = Job(
            Pipeline(Command(*driver.build_cmd())), "dnascope", self.cores
        )

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            DNAModelApply(
                model,
                dnascope_tmp,
                self.output_vcf,
            )
        )
        dnamodelapply_job = Job(
            Pipeline(Command(*driver.build_cmd())), "model-apply", self.cores
        )

        return (dnascope_job, dnamodelapply_job)

    def build_ploidy_job(
        self,
        ploidy_json: pathlib.Path,
        deduped_bam: List[pathlib.Path],
    ) -> Job:
        """Estimate sample ploidy and sex"""
        estimate_ploidy = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("estimate_ploidy.py"))
        ).resolve()
        ploidy_job = Job(
            cmds.cmd_estimate_ploidy(
                ploidy_json,
                deduped_bam,
                estimate_ploidy,
            ),
            "estimate-ploidy",
            0,
        )
        return ploidy_job

    def get_sex(self) -> None:
        """Retrieve the sample sex"""
        if self.dry_run:
            self.logger.info("Setting sample sex to MALE for dry-run")
            self.sample_sex = SampleSex.MALE
            return
        with open(self.ploidy_json) as fh:
            data = json.load(fh)
            sex = data["sex"]
            if sex == "female":
                self.sample_sex = SampleSex.FEMALE
            elif sex == "male":
                self.sample_sex = SampleSex.MALE
            else:
                self.sample_sex = SampleSex.UNKNOWN

    def build_t1k_job(
        self,
        out_basename: pathlib.Path,
        deduped_bam: pathlib.Path,
        gene_seq: pathlib.Path,
        gene_coord: pathlib.Path,
        preset: str,
    ) -> Job:
        """Genotype diverse sequences with T1K"""
        job = Job(
            cmds.cmd_t1k(
                out_basename,
                deduped_bam,
                gene_seq=gene_seq,
                gene_coord=gene_coord,
                preset=preset,
                threads=1,
            ),
            "T1K-HLA-KIR",
            0,
        )
        return job

    def build_segdup_job(
        self,
        out_segdup: pathlib.Path,
        realigned_cram: pathlib.Path,
        genes: str,
    ) -> Job:
        """Call variants in difficult SegDups"""
        assert self.reference
        assert self.model_bundle

        segdup_job = Job(
            cmds.cmd_segdup_caller(
                out_segdup,
                realigned_cram,
                reference=self.reference,
                sr_bundle=self.model_bundle,
                genes=genes,
            ),
            "segdup-caller",
            self.cores,
        )
        return segdup_job

    def build_second_dag(self) -> DAG:
        """Build the second DAG for the pangenome pipeline"""
        self.logger.info("Building the second pangenome DAG")
        dag = DAG()

        # Output files
        out_cnv = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_cnv.vcf.gz")
        )
        out_expansions = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_expansion")
        )

        # Intermediate file paths
        cnvscope_tmp = self.tmp_dir.joinpath("sample_cnvcope_tmp.vcf.gz")
        cnvscope_autosomes = self.tmp_dir.joinpath(
            "sample_cnvcope_autosomes.vcf.gz"
        )
        cnvscope_haploid_tmp = self.tmp_dir.joinpath(
            "sample_cnvcope_haploid_tmp.vcf.gz"
        )
        cnvscope_haploid = self.tmp_dir.joinpath(
            "sample_cnvcope_haploid.vcf.gz"
        )

        # CNV calling
        (
            cnvscope_job,
            cnvapply_job,
            haploid_cnv_job,
            haploid_apply_job,
            concat_job,
        ) = self.build_cnv_jobs(
            out_cnv,
            cnvscope_tmp,
            cnvscope_autosomes,
            cnvscope_haploid_tmp,
            cnvscope_haploid,
            self.deduped_cram,
        )
        dag.add_job(cnvscope_job)
        dag.add_job(cnvapply_job, {cnvscope_job})
        if haploid_cnv_job and haploid_apply_job and concat_job:
            dag.add_job(haploid_cnv_job)
            dag.add_job(haploid_apply_job, {haploid_cnv_job})
            dag.add_job(concat_job, {cnvapply_job, haploid_apply_job})

        # Repeat expansions
        if self.expansion_catalog:
            expansion_job = self.build_expansion_job(
                out_expansions,
                self.deduped_cram,
                self.expansion_catalog,
            )
            dag.add_job(expansion_job)

        return dag

    def build_cnv_jobs(
        self,
        out_cnv: pathlib.Path,
        cnvscope_tmp: pathlib.Path,
        cnvscope_autosomes: pathlib.Path,
        cnvscope_haploid_tmp: pathlib.Path,
        cnvscope_haploid: pathlib.Path,
        realigned_cram: pathlib.Path,
    ) -> Tuple[Job, Job, Optional[Job], Optional[Job], Optional[Job]]:
        """Call CNVs with CNVscope"""
        assert self.model_bundle

        model = self.model_bundle.joinpath("cnv.model")

        # First pass will be autosomes for males
        # or autosomes + chrX for females
        first_contigs = list(["chr" + str(x) for x in range(1, 23)])
        if self.sample_sex != SampleSex.MALE:
            first_contigs.append("chrX")

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[realigned_cram],
            interval=",".join(first_contigs),
        )
        driver.add_algo(
            CNVscope(
                cnvscope_tmp,
                model=model,
            )
        )
        cnvscope_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "cnvscope-autosomes",
            self.cores,
        )

        # Autosome model apply
        outfile = (
            cnvscope_autosomes
            if self.sample_sex == SampleSex.MALE
            else out_cnv
        )
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            CNVModelApply(
                outfile,
                model,
                cnvscope_tmp,
            )
        )
        cnvapply_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "cnvmodelapply-autosomes",
            self.cores,
        )

        if self.sample_sex != SampleSex.MALE:
            return (cnvscope_job, cnvapply_job, None, None, None)

        # Sex chrom CNVscope
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[realigned_cram],
            interval="chrX,chrY",
        )
        driver.add_algo(
            CNVscope(
                cnvscope_haploid_tmp,
                model=model,
            )
        )
        haploid_cnv_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "cnvscope-haploid",
            self.cores,
        )

        # Sex chrom CNVModelApply
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            CNVModelApply(
                cnvscope_haploid,
                model,
                cnvscope_haploid_tmp,
            )
        )
        haploid_apply_job = Job(
            Pipeline(Command(*driver.build_cmd())),
            "cnvmodelapply-haploid",
            self.cores,
        )

        # Concate autosome and sex chromsome VCFs
        concat_job = Job(
            cmds.bcftools_concat(
                out_cnv, [cnvscope_autosomes, cnvscope_haploid]
            ),
            "concat-calls",
            0,
        )

        return (
            cnvscope_job,
            cnvapply_job,
            haploid_cnv_job,
            haploid_apply_job,
            concat_job,
        )

    def build_expansion_job(
        self,
        out_expansions: pathlib.Path,
        realigned_cram: pathlib.Path,
        expansion_catalog: pathlib.Path,
    ) -> Job:
        """Identify repeat expansions"""
        assert self.reference

        expansion_job = Job(
            cmds.cmd_expansion_hunter(
                out_expansions,
                realigned_cram,
                reference=self.reference,
                variant_catalog=expansion_catalog,
                sex="male" if self.sample_sex == SampleSex.MALE else "female",
                threads=self.cores,
            ),
            "expansion-hunter",
            0,
        )
        return expansion_job
