"""
A base class for pangenome pipelines
"""

import copy
from enum import Enum
import json
import pathlib
import sys
from typing import List, Optional, Tuple

from importlib_resources import files

from . import command_strings as cmds
from .dag import DAG
from .driver import CNVscope, CNVModelApply, Driver
from .job import Job
from .pipeline import BasePipeline
from .shell_pipeline import Command, Pipeline
from .util import path_arg


class SampleSex(Enum):
    FEMALE = 1
    MALE = 2
    UNKNOWN = 3


class BasePangenome(BasePipeline):
    """A pipeline base class for short reads"""

    params = copy.deepcopy(BasePipeline.params)
    params.update(
        {
            # Required arguments
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
            # Additional arguments
            "bam_format": {
                "help": (
                    "Use the BAM format instead of CRAM for output aligned "
                    "files."
                ),
                "action": "store_true",
            },
            "dbsnp": {
                "flags": ["-d", "--dbsnp"],
                "help": (
                    "dbSNP vcf file Supplying this file will annotate "
                    "variants with their dbSNP refSNP ID numbers."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "expansion_catalog": {
                "help": (
                    "An ExpansionHunter variant catalog. Required for "
                    "expansion calling."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "kmer_memory": {
                "help": "Memory limit for KMC in GB.",
                "default": 30,
                "type": int,
            },
            "segdup_caller_genes": {
                "type": str,
                "help": (
                    "Genes for SegDup calling. Ex: "
                    "'CFH,CFHR3,CYP11B1,CYP2D6,GBA,NCF1,PMS2,SMN1,STRC'. "
                    "Required for SegDup calling."
                ),
            },
            "skip_cnv": {
                "help": "Skip CNV calling.",
                "action": "store_true",
            },
            "t1k_hla_seq": {
                "help": (
                    "The DNA HLA seq FASTA file for T1K. Required for HLA "
                    "calling."
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
                    "The DNA KIR seq FASTA file for T1K. Required for KIA "
                    "calling."
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
        }
    )

    positionals = BasePipeline.positionals

    def __init__(self) -> None:
        super().__init__()
        self.gbz: Optional[pathlib.Path] = None
        self.hapl: Optional[pathlib.Path] = None
        self.model_bundle: Optional[pathlib.Path] = None
        self.r1_fastq: List[pathlib.Path] = []
        self.r2_fastq: List[pathlib.Path] = []
        self.bam_format = False
        self.dbsnp: Optional[pathlib.Path] = None
        self.expansion_catalog: Optional[pathlib.Path] = None
        self.kmer_memory = 30
        self.segdup_caller_genes: Optional[str] = None
        self.skip_cnv = False
        self.t1k_hla_seq: Optional[pathlib.Path] = None
        self.t1k_hla_coord: Optional[pathlib.Path] = None
        self.t1k_kir_seq: Optional[pathlib.Path] = None
        self.t1k_kir_coord: Optional[pathlib.Path] = None

    def validate_t1k(self) -> None:
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

    def build_kmc_job(
        self, kmer_prefix: pathlib.Path, job_threads: int
    ) -> Job:
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
            job_threads,
        )

        return kmc_job

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

    def get_sex(self, ploidy_json: pathlib.Path) -> None:
        """Retrieve the sample sex"""
        if self.dry_run:
            self.logger.info("Setting sample sex to MALE for dry-run")
            self.sample_sex = SampleSex.MALE
            return
        with open(ploidy_json) as fh:
            data = json.load(fh)
            sex = data["sex"]
            if sex == "female":
                self.sample_sex = SampleSex.FEMALE
            elif sex == "male":
                self.sample_sex = SampleSex.MALE
            else:
                self.sample_sex = SampleSex.UNKNOWN

    def build_second_dag(self, input_aln: pathlib.Path) -> DAG:
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
        if not self.skip_cnv:
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
                input_aln,
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
                input_aln,
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
