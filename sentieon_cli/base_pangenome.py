"""
A base class for pangenome pipelines
"""

import copy
import pathlib
import sys
from typing import List, Optional

from . import command_strings as cmds
from .job import Job
from .pipeline import BasePipeline
from .util import path_arg


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
            "kmer_memory": {
                "help": "Memory limit for KMC in GB.",
                "default": 30,
                "type": int,
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
        self.kmer_memory = 30
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
