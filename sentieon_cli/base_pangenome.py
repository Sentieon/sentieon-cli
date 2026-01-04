"""
A base class for pangenome pipelines
"""

import copy
from enum import Enum
import json
import pathlib
from typing import List, Optional

from importlib_resources import files

from . import command_strings as cmds
from .job import Job
from .pipeline import BasePipeline
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
