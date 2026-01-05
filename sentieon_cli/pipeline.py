"""
A pipeline class
"""

from abc import ABC, abstractmethod
import argparse
import multiprocessing as mp
import os
import pathlib
import shutil
import sys
from typing import Any, Dict, List, Optional

import packaging.version

from . import command_strings as cmds
from .dag import DAG
from .executor import BaseExecutor, DryRunExecutor, LocalExecutor
from .job import Job
from .logging import get_logger
from .scheduler import ThreadScheduler
from .util import __version__, check_version, path_arg, tmp

MULTIQC_MIN_VERSION = {
    "multiqc": packaging.version.Version("1.18"),
}


class BasePipeline(ABC):
    """A pipeline base class"""

    params: Dict[str, Dict[str, Any]] = {
        # Required arguments
        "reference": {
            "flags": ["-r", "--reference"],
            "required": True,
            "help": "fasta for reference genome.",
            "type": path_arg(exists=True, is_file=True),
        },
        # Additional arguments
        "cores": {
            "flags": ["-t", "--cores"],
            "help": (
                "Number of threads/processes to use. Defaults to all "
                "available."
            ),
            "default": mp.cpu_count(),
        },
        "dry_run": {
            "help": "Print the commands without running them.",
            "action": "store_true",
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

    @classmethod
    def add_arguments(cls, parser: argparse.ArgumentParser):
        for k, kwargs in cls.params.items():
            flags = ["--" + k]
            if "default" in kwargs and "type" not in kwargs:
                kwargs["type"] = type(kwargs["default"])
            if "flags" in kwargs:
                flags = kwargs["flags"]
                del kwargs["flags"]
            parser.add_argument(*flags, **kwargs)

        for k, kwargs in cls.positionals.items():
            parser.add_argument(k, **kwargs)

    def handle_arguments(self, args: argparse.Namespace):
        """Update self using the argparse object"""
        for k in self.params.keys():
            assert k in self.__dict__
            if k in args.__dict__:
                val = getattr(args, k)
                if val is not None:
                    setattr(self, k, val)
        for k in self.positionals.keys():
            assert k in self.__dict__
            if k in args.__dict__:
                setattr(self, k, getattr(args, k))

    def setup_logging(self, args: argparse.Namespace) -> None:
        self.logger = get_logger(__name__)
        assert self.logger.parent
        self.logger.parent.setLevel(args.loglevel)
        self.logger.info("Starting sentieon-cli version: %s", __version__)

    def __init__(self) -> None:
        self.reference: Optional[pathlib.Path] = None
        self.cores = mp.cpu_count()
        self.numa_nodes: List[str] = []
        self.dry_run = False
        self.retain_tmpdir = False
        self.skip_version_check = False
        self.output_vcf: Optional[pathlib.Path] = None

    def main(self, args: argparse.Namespace) -> None:
        """Run the DNAscope pipeline"""
        self.handle_arguments(args)
        self.setup_logging(args)
        self.validate()
        self.configure()

        tmp_dir_str = tmp()
        self.tmp_dir = pathlib.Path(tmp_dir_str)

        dag = self.build_dag()
        executor = self.run(dag)

        if not self.retain_tmpdir:
            shutil.rmtree(tmp_dir_str)

        self.check_execution(dag, executor)

    def check_execution(
        self,
        dag: DAG,
        executor: BaseExecutor,
    ):
        """Check the DAG and executor after a run"""
        if executor.jobs_with_errors:
            raise ValueError("Execution failed")

        if len(dag.waiting_jobs) > 0 or len(dag.ready_jobs) > 0:
            raise ValueError(
                "The DAG has some unexecuted jobs\n"
                f"Waiting jobs: {dag.waiting_jobs}\n"
                f"Ready jobs: {dag.ready_jobs}\n"
            )

    @abstractmethod
    def validate(self) -> None:
        pass

    def validate_output_vcf(self) -> None:
        if not str(self.output_vcf).endswith(".vcf.gz"):
            self.logger.error("The output file should end with '.vcf.gz'")
            sys.exit(2)

    def validate_ref(self) -> None:
        # Confirm the presence of the reference index file
        fai_file = str(self.reference) + ".fai"
        if not os.path.isfile(fai_file):
            self.logger.error(
                "Fasta index file %s does not exist. Please index the "
                "reference genome with 'samtools faidx'",
                fai_file,
            )
            sys.exit(2)

    @abstractmethod
    def configure(self) -> None:
        pass

    @abstractmethod
    def build_dag(self) -> DAG:
        pass

    def multiqc(self) -> Optional[Job]:
        """Run MultiQC on the metrics files"""

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

    def run(self, dag: DAG) -> BaseExecutor:
        """Execute the DAG"""
        self.logger.debug("Creating the scheduler")
        resources: Dict[str, int] = {}
        for i in range(len(self.numa_nodes)):
            resources[f"node{i}"] = 1
        scheduler = ThreadScheduler(
            dag,
            self.cores,
            resources,
        )

        self.logger.debug("Creating the executor")
        Executor = DryRunExecutor if self.dry_run else LocalExecutor
        executor = Executor(scheduler)

        self.logger.info("Starting execution")
        executor.execute()
        return executor
