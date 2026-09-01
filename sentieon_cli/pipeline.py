"""
A pipeline class
"""

from abc import ABC, abstractmethod
import argparse
import json
import multiprocessing as mp
import os
import pathlib
import shutil
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

import packaging.version

from importlib.resources import files

from . import command_strings as cmds
from .dag import DAG
from .exceptions import DagExecutionError
from .executor import BaseExecutor, DryRunExecutor, LocalExecutor
from .job import Job
from .logging import get_logger, set_console_level
from .run_logs import RunLogs
from .scheduler import ThreadScheduler
from .shard import (
    PloidyContigs,
    detect_reference_build,
    par_bed_for_build,
)
from .util import (
    SampleSex,
    __version__,
    check_version,
    cnvscope_sex_args,
    path_arg,
    tmp,
)

MULTIQC_MIN_VERSION = {
    "multiqc": packaging.version.Version("1.18"),
}

BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")


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
        "log_dir": {
            "help": (
                "Directory for the run's log files. Defaults to the output "
                "VCF with the '.vcf.gz' suffix replaced by '_logs'."
            ),
            "type": path_arg(),
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
        # Build a fresh kwargs dict per argument so repeated calls do not
        # mutate the shared class-level params/positionals specs.
        for k, spec in cls.params.items():
            kwargs = dict(spec)
            flags = kwargs.pop("flags", ["--" + k])
            if "default" in kwargs and "type" not in kwargs:
                kwargs["type"] = type(kwargs["default"])
            parser.add_argument(*flags, **kwargs)

        for k, spec in cls.positionals.items():
            parser.add_argument(k, **dict(spec))

    def handle_arguments(self, args: argparse.Namespace):
        """Update self using the argparse object"""
        for k in self.params.keys():
            if k not in self.__dict__:
                raise ValueError(
                    f"Parameter '{k}' is not an attribute of "
                    f"{self.__class__.__name__}"
                )
            if k in args.__dict__:
                val = getattr(args, k)
                if val is not None:
                    setattr(self, k, val)
        for k in self.positionals.keys():
            if k not in self.__dict__:
                raise ValueError(
                    f"Positional '{k}' is not an attribute of "
                    f"{self.__class__.__name__}"
                )
            if k in args.__dict__:
                setattr(self, k, getattr(args, k))

    def __init__(self) -> None:
        self.reference: Optional[pathlib.Path] = None
        self.cores = mp.cpu_count()
        self.numa_nodes: List[str] = []
        self.dry_run = False
        self.retain_tmpdir = False
        self.skip_version_check = False
        self.log_dir: Optional[pathlib.Path] = None
        self.output_vcf: Optional[pathlib.Path] = None
        self.run_logs: Optional[RunLogs] = None
        self.sample_sex: Optional[SampleSex] = None
        self.cnv_par_bed: Optional[pathlib.Path] = None
        self.reference_build: Optional[str] = None

    def setup_logging(self, args: argparse.Namespace) -> None:
        """Configure console logging"""
        self.logger = get_logger(__name__)
        set_console_level(args.loglevel)

    def start_run_logs(self) -> None:
        """Create the run's log directory and start writing `run.log`.

        Called after `validate`, so a run rejected for an invalid output path
        neither creates directories nor clobbers a previous run's logs.
        """
        # File logging is skipped for dry-runs and when there is nothing to
        # derive a log directory from (a bare pipeline, as used by tests).
        if self.dry_run or (self.log_dir is None and self.output_vcf is None):
            self.logger.info("Starting sentieon-cli version: %s", __version__)
            return

        log_dir = self.log_dir
        if log_dir is None:
            # Defensive: `validate` has already checked the output path, but
            # not every pipeline's validation covers the suffix.
            self.validate_output_suffix()
            log_dir = pathlib.Path(
                str(self.output_vcf).removesuffix(".vcf.gz") + "_logs"
            )
        run_logs = RunLogs(log_dir)
        try:
            run_logs.setup()
        except OSError as exc:
            self.logger.error(
                "Could not prepare the log directory %s: %s", log_dir, exc
            )
            sys.exit(2)
        self.run_logs = run_logs

        # After the file handler is attached, so the banner reaches run.log
        self.logger.info("Starting sentieon-cli version: %s", __version__)
        self.logger.info("Writing logs to: %s", self.run_logs.log_dir)

    def log_completion(self, success: bool, start_time: float) -> None:
        """Report the outcome and duration of the run"""
        self.logger.info(
            "Finished sentieon-cli (status: %s, elapsed: %.1fs)",
            "succeeded" if success else "failed",
            time.monotonic() - start_time,
        )
        if not success and self.run_logs:
            self.logger.info(
                "Logs from this run are in: %s", self.run_logs.log_dir
            )

    def main(self, args: argparse.Namespace) -> None:
        """Run the DNAscope pipeline"""
        self.handle_arguments(args)
        self.setup_logging(args)
        start_time = time.monotonic()
        success = False
        try:
            self.validate()
            self.start_run_logs()
            self.configure()

            tmp_dir_str = tmp()
            self.tmp_dir = pathlib.Path(tmp_dir_str)

            try:
                dag = self.build_dag()
                executor = self.run(dag)
                self.check_execution(dag, executor)

                # Jobs that depend on the results of the first DAG, such
                # as sex-aware calling after ploidy estimation
                second_dag = self.build_second_dag()
                if second_dag is not None:
                    executor = self.run(second_dag)
                    self.check_execution(second_dag, executor)
            finally:
                if not self.retain_tmpdir:
                    shutil.rmtree(tmp_dir_str)
            success = True
        finally:
            self.log_completion(success, start_time)

    def check_execution(
        self,
        dag: DAG,
        executor: BaseExecutor,
    ):
        """Check the DAG and executor after a run"""
        if executor.jobs_with_errors:
            failed = ", ".join(str(job) for job in executor.jobs_with_errors)
            message = f"Execution failed for jobs: {failed}"
            if self.run_logs:
                message += f"\nTask logs are in: {self.run_logs.task_logs}"
            raise DagExecutionError(message)

        if len(dag.waiting_jobs) > 0 or len(dag.ready_jobs) > 0:
            raise DagExecutionError(
                "The DAG has some unexecuted jobs\n"
                f"Waiting jobs: {dag.waiting_jobs}\n"
                f"Ready jobs: {dag.ready_jobs}\n"
            )

    @abstractmethod
    def validate(self) -> None:
        pass

    def validate_output_suffix(self) -> None:
        """Confirm the output VCF file name ends in '.vcf.gz'"""
        if not str(self.output_vcf).endswith(".vcf.gz"):
            self.logger.error("The output file should end with '.vcf.gz'")
            sys.exit(2)

    def validate_output_vcf(self) -> None:
        self.validate_output_suffix()
        assert self.output_vcf is not None
        parent = self.output_vcf.resolve().parent
        if not parent.is_dir():
            self.logger.error(
                "The parent directory of the output VCF does not exist: %s",
                parent,
            )
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

    def validate_bwa_index(self) -> None:
        """Confirm the presence of BWA index files next to the reference."""
        missing = [
            str(self.reference) + suf
            for suf in BWA_INDEX_SUFFIXES
            if not os.path.isfile(str(self.reference) + suf)
        ]
        if missing:
            self.logger.error(
                "BWA index files are missing for reference %s: %s. Please "
                "index the reference with 'sentieon bwa index %s'.",
                self.reference,
                ", ".join(missing),
                self.reference,
            )
            sys.exit(2)

    @abstractmethod
    def configure(self) -> None:
        pass

    @abstractmethod
    def build_dag(self) -> DAG:
        pass

    def build_second_dag(self) -> Optional[DAG]:
        """Build a second DAG, run after the first one has finished.

        Pipelines with jobs that depend on the results of the first DAG,
        such as sex-aware calling after ploidy estimation, override this
        hook. Returning `None` runs a single DAG.
        """
        return None

    def build_ploidy_job(
        self,
        ploidy_json: pathlib.Path,
        deduped_bam: List[pathlib.Path],
        ploidy_contigs: Optional[PloidyContigs] = None,
    ) -> Job:
        """Estimate sample ploidy and sex"""
        estimate_ploidy = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("estimate_ploidy.py"))
        ).resolve()
        ploidy_contigs = ploidy_contigs or PloidyContigs()
        ploidy_job = Job(
            cmds.cmd_estimate_ploidy(
                ploidy_json,
                deduped_bam,
                estimate_ploidy,
                contigs=ploidy_contigs.contigs,
                autosomes=ploidy_contigs.autosomes,
                x_contig=ploidy_contigs.x_contig,
                y_contig=ploidy_contigs.y_contig,
            ),
            "estimate-ploidy",
            0,
            task_name="ploidy",
        )
        return ploidy_job

    def get_sex(self, ploidy_json: pathlib.Path) -> None:
        """Retrieve the sample sex"""
        if self.sample_sex is not None:
            # Supplied through `--sample_sex`
            return
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

    def resolve_cnv_par_bed(
        self,
        fai_data: Dict[str, Dict[str, int]],
        par_bed: Optional[pathlib.Path] = None,
        cnv_will_run: bool = True,
    ) -> None:
        """Identify the reference build and select the PAR BED file.

        `self.reference_build` is always identified, as the ploidy
        estimation contigs follow it. The PAR BED file is only used by
        CNV calling, so it is looked up only when CNVs will be called.
        """
        self.reference_build = detect_reference_build(fai_data)
        if not cnv_will_run:
            return
        if par_bed is not None:
            self.cnv_par_bed = par_bed
            return
        self.cnv_par_bed = par_bed_for_build(self.reference_build)

    def validate_cnv_par(self, cnv_will_run: bool) -> None:
        """Confirm a PAR BED file is available for CNV calling.

        The check runs during validation, before any job starts, and does
        not depend on the sample sex.
        """
        if not cnv_will_run or self.cnv_par_bed is not None:
            return

        self.logger.error(
            "CNV calling uses a BED file of the pseudo-autosomal regions "
            "(PAR), and no PAR BED file is available for this reference "
            "genome. Please supply the `--par_bed` argument."
        )
        sys.exit(2)

    def cnv_sex_args(self) -> Tuple[Optional[str], Optional[pathlib.Path]]:
        """The CNVscope `--sex` and `--par` arguments for this run"""
        return cnvscope_sex_args(self.sample_sex, self.cnv_par_bed)

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
            task_name="multiqc",
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
        executor: BaseExecutor
        if self.dry_run:
            executor = DryRunExecutor(scheduler)
        else:
            # Handle Ctrl-C/SIGTERM by terminating running jobs gracefully;
            # the handlers are installed only for the duration of the run.
            executor = LocalExecutor(
                scheduler,
                install_signal_handlers=True,
                run_logs=self.run_logs,
            )

        self.logger.info("Starting execution")
        executor.execute()
        return executor
