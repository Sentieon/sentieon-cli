"""
Annotation transfer from a population VCF
"""

from dataclasses import dataclass
import pathlib
import re
import subprocess as sp
import tempfile
from typing import Dict, Iterable, List, Optional, Set, Tuple

from importlib.resources import files

from .. import command_strings as cmds
from ..dag import DAG
from ..job import Job
from ..logging import get_logger
from ..shard import Shard
from .base import Stage, StageResult

logger = get_logger(__name__)


def build_transfer_jobs(
    out_vcf: pathlib.Path,
    pop_vcf: pathlib.Path,
    raw_vcf: pathlib.Path,
    base_tmp_dir: pathlib.Path,
    shards: List[Shard],
    pop_vcf_contigs: Dict[str, Optional[int]],
    fai_data: Dict[str, Dict[str, int]],
    dry_run: bool = False,
    cores: int = 1,
    *,
    tag: Optional[str] = None,
) -> Tuple[List[Job], Job]:
    """Transfer annotations from the pop_vcf to the raw_vcf.

    ``tag`` is inserted into the job names (``merge-trim-{tag}-...``) so a
    pipeline that transfers annotations more than once in a run can tell
    the two fan-outs apart.
    """

    name_base = "merge-trim" if tag is None else f"merge-trim-{tag}"

    # Get a unique tmpdir
    tmp_dir_str = tempfile.mkdtemp(dir=base_tmp_dir)
    tmp_dir = pathlib.Path(tmp_dir_str)

    # Generate merge rules from the population VCF
    merge_rules = "AC_v20:sum,AF_v20:sum,AC_genomes:sum,AF_genomes:sum"
    if not dry_run:
        kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
        cmd = ["bcftools", "view", "-h", str(pop_vcf)]
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
    for i, shard in enumerate(shards):
        # Use a BED file for unusual contig names
        subset_bed = tmp_dir.joinpath(
            f"sample-dnascope_transfer-subset{i}.bed"
        )

        # Extract contigs not in the pop vcf as merge will fail
        if shard.contig not in pop_vcf_contigs:
            if shard.contig in seen_contigs:
                continue
            logger.info("Skipping transfer for contig: %s", shard.contig)
            seen_contigs.add(shard.contig)
            subset_vcf = tmp_dir.joinpath(
                f"sample-dnascope_transfer-subset{i}.vcf.gz"
            )

            ctg_len = fai_data[shard.contig]["length"]
            if not dry_run:
                with open(subset_bed, "w") as fh:
                    print(f"{shard.contig}\t0\t{ctg_len}", file=fh)

            view_job = Job(
                cmds.cmd_bcftools_view_regions(
                    subset_vcf,
                    raw_vcf,
                    regions_file=subset_bed,
                ),
                f"{name_base}-extra",
                1,
                task_name="annotation-transfer",
            )
            sharded_merge_jobs.append(view_job)
            sharded_vcfs.append(subset_vcf)
        else:
            if not dry_run:
                with open(subset_bed, "w") as fh:
                    print(
                        f"{shard.contig}\t{shard.start}\t{shard.stop}",
                        file=fh,
                    )

            logger.debug("Transferring shard: %s", shard)
            shard_vcf = tmp_dir.joinpath(
                f"sample-dnascope_transfer-shard{i}.vcf.gz"
            )
            merge_job = Job(
                cmds.cmd_bcftools_merge_trim(
                    shard_vcf,
                    raw_vcf,
                    pop_vcf,
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
                f"{name_base}-{i}",
                1,
                task_name="annotation-transfer",
            )
            sharded_merge_jobs.append(merge_job)
            sharded_vcfs.append(shard_vcf)

    # Concat all shards
    concat_job = Job(
        cmds.bcftools_concat(
            out_vcf,
            sharded_vcfs,
            xargs=["--no-version", "--threads", str(cores)],
        ),
        f"{name_base}-concat",
        cores,
        task_name="annotation-transfer",
    )
    return (sharded_merge_jobs, concat_job)


@dataclass(frozen=True)
class TransferConfig:
    """The run-wide inputs the annotation transfer needs.

    Every pipeline that transfers annotations builds the same fan-out from
    the same four values, so they travel together instead of being threaded
    through each call site.
    """

    pop_vcf: pathlib.Path
    shards: List[Shard]
    pop_vcf_contigs: Dict[str, Optional[int]]
    fai_data: Dict[str, Dict[str, int]]


@dataclass(kw_only=True)
class TransferResult(StageResult):
    """The jobs and output of `TransferStage`"""

    shard_jobs: List[Job]
    concat_job: Job
    out_vcf: pathlib.Path


@dataclass(kw_only=True)
class TransferStage(Stage):
    """Transfer annotations from the population VCF onto a raw VCF.

    One job per genomic shard merges the shard's annotations, then a
    single concat job assembles ``out_vcf``. ``tag`` disambiguates the
    job names when a pipeline transfers more than once in a run.
    """

    config: TransferConfig
    raw_vcf: pathlib.Path
    out_vcf: pathlib.Path
    tag: Optional[str] = None

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> TransferResult:
        deps = set(upstream)

        shard_jobs, concat_job = build_transfer_jobs(
            self.out_vcf,
            self.config.pop_vcf,
            self.raw_vcf,
            self.ctx.tmp_dir,
            self.config.shards,
            self.config.pop_vcf_contigs,
            self.config.fai_data,
            self.ctx.dry_run,
            self.ctx.cores,
            tag=self.tag,
        )

        for job in shard_jobs:
            dag.add_job(job, deps)
        dag.add_job(concat_job, set(shard_jobs))

        return TransferResult(
            jobs=[*shard_jobs, concat_job],
            terminal={concat_job},
            shard_jobs=shard_jobs,
            concat_job=concat_job,
            out_vcf=self.out_vcf,
        )
