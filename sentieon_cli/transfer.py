"""
Annotation transfer functionality
"""

import pathlib
import re
import subprocess as sp
from typing import Any, Dict, List, Optional, Set, Tuple

from importlib_resources import files

from . import command_strings as cmds
from .job import Job
from .logging import get_logger
from .shard import Shard

logger = get_logger(__name__)


def build_transfer_jobs(
    out_vcf: pathlib.Path,
    pop_vcf: pathlib.Path,
    raw_vcf: pathlib.Path,
    tmp_dir: pathlib.Path,
    shards: List[Shard],
    pop_vcf_contigs: Dict[str, Optional[int]],
    fai_data: Dict[str, Dict[str, int]],
    dry_run: bool = False,
    cores: int = 1,
) -> Tuple[List[Job], Job]:
    """Transfer annotations from the pop_vcf to the raw_vcf"""

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
                "merge-trim-extra",
                1,
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
            xargs=["--no-version", "--threads", str(cores)],
        ),
        "merge-trim-concat",
        cores,
    )
    return (sharded_merge_jobs, concat_job)
