"""
HLA and KIR genotyping with T1K
"""

from dataclasses import dataclass
import pathlib
from typing import Dict, Iterable, List, Optional

import packaging.version

from .. import command_strings as cmds
from ..dag import DAG
from ..job import Job
from .alignment import ReadWriterStage
from .base import Stage, StageResult

T1K_MIN_VERSIONS: Dict[str, Optional[packaging.version.Version]] = {
    "run-t1k": None,
}

DEFAULT_T1K_HLA_LOCUS = "chr6:28510020-33480577"
DEFAULT_T1K_KIR_LOCUS = "chr19:53100000-55800000"


@dataclass(kw_only=True)
class T1KResult(StageResult):
    """The jobs and outputs of `T1KStage`"""

    extract_job: Job
    t1k_job: Job
    extracted_bam: pathlib.Path
    output_dir: pathlib.Path


@dataclass(kw_only=True)
class T1KStage(Stage):
    """Extract the reads at a gene locus and genotype them with T1K.

    ``tag`` names the gene group ("hla" or "kir"); it labels the jobs and
    the extracted BAM. ``preset`` is T1K's own `--preset` ("hla-wgs" or
    "kir-wgs") and ``locus`` the region the reads are extracted from.
    """

    inputs: List[pathlib.Path]
    gene_seq: pathlib.Path
    gene_coord: pathlib.Path
    locus: str
    preset: str
    tag: str
    output_dir: pathlib.Path
    task_name: str = "t1k"

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> T1KResult:
        deps = set(upstream)

        # Extract reads overlapping the locus to a BAM file with ReadWriter
        extracted_bam = self.ctx.tmp_dir.joinpath(f"sample_{self.tag}.bam")
        extract_job = (
            ReadWriterStage(
                ctx=self.ctx,
                inputs=self.inputs,
                output=extracted_bam,
                name=f"t1k-{self.tag}-extract",
                task_name=self.task_name,
                interval=self.locus,
            )
            .build()
            .jobs[0]
        )
        dag.add_job(extract_job, deps)

        t1k_job = Job(
            cmds.cmd_t1k(
                self.output_dir,
                extracted_bam,
                gene_seq=self.gene_seq,
                gene_coord=self.gene_coord,
                preset=self.preset,
                threads=self.ctx.cores,
            ),
            f"t1k-{self.tag}",
            self.ctx.cores,
            task_name=self.task_name,
        )
        dag.add_job(t1k_job, {extract_job})

        return T1KResult(
            jobs=[extract_job, t1k_job],
            terminal={t1k_job},
            extract_job=extract_job,
            t1k_job=t1k_job,
            extracted_bam=extracted_bam,
            output_dir=self.output_dir,
        )
