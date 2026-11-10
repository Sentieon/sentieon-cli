"""
UltimaPE model bundle support
"""

from dataclasses import dataclass
import logging
import pathlib
import re
import sys
from typing import Dict, Iterable, List, Optional

import packaging.version

from .. import command_strings as cmds
from ..dag import DAG
from ..driver import LuaScript, RepeatModel
from ..job import Job
from ..util import require_versions
from .base import Stage, StageResult, driver_job

ULTIMA_PE_PLATFORM = "ULTIMAPE"

ULTIMA_PE_MIN_VERSIONS = {
    "sentieon driver": packaging.version.Version("202503.05.rc4"),
}

ULTIMA_PE_READ_FILTER = "UltimaReadFilter,remove_t0_tp=0,keep_oq=true"

ULTIMA_PE_BUNDLE_MEMBERS = (
    "dnascope.model",
    "read_tables.tsv",
    "ultima.model",
)

# Contigs excluded from the primary-contig index
NON_PRIMARY_CONTIG = re.compile(r"_alt$|^HLA-|_decoy$|^chrEBV$")

# Contigs of the hap-length table
HAPLEN_CONTIG = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$")


def validate_ultima_pe(
    logger: logging.Logger,
    bundle_members: Iterable[str],
    *,
    gvcf: bool,
    bed: Optional[pathlib.Path],
    skip_version_check: bool,
) -> None:
    """Check the run's inputs for an UltimaPE bundle"""
    members = set(bundle_members)
    missing = [m for m in ULTIMA_PE_BUNDLE_MEMBERS if m not in members]
    if missing:
        logger.error(
            "Expected model files not found in the model bundle file: %s",
            ", ".join(missing),
        )
        sys.exit(2)
    if gvcf:
        logger.error(
            "The `--gvcf` argument is not supported with an UltimaPE model "
            "bundle."
        )
        sys.exit(2)
    if bed is None:
        logger.warning(
            "A BED file is recommended with an UltimaPE model bundle. "
            "Without one, the repeat model also counts the alt and decoy "
            "contigs."
        )
    require_versions(ULTIMA_PE_MIN_VERSIONS, skip=skip_version_check)


def ultima_pe_repeat_model(output: pathlib.Path) -> RepeatModel:
    """RepeatModel for the UltimaPE DNAscope pass"""
    return RepeatModel(
        output,
        min_map_qual=1,
        min_group_count=5000,
        read_flag_mask="drop=supplementary",
        repeat_extension=5,
        max_repeat_unit_size=2,
        min_repeat_count=6,
    )


def ultima_pe_env(bundle: pathlib.Path) -> Dict[str, str]:
    """Environment for the UltimaPE DNAscope pass"""
    return {"HMGL_TABLES": str(bundle)}


def primary_contigs(fai_data: Dict[str, Dict[str, int]]) -> List[str]:
    """Primary contigs of the reference"""
    return [c for c in fai_data if not NON_PRIMARY_CONTIG.search(c)]


def haplen_intervals(fai_data: Dict[str, Dict[str, int]]) -> List[str]:
    """BED lines for the hap-length table"""
    return [
        f"{c}\t0\t{v['length']}\n"
        for c, v in fai_data.items()
        if HAPLEN_CONTIG.match(c)
    ]


@dataclass(kw_only=True)
class PrimaryIndexResult(StageResult):
    """The jobs and the index of `PrimaryIndexStage`"""

    index: pathlib.Path


@dataclass(kw_only=True)
class PrimaryIndexStage(Stage):
    """minimap2 index of the reference's primary contigs"""

    fai_data: Dict[str, Dict[str, int]]
    prefix: pathlib.Path

    def add_to(
        self, dag: DAG, upstream: Iterable[Job] = ()
    ) -> PrimaryIndexResult:
        deps = set(upstream)
        contigs = pathlib.Path(str(self.prefix) + ".contigs")
        fasta = pathlib.Path(str(self.prefix) + ".fa")
        index = pathlib.Path(str(self.prefix) + ".asm5.mmi")
        if not self.ctx.dry_run:
            contigs.write_text(
                "".join(f"{c}\n" for c in primary_contigs(self.fai_data))
            )

        fasta_job = Job(
            cmds.cmd_samtools_faidx_regions(
                fasta, self.ctx.reference, contigs
            ),
            "primary-fasta",
            0,
            task_name="hap-table",
        )
        index_job = Job(
            cmds.cmd_minimap2_index(
                index, fasta, "asm5", threads=self.ctx.cores
            ),
            "primary-index",
            0,
            task_name="hap-table",
        )
        dag.add_job(fasta_job, deps)
        dag.add_job(index_job, {fasta_job})
        return PrimaryIndexResult(
            jobs=[fasta_job, index_job],
            terminal={index_job},
            index=index,
        )


@dataclass(kw_only=True)
class HapTableResult(StageResult):
    """The jobs and the table of `HapTableStage`"""

    table: pathlib.Path


@dataclass(kw_only=True)
class HapTableStage(Stage):
    """Hap-length table from the sample's haplotypes"""

    script: pathlib.Path
    haplotypes: pathlib.Path
    index: pathlib.Path
    sample: str
    fai_data: Dict[str, Dict[str, int]]
    haps_bam: pathlib.Path
    table: pathlib.Path

    def add_to(self, dag: DAG, upstream: Iterable[Job] = ()) -> HapTableResult:
        deps = set(upstream)
        bed = pathlib.Path(str(self.table).removesuffix(".tsv.gz") + ".bed")
        if not self.ctx.dry_run:
            bed.write_text("".join(haplen_intervals(self.fai_data)))

        haps_job = Job(
            cmds.cmd_hapslice_minimap2(
                self.haps_bam,
                self.script,
                self.haplotypes,
                self.index,
                self.sample,
                threads=self.ctx.cores,
            ),
            "haps-bam",
            self.ctx.cores,
            task_name="hap-table",
        )
        haplen_job = driver_job(
            self.ctx,
            [LuaScript(self.script, [self.table], entry="haplen")],
            name="haplen",
            task_name="hap-table",
            inputs=[self.haps_bam],
            interval=bed,
        )
        dag.add_job(haps_job, deps)
        dag.add_job(haplen_job, {haps_job})
        return HapTableResult(
            jobs=[haps_job, haplen_job],
            terminal={haplen_job},
            table=self.table,
        )
