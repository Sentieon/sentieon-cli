"""
Sharding functionality for Sentieon pipelines
"""

import pathlib
import re
import subprocess as sp
from typing import Dict, List, NamedTuple, Optional

from .logging import get_logger

logger = get_logger(__name__)


GRCH38_CONTIGS: Dict[str, int] = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}


class Shard(NamedTuple):
    contig: str
    start: int
    stop: int

    def __str__(self) -> str:
        return f"{self.contig}:{self.start}-{self.stop}"

    def bcftools_str(self) -> str:
        return f"{{{self.contig}}}:{self.start}-{self.stop}"


def parse_fai(ref_fai: pathlib.Path) -> Dict[str, Dict[str, int]]:
    """Parse a faidx index"""
    contigs: Dict[str, Dict[str, int]] = {}
    with open(ref_fai) as fh:
        for line in fh:
            try:
                chrom, length, offset, lb, lw = line.rstrip().split()
            except ValueError as err:
                logger.error(
                    "Reference fasta index (.fai) does not have the expected "
                    "format"
                )
                raise err
            contigs[chrom] = {
                "length": int(length),
                "offset": int(offset),
                "linebases": int(lb),
                "linewidth": int(lw),
            }
    return contigs


def determine_shards_from_fai(
    fai_data: Dict[str, Dict[str, int]], step: int
) -> List[Shard]:
    """Generate shards of the genome from the fasta index"""
    shards: List[Shard] = []
    for ctg, d in fai_data.items():
        pos = 1
        length = d["length"]
        while pos <= length:
            end = pos + step - 1
            end = end if end < length else length
            shards.append(Shard(ctg, pos, end))
            pos = end + 1
    return shards


def vcf_contigs(
    in_vcf: pathlib.Path, dry_run=False
) -> Dict[str, Optional[int]]:
    """Report the contigs in the input VCF"""
    if dry_run:
        return {
            "chr1": 100,
            "chr2": 200,
            "chr3": 300,
        }
    kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
    cmd = ["bcftools", "view", "-h", str(in_vcf)]
    p = sp.run(cmd, capture_output=True, text=True)
    contigs: Dict[str, Optional[int]] = {}
    for line in p.stdout.split("\n"):
        if not line.startswith("##contig"):
            continue
        s = line.index("<")
        e = line.index(">")
        d = dict(kvpat.findall(line[s + 1 : e]))  # noqa: E203
        ctg: str = d["ID"]
        length: Optional[str] = d.get("length", None)
        contigs[ctg] = int(length) if length else None
    return contigs
