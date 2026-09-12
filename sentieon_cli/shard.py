"""
Sharding functionality for Sentieon pipelines
"""

import pathlib
import re
from typing import Dict, Iterable, List, NamedTuple, Optional

from importlib.resources import files

from .logging import get_logger
from .util import VcfHeaderError, read_vcf_header

logger = get_logger(__name__)


# Signature contigs used to recognize the common human reference builds.
# A build is only recognized when every signature contig is present with
# the exact expected length.
BUILD_SIGNATURES: Dict[str, Dict[str, int]] = {
    "hg38": {
        "chr1": 248956422,
        "chr2": 242193529,
        "chrX": 156040895,
    },
    "hg19": {
        "chr1": 249250621,
        "chr2": 243199373,
        "chrX": 155270560,
    },
    "b37": {
        "1": 249250621,
        "2": 243199373,
        "X": 155270560,
    },
    "chm13": {
        "chr1": 248387328,
        "chr2": 242696752,
        "chrX": 154259566,
    },
}

DEFAULT_PLOIDY_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
B37_PLOIDY_CONTIGS = [str(i) for i in range(1, 23)] + ["X", "Y"]


class PloidyContigs(NamedTuple):
    """Contig names used for ploidy estimation.

    Fields are ``None`` when the defaults of `estimate_ploidy.py` apply.
    """

    contigs: Optional[List[str]] = None
    autosomes: Optional[List[str]] = None
    x_contig: Optional[str] = None
    y_contig: Optional[str] = None


class Shard(NamedTuple):
    contig: str
    start: int
    stop: int

    def __str__(self) -> str:
        return f"{self.contig}:{self.start}-{self.stop}"


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


def detect_reference_build(
    fai_data: Dict[str, Dict[str, int]],
) -> Optional[str]:
    """Identify the reference build from the fasta index.

    A build is reported only when all of its signature contigs are present
    with the expected lengths. Unrecognized references return `None`; the
    `--par_bed` argument covers those references.
    """
    for build, signature in BUILD_SIGNATURES.items():
        if all(
            fai_data.get(ctg, {}).get("length") == length
            for ctg, length in signature.items()
        ):
            logger.debug("Identified the reference build as '%s'", build)
            return build
    logger.debug("Could not identify the reference build")
    return None


def par_bed_for_build(build: Optional[str]) -> Optional[pathlib.Path]:
    """The packaged pseudo-autosomal region (PAR) BED for a build.

    Returns `None` when the build is unknown or the packaged BED file is
    absent. A zero-length data file is treated the same as a missing one,
    so the wiring stays safe until the BED file contents are shipped.
    """
    if build is None:
        return None

    try:
        par_bed = pathlib.Path(
            str(files("sentieon_cli.data").joinpath(f"par_{build}.bed"))
        )
    except ModuleNotFoundError:
        par_bed = None  # type: ignore[assignment]

    if par_bed is None or not par_bed.is_file() or par_bed.stat().st_size == 0:
        logger.warning(
            "No pseudo-autosomal region (PAR) BED file is packaged for the "
            "'%s' reference build.",
            build,
        )
        return None
    return par_bed


def ploidy_contigs_for_build(build: Optional[str]) -> PloidyContigs:
    """Contig names used for ploidy estimation with a reference build"""
    if build == "b37":
        return PloidyContigs(
            contigs=list(B37_PLOIDY_CONTIGS),
            autosomes=list(B37_PLOIDY_CONTIGS[:-2]),
            x_contig="X",
            y_contig="Y",
        )
    # The `estimate_ploidy.py` defaults are the chr-prefixed contig names
    return PloidyContigs()


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


def parse_vcf_contigs(
    header: Iterable[str],
) -> Dict[str, Optional[int]]:
    """The ID and length of every `##contig` line of a VCF header.

    A `##contig` line without a `length` maps to `None`.
    """
    kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
    contigs: Dict[str, Optional[int]] = {}
    for line in header:
        if not line.startswith("##contig"):
            continue
        s = line.index("<")
        e = line.index(">")
        d = dict(kvpat.findall(line[s + 1 : e]))  # noqa: E203
        ctg: str = d["ID"]
        length: Optional[str] = d.get("length", None)
        contigs[ctg] = int(length) if length else None
    return contigs


def vcf_contigs(in_vcf: pathlib.Path) -> Dict[str, Optional[int]]:
    """Report the contigs in the input VCF.

    Returns an empty dictionary, and logs why, when the header could not
    be read.
    """
    try:
        header = read_vcf_header(in_vcf)
    except (VcfHeaderError, OSError) as err:
        logger.error("%s", err)
        return {}
    return parse_vcf_contigs(header)
