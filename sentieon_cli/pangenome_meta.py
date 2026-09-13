"""
Read pangenome graph metadata without loading the graph.

A GBZ file is a `simple-sds` serialization of::

    GBZ header (16 B) | GBZ tags | GBWT | GBWTGraph sequences ...

and the GBWT inside it is::

    GBWT header (48 B) | GBWT tags | BWT | DA samples | Option<Metadata>

Every component carries its own length, so the multi-GB BWT body is
skipped with a seek and only the few MB holding the path names are read.
Parsing a whole-genome HPRC graph this way takes well under a second.

The metadata answers the two questions the pangenome pipelines ask before
they run:

* `pangenome_ref_name`, the sample name that `vg haplotypes
  --set-reference` and `vg convert -Q` expect. A wrong name is not an
  error for `vg`: it warns, writes a sampled graph with an empty
  `reference_samples` tag and no reference paths, and then converts a GFA
  with no `SN:Z:` rGFA tags. `pgutil lift` exits 0 with every read
  unmapped, so the run completes with silently ruined accuracy.
* `pangenome_contig_prefix`, the `<sample>#<phase>#` prefix of the
  reference path names, stripped by `pgutil lift` and by the `PangenomeSV`
  and `PGHapUpdate` algos.

The backbone reference is the sample in the `reference_samples` tag whose
paths are unfragmented: one path per contig, all with fragment 0. Any
other listed sample is chopped into fragments by the graph construction.

The module also runs as a script, for the run-time checks the pipelines
wire into their DAGs::

    python -m sentieon_cli.pangenome_meta detect --gbz graph.gbz
    python -m sentieon_cli.pangenome_meta check-gbz --gbz sample.gbz
        --reference_name GRCh38 --contig_prefix 'GRCh38#0#'
    python -m sentieon_cli.pangenome_meta check-gfa --gfa sample.gfa
        --reference_name GRCh38 --contig_prefix 'GRCh38#0#'
"""

from __future__ import annotations

import argparse
import pathlib
import struct
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import BinaryIO, Dict, List, Optional, Sequence, Tuple, Union

from .logging import get_logger

# Spelled out rather than taken from `__name__`: run as a script with
# `python -m sentieon_cli.pangenome_meta`, `__name__` is "__main__" and the
# logger would no longer be a child of the package logger, so its messages
# would miss the console handler and the run log.
logger = get_logger("sentieon_cli.pangenome_meta")

# Bytes in a `simple-sds` word
WORD = 8

# Magic numbers of the serialized structures
GBZ_TAG = 0x205A4247
GBWT_TAG = 0x6B376B37
METADATA_TAG = 0x6B375E7A
HAPL_MAGIC = 0x4C504148

# How much of a GFA `check_gfa` reads before giving up on the `SN:Z:` tags
DEFAULT_MAX_GFA_LINES = 100000

# A GBWT path name record: (sample, contig, phase, fragment)
PathName = Tuple[int, int, int, int]

PathLike = Union[str, pathlib.Path]


class PangenomeMetadataError(ValueError):
    """A pangenome file could not be parsed, or does not describe one
    unambiguous reference sample."""


def unpack_bits(words: bytes, count: int, width: int) -> List[int]:
    """Unpack `count` little-endian `width`-bit integers from 64-bit words.

    `simple-sds` packs an integer vector into a stream of 64-bit words with
    no padding, so an element may straddle a word boundary. The running
    buffer keeps the decode linear in `count`.
    """
    mask = (1 << width) - 1
    out: List[int] = []
    buf = 0
    nbuf = 0
    offset = 0
    while len(out) < count:
        while nbuf < width:
            end = offset + WORD
            buf |= int.from_bytes(words[offset:end], "little") << nbuf
            nbuf += 64
            offset = end
        out.append(buf & mask)
        buf >>= width
        nbuf -= width
    return out


class _Reader:
    """A cursor over the `simple-sds` serialization of a GBZ or GBWT.

    Every method reads (or skips) exactly one serialized structure and
    leaves the cursor on the next one.
    """

    def __init__(self, fh: BinaryIO, path: PathLike) -> None:
        self.fh = fh
        self.path = str(path)
        self.size = fh.seek(0, 2)
        fh.seek(0)

    def _fail(self, msg: str) -> PangenomeMetadataError:
        return PangenomeMetadataError(f"{self.path}: {msg}")

    def _read(self, nbytes: int) -> bytes:
        data = self.fh.read(nbytes)
        if len(data) != nbytes:
            raise self._fail("the file ends before its contents do")
        return data

    def _checked(self, nbytes: int) -> int:
        """Reject a length that cannot fit in the file"""
        if nbytes < 0 or nbytes > self.size:
            raise self._fail(
                f"implausible length {nbytes} in the file structure; the "
                "file may be truncated or may not be a pangenome file"
            )
        return nbytes

    # --- simple-sds primitives -------------------------------------------
    def u64(self) -> int:
        return int(struct.unpack("<Q", self._read(WORD))[0])

    def u32x2(self) -> Tuple[int, int]:
        low, high = struct.unpack("<II", self._read(WORD))
        return int(low), int(high)

    def skip(self, nbytes: int) -> None:
        self.fh.seek(self._checked(nbytes), 1)

    def vec_u8(self) -> bytes:
        count = self._checked(self.u64())
        data = self._read(count)
        self._read((-count) % WORD)  # padding to a whole word
        return data

    def skip_vec_u8(self) -> None:
        count = self.u64()
        self.skip(count + (-count) % WORD)

    def vec_u64_words(self) -> bytes:
        return self._read(self._checked(self.u64() * WORD))

    def skip_vec_u64(self) -> None:
        self.skip(self.u64() * WORD)

    def skip_option(self) -> None:
        """Skip an `Option<...>` stored as a word count plus its words"""
        self.skip(self.u64() * WORD)

    def raw_vector(self) -> Tuple[int, bytes]:
        nbits = self.u64()
        return nbits, self.vec_u64_words()

    def skip_raw_vector(self) -> None:
        self.u64()
        self.skip_vec_u64()

    def int_vector(self) -> List[int]:
        count = self._checked(self.u64())
        width = self.u64()
        if not 0 < width <= 64:
            raise self._fail(f"unsupported integer width {width}")
        nbits, words = self.raw_vector()
        if nbits != count * width:
            raise self._fail("inconsistent integer vector")
        return unpack_bits(words, count, width)

    def skip_int_vector(self) -> None:
        self.u64()
        self.u64()
        self.skip_raw_vector()

    def bit_vector(self) -> Tuple[int, bytes]:
        self.u64()  # the number of set bits
        nbits, words = self.raw_vector()
        self.skip_option()  # rank support
        self.skip_option()  # select support
        self.skip_option()  # select-zero support
        return nbits, words

    def skip_bit_vector(self) -> None:
        self.u64()
        self.skip_raw_vector()
        self.skip_option()
        self.skip_option()
        self.skip_option()

    def sparse_vector_ones(self) -> List[int]:
        """The positions of the set bits of an Elias-Fano sparse vector"""
        self.u64()  # universe size
        _nbits, words = self.bit_vector()
        count = self._checked(self.u64())
        width = self.u64()
        lbits, lwords = self.raw_vector()
        if lbits != count * width:
            raise self._fail("inconsistent sparse vector")
        low = unpack_bits(lwords, count, width)

        # The high part is unary-coded: the i-th set bit at position `p`
        # contributes the bucket `p - i`.
        result: List[int] = []
        index = 0
        for offset in range(0, len(words), WORD):
            end = offset + WORD
            word = int.from_bytes(words[offset:end], "little")
            base = offset * 8
            while word:
                lsb = word & -word
                position = base + lsb.bit_length() - 1
                word ^= lsb
                bucket = position - index
                result.append((bucket << width) | low[index])
                index += 1
                if index == count:
                    return result
        return result

    def skip_sparse_vector(self) -> None:
        self.u64()
        self.skip_bit_vector()
        self.skip_int_vector()

    # --- gbwt-rs structures ----------------------------------------------
    def string_array(self) -> List[str]:
        """A `StringArray`: offsets, a byte alphabet and packed indices"""
        starts = self.sparse_vector_ones()
        alphabet = self.vec_u8()
        packed = self.int_vector()
        try:
            data = bytes(alphabet[index] for index in packed)
        except IndexError as err:
            raise self._fail(
                "a string array indexes its alphabet out of range"
            ) from err
        starts.append(len(data))
        return [
            data[begin:end].decode("utf-8", "replace")
            for begin, end in zip(starts, starts[1:])
        ]

    def tags(self) -> Dict[str, str]:
        """A `Tags`: a string array of alternating keys and values"""
        array = self.string_array()
        return {array[2 * i]: array[2 * i + 1] for i in range(len(array) // 2)}

    def dictionary(self) -> List[str]:
        """A `Dictionary`: a string array plus a sorted-id index"""
        strings = self.string_array()
        self.skip_int_vector()  # sorted_ids
        return strings


@dataclass
class GBZMetadata:
    """The headers, tags and path names of a GBZ file"""

    gbz_version: int
    gbz_tags: Dict[str, str]
    gbwt_version: int
    gbwt_tags: Dict[str, str]
    sample_names: List[str]
    contig_names: List[str]
    path_names: List[PathName]
    sample_count: int = 0
    haplotype_count: int = 0
    contig_count: int = 0
    sequence_count: int = 0
    bytes_read: int = 0

    @property
    def reference_samples(self) -> List[str]:
        """The sample names listed by the `reference_samples` GBWT tag"""
        return self.reference_samples_tag.split()

    @property
    def reference_samples_tag(self) -> str:
        """The raw `reference_samples` GBWT tag, `''` when it is absent"""
        return self.gbwt_tags.get("reference_samples", "")

    def sample_of(self, path: PathName) -> str:
        """The sample name of a path record"""
        return self.sample_names[path[0]]

    def contig_of(self, path: PathName) -> str:
        """The contig name of a path record"""
        return self.contig_names[path[1]]

    def full_path_name(self, path: PathName) -> str:
        """The `<sample>#<phase>#<contig>` name of a path record"""
        return f"{self.sample_of(path)}#{path[2]}#{self.contig_of(path)}"


@dataclass
class HaplHeader:
    """The header of a `vg haplotypes` (.hapl) file"""

    version: int
    top_level_chains: int
    construction_jobs: int
    total_subchains: int
    total_kmers: int
    k: int
    tags: Dict[str, str] = field(default_factory=dict)


@dataclass
class PangenomeReference:
    """The backbone reference sample of a pangenome graph"""

    ref_name: str
    contig_prefix: str
    contigs: List[str]
    reference_samples: List[str]
    n_paths: int


@dataclass
class _SampleSummary:
    """What the path names say about one sample"""

    n_paths: int
    phases: List[int]
    contigs: List[str]
    fragmented: bool


def read_gbz_metadata(path: PathLike) -> GBZMetadata:
    """Read the tags and path names of a GBZ file.

    Only the file's headers, tags and metadata block are read; the BWT and
    the graph sequences are skipped with seeks.

    Raises:
        PangenomeMetadataError: the file is not a GBZ, is truncated, or
            carries no GBWT metadata.
    """
    try:
        with open(path, "rb") as fh:
            return _read_gbz_metadata(fh, path)
    except PangenomeMetadataError:
        raise
    except (struct.error, IndexError, ValueError, EOFError) as err:
        raise PangenomeMetadataError(
            f"{path}: the file could not be parsed as a GBZ pangenome "
            f"({err})"
        ) from err


def _read_gbz_metadata(fh: BinaryIO, path: PathLike) -> GBZMetadata:
    reader = _Reader(fh, path)

    tag, gbz_version = reader.u32x2()
    if tag != GBZ_TAG:
        raise PangenomeMetadataError(f"{path}: not a GBZ pangenome file")
    reader.u64()  # GBZ flags
    gbz_tags = reader.tags()

    tag, gbwt_version = reader.u32x2()
    if tag != GBWT_TAG:
        raise PangenomeMetadataError(
            f"{path}: the GBWT header was not found where the GBZ format "
            "puts it"
        )
    sequence_count = reader.u64()
    reader.u64()  # size
    reader.u64()  # offset
    reader.u64()  # alphabet size
    reader.u64()  # GBWT flags
    gbwt_tags = reader.tags()

    # The BWT body: a sparse-vector index and a (possibly compressed)
    # Vec<u8>. Both are skipped without reading their contents.
    reader.skip_sparse_vector()
    reader.skip_vec_u8()
    reader.skip_vec_u64()  # the document array samples

    # Option<Metadata>
    if reader.u64() == 0:
        raise PangenomeMetadataError(
            f"{path}: the GBWT carries no metadata, so its path names are "
            "unknown"
        )
    tag, _meta_version = reader.u32x2()
    if tag != METADATA_TAG:
        raise PangenomeMetadataError(
            f"{path}: the GBWT metadata header was not found"
        )
    sample_count = reader.u64()
    haplotype_count = reader.u64()
    contig_count = reader.u64()
    reader.u64()  # metadata flags
    n_paths = reader._checked(reader.u64())
    raw = reader._read(n_paths * 16)
    path_names: List[PathName] = [
        (rec[0], rec[1], rec[2], rec[3])
        for rec in struct.iter_unpack("<IIII", raw)
    ]
    sample_names = reader.dictionary()
    contig_names = reader.dictionary()

    return GBZMetadata(
        gbz_version=gbz_version,
        gbz_tags=gbz_tags,
        gbwt_version=gbwt_version,
        gbwt_tags=gbwt_tags,
        sample_names=sample_names,
        contig_names=contig_names,
        path_names=path_names,
        sample_count=sample_count,
        haplotype_count=haplotype_count,
        contig_count=contig_count,
        sequence_count=sequence_count,
        bytes_read=fh.tell(),
    )


def read_hapl_header(path: PathLike) -> HaplHeader:
    """Read the header of a `vg haplotypes` (.hapl) file.

    Raises:
        PangenomeMetadataError: the file is not a `vg haplotypes` file or
            is too short to hold a header.
    """
    try:
        with open(path, "rb") as fh:
            reader = _Reader(fh, path)
            magic, version = reader.u32x2()
            if magic != HAPL_MAGIC:
                raise PangenomeMetadataError(
                    f"{path}: not a vg haplotypes file"
                )
            return HaplHeader(
                version=version,
                top_level_chains=reader.u64(),
                construction_jobs=reader.u64(),
                total_subchains=reader.u64(),
                total_kmers=reader.u64(),
                k=reader.u64(),
                # Tags were added to the header in version 6
                tags=reader.tags() if version >= 6 else {},
            )
    except PangenomeMetadataError:
        raise
    except (struct.error, IndexError, ValueError, EOFError) as err:
        raise PangenomeMetadataError(
            f"{path}: the file could not be parsed as a vg haplotypes "
            f"file ({err})"
        ) from err


def _summarize_samples(meta: GBZMetadata) -> Dict[str, _SampleSummary]:
    """Group the path names of a GBZ by sample.

    A sample is "fragmented" when any of its paths carries a non-zero
    fragment, or when one (contig, phase) pair has more than one path.
    """
    by_sample: Dict[int, List[PathName]] = defaultdict(list)
    for path in meta.path_names:
        by_sample[path[0]].append(path)

    summaries: Dict[str, _SampleSummary] = {}
    for sample_id, paths in by_sample.items():
        if sample_id >= len(meta.sample_names):
            continue
        per_contig: Dict[Tuple[int, int], int] = defaultdict(int)
        for path in paths:
            per_contig[(path[1], path[2])] += 1
        summaries[meta.sample_names[sample_id]] = _SampleSummary(
            n_paths=len(paths),
            phases=sorted({path[2] for path in paths}),
            contigs=sorted({meta.contig_of(path) for path in paths}),
            fragmented=(
                any(path[3] != 0 for path in paths)
                or any(count > 1 for count in per_contig.values())
            ),
        )
    return summaries


def detect_pangenome_reference(gbz_path: PathLike) -> PangenomeReference:
    """Detect the backbone reference sample of a pangenome graph.

    The GBWT `reference_samples` tag lists the graph's reference samples,
    for example `"GRCh38 CHM13"`. Exactly one of them is the backbone: its
    paths are unfragmented, one per contig. The others were chopped into
    fragments when the graph was built.

    Raises:
        PangenomeMetadataError: the file is not a GBZ, has no metadata,
            names no reference sample, or names no single unfragmented
            one.
    """
    meta = read_gbz_metadata(gbz_path)
    ref_samples = meta.reference_samples
    if not ref_samples:
        raise PangenomeMetadataError(
            f"{gbz_path}: the GBWT 'reference_samples' tag is missing or "
            "empty, so the graph names no reference sample"
        )

    summaries = _summarize_samples(meta)
    present = [name for name in ref_samples if name in summaries]
    if not present:
        raise PangenomeMetadataError(
            f"{gbz_path}: no sample named by the 'reference_samples' tag "
            f"({' '.join(ref_samples)}) has any path in the graph"
        )

    backbone = [name for name in present if not summaries[name].fragmented]
    if not backbone:
        raise PangenomeMetadataError(
            f"{gbz_path}: every reference sample ({', '.join(present)}) is "
            "split into fragments, so the backbone reference cannot be "
            "identified"
        )
    if len(backbone) > 1:
        raise PangenomeMetadataError(
            f"{gbz_path}: more than one reference sample is unfragmented "
            f"({', '.join(backbone)}), so the backbone reference is "
            "ambiguous"
        )

    name = backbone[0]
    summary = summaries[name]
    if len(summary.phases) != 1:
        phases = ", ".join(str(phase) for phase in summary.phases)
        raise PangenomeMetadataError(
            f"{gbz_path}: the reference sample '{name}' has paths on more "
            f"than one haplotype ({phases}), so its contig prefix is "
            "ambiguous"
        )

    reference = PangenomeReference(
        ref_name=name,
        contig_prefix=f"{name}#{summary.phases[0]}#",
        contigs=summary.contigs,
        reference_samples=ref_samples,
        n_paths=summary.n_paths,
    )
    logger.debug(
        "Detected pangenome reference '%s' (contig prefix '%s', %d paths) "
        "in '%s'",
        reference.ref_name,
        reference.contig_prefix,
        reference.n_paths,
        gbz_path,
    )
    return reference


def check_sample_gbz(
    gbz_path: PathLike,
    ref_name: str,
    contig_prefix: str,
) -> None:
    """Check that a sampled pangenome kept its reference paths.

    `vg haplotypes --set-reference NAME` only warns when `NAME` is not a
    reference sample of the input graph. The sampled graph it writes then
    has an empty `reference_samples` tag and no reference paths at all,
    which every downstream step turns into empty output.

    Raises:
        PangenomeMetadataError: the file could not be read, its
            `reference_samples` tag is not exactly `ref_name`, or the
            reference sample has no path under `contig_prefix`.
    """
    meta = read_gbz_metadata(gbz_path)
    observed = meta.reference_samples
    if observed != [ref_name]:
        raise PangenomeMetadataError(
            f"{gbz_path}: the sampled pangenome has "
            f"reference_samples='{meta.reference_samples_tag}', expected "
            f"'{ref_name}'. The pangenome reference name may be wrong: "
            "`vg haplotypes --set-reference` drops every reference path "
            "when the name it is given is not a reference sample of the "
            "input graph."
        )

    reference_paths = [
        path for path in meta.path_names if meta.sample_of(path) == ref_name
    ]
    if not any(
        meta.full_path_name(path).startswith(contig_prefix)
        for path in reference_paths
    ):
        raise PangenomeMetadataError(
            f"{gbz_path}: none of the {len(reference_paths)} paths of the "
            f"reference sample '{ref_name}' starts with the contig prefix "
            f"'{contig_prefix}'. The pangenome reference name or contig "
            "prefix may be wrong."
        )
    logger.info(
        "The sampled pangenome '%s' has %d reference paths for '%s'",
        gbz_path,
        len(reference_paths),
        ref_name,
    )


def _gfa_header_reference_samples(header: str) -> List[str]:
    """The samples named by the `RS:Z:` tag of a GFA header line"""
    for gfa_field in header.rstrip("\n").split("\t")[1:]:
        if gfa_field.startswith("RS:Z:"):
            return gfa_field.removeprefix("RS:Z:").split()
    return []


def check_gfa(
    gfa_path: PathLike,
    ref_name: str,
    contig_prefix: str,
    max_lines: int = DEFAULT_MAX_GFA_LINES,
) -> None:
    """Check that a GFA carries the rGFA tags of the reference.

    `vg convert -f -Q NAME` writes an rGFA whose `S` lines carry
    `SN:Z:<sample>#<phase>#<contig>`, and lists the graph's remaining
    reference samples in the header's `RS:Z:` tag. When `NAME` does not
    match a reference sample, `vg` writes no `SN:Z:` tags at all and
    leaves `NAME` in `RS:Z:`; `pgutil lift` then leaves every read
    unmapped.

    Only the first `max_lines` lines are read: a whole-genome GFA runs to
    many GB, and its first segment lines already carry the tags.

    Raises:
        PangenomeMetadataError: the header lists `ref_name` among its
            reference samples, or no segment line in the head of the file
            carries an `SN:Z:<contig_prefix>` tag.
    """
    sn_tag = f"SN:Z:{contig_prefix}"
    header_checked = False
    lines_read = 0
    with open(gfa_path, "r") as fh:
        for line in fh:
            if lines_read >= max_lines:
                break
            lines_read += 1
            if line.startswith("H\t") and not header_checked:
                header_checked = True
                rs_samples = _gfa_header_reference_samples(line)
                if ref_name in rs_samples:
                    raise PangenomeMetadataError(
                        f"{gfa_path}: the GFA header still lists "
                        f"'{ref_name}' among its reference samples "
                        f"(RS:Z:{' '.join(rs_samples)}), so the rGFA tags "
                        "were written for a different reference. The "
                        "pangenome reference name may be wrong."
                    )
            elif line.startswith("S\t") and sn_tag in line:
                logger.info(
                    "The GFA '%s' carries '%s' rGFA tags", gfa_path, sn_tag
                )
                return

    raise PangenomeMetadataError(
        f"{gfa_path}: no segment line in the first {lines_read} lines "
        f"carries an '{sn_tag}' tag. `vg convert -Q` writes no rGFA tags "
        "when its reference name does not match a reference sample of the "
        "graph, and `pgutil lift` then leaves every read unmapped."
    )


def _build_parser() -> argparse.ArgumentParser:
    """The argument parser of the `pangenome_meta` script"""
    parser = argparse.ArgumentParser(
        prog="python -m sentieon_cli.pangenome_meta",
        description="Inspect and check pangenome graph metadata",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    detect = subparsers.add_parser(
        "detect",
        help="Print the reference sample detected in a pangenome GBZ",
    )
    detect.add_argument(
        "--gbz", required=True, type=pathlib.Path, help="The GBZ file"
    )

    check_gbz = subparsers.add_parser(
        "check-gbz",
        help="Check that a sampled GBZ kept its reference paths",
    )
    check_gbz.add_argument(
        "--gbz",
        required=True,
        type=pathlib.Path,
        help="The sampled GBZ file",
    )
    check_gbz.add_argument(
        "--reference_name",
        required=True,
        help="The expected pangenome reference sample name",
    )
    check_gbz.add_argument(
        "--contig_prefix",
        required=True,
        help="The expected pangenome contig prefix",
    )

    check_gfa_parser = subparsers.add_parser(
        "check-gfa",
        help="Check that a GFA carries the reference's rGFA tags",
    )
    check_gfa_parser.add_argument(
        "--gfa", required=True, type=pathlib.Path, help="The GFA file"
    )
    check_gfa_parser.add_argument(
        "--reference_name",
        required=True,
        help="The expected pangenome reference sample name",
    )
    check_gfa_parser.add_argument(
        "--contig_prefix",
        required=True,
        help="The expected pangenome contig prefix",
    )
    check_gfa_parser.add_argument(
        "--max_lines",
        type=int,
        default=DEFAULT_MAX_GFA_LINES,
        help=(
            "Stop looking for the rGFA tags after this many lines "
            f"(default: {DEFAULT_MAX_GFA_LINES})"
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run a pangenome metadata check; return 0 on success and 1 on
    failure"""
    args = _build_parser().parse_args(argv)
    try:
        if args.command == "detect":
            reference = detect_pangenome_reference(args.gbz)
            print(f"pangenome_ref_name\t{reference.ref_name}")
            print(f"pangenome_contig_prefix\t{reference.contig_prefix}")
            print(
                "reference_samples\t" + " ".join(reference.reference_samples)
            )
            print(f"n_paths\t{reference.n_paths}")
            print(f"n_contigs\t{len(reference.contigs)}")
            print("contigs\t" + ",".join(reference.contigs))
        elif args.command == "check-gbz":
            check_sample_gbz(args.gbz, args.reference_name, args.contig_prefix)
        else:
            check_gfa(
                args.gfa,
                args.reference_name,
                args.contig_prefix,
                max_lines=args.max_lines,
            )
    except (PangenomeMetadataError, OSError) as err:
        logger.error("%s", err)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
