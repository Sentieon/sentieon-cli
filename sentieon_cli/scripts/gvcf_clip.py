#!/usr/bin/env python
"""Restrict a gVCF to the inside, or the outside, of a set of BED regions

The DNAscope hybrid pipeline splices its two calling passes together: the
first pass covers the positions outside the realigned regions and the
second pass the positions inside them. Variant records are kept or dropped
by their start position, as with `bcftools view -T`. Reference blocks are
clipped at the region boundaries. A dropped variant record can span
positions on the kept side, such as a deletion that starts inside a region
and ends past it; a reference block with no genotype confidence covers
whichever of those positions no kept record covers. So the spliced gVCF
has neither overlapping blocks nor gaps.

An output of `-` writes uncompressed VCF to stdout.
"""

import argparse
import bisect
import collections
import functools
import io
import multiprocessing as mp
import re
import sys

import vcflib

REF_BLOCK_ALTS = (["<NON_REF>"], ["<*>"])
# The FORMAT fields of a fill block, in this order, when the header has them
FILL_KEYS = ("GT", "DP", "GQ", "MIN_DP", "PL")
GT_SEP = re.compile(r"[/|]")

Contig = collections.namedtuple("Contig", "length offset width skip")


class Reference:
    """Random access to an uncompressed, faidx-indexed FASTA file"""

    def __init__(self, path):
        self.path = path
        self.index = collections.OrderedDict()
        with io.open(self.path + ".fai", "rb") as fp:
            for line in fp:
                flds = line.rstrip().decode().split("\t")
                if len(flds) != 5:
                    raise RuntimeError("Corrupted index")
                self.index[flds[0]] = Contig._make(map(int, flds[1:]))
        self.fp = io.open(path, "rb")

    def __getstate__(self):
        odict = self.__dict__.copy()
        odict["fp"] = self.fp.tell()
        return odict

    def __setstate__(self, ndict):
        fp = io.open(ndict["path"], "rb")
        fp.seek(ndict.pop("fp"))
        ndict["fp"] = fp
        self.__dict__.update(ndict)

    def get(self, c, s, e):
        """The 0-based, half-open reference sequence [s, e) of contig c"""
        ci = self.index[c]
        e = min(e, ci.length)
        if s > e:
            raise ValueError
        seq = b""
        while s < e:
            o = ci.offset + s // ci.width * ci.skip + s % ci.width
            n = min(ci.width - s % ci.width, e - s)
            self.fp.seek(o)
            seq += self.fp.read(n)
            s += n
        return seq.decode()


def load_regions(path):
    """Read a BED file into sorted, merged regions per contig

    Each contig maps to a pair of lists, the 0-based region starts and the
    exclusive region ends.
    """
    raw = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if not cols or cols[0] in ("track", "browser"):
                continue
            if cols[0].startswith("#"):
                continue
            raw[cols[0]].append((int(cols[1]), int(cols[2])))

    regions = {}
    for chrom, intervals in raw.items():
        starts, ends = [], []
        for s, e in sorted(intervals):
            if starts and s <= ends[-1]:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        regions[chrom] = (starts, ends)
    return regions


def in_regions(regions, chrom, pos):
    """Whether the 0-based position falls inside a region"""
    starts, ends = regions.get(chrom, ([], []))
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < ends[i]


def pieces(regions, chrom, start, end, include):
    """The parts of [start, end) inside, or outside, the regions"""
    starts, ends = regions.get(chrom, ([], []))
    out = []
    cur = start
    i = bisect.bisect_right(ends, start)  # first region ending after start
    while i < len(starts) and starts[i] < end:
        s, e = max(starts[i], start), min(ends[i], end)
        if include:
            out.append((s, e))
        elif cur < s:
            out.append((cur, s))
        cur = e
        i += 1
    if not include and cur < end:
        out.append((cur, end))
    return out


def set_end(info, end):
    """Set the INFO/END value of an INFO column"""
    flds = [] if info == "." else info.split(";")
    for i, fld in enumerate(flds):
        if fld.startswith("END="):
            flds[i] = "END=%d" % end
            break
    else:
        flds.append("END=%d" % end)
    return ";".join(flds)


def clip_block(v, start, end, reference):
    """A copy of the reference block `v` covering [start, end)"""
    if start == v.pos:
        ref = v.ref[0]
    else:
        ref = reference.get(v.chrom, start, start + 1).upper()
    cols = v.line.split("\t")
    cols[1] = str(start + 1)
    cols[3] = ref
    cols[7] = set_end(cols[7], end)
    info = dict(v.info, END=end)
    return vcflib.Variant(
        v.chrom,
        start,
        v.id,
        ref,
        v.alt,
        v.qual,
        v.filter,
        info,
        v.samples,
        end,
        "\t".join(cols),
    )


def fill_block(v, start, end, reference, alt, keys):
    """A reference block over [start, end) for the dropped record `v`

    The block claims no genotype confidence: GQ is 0 and the likelihoods
    are flat. DP and MIN_DP are the depth of `v`.
    """
    ref = reference.get(v.chrom, start, start + 1).upper()
    cols = [v.chrom, str(start + 1), ".", ref, alt, ".", ".", "END=%d" % end]
    if v.samples:
        cols.append(":".join(keys))
    for sample in v.samples:
        ploidy = len(GT_SEP.split(sample.get("GT") or "0/0"))
        dp = sample.get("DP")
        values = {
            "GT": "/".join(["0"] * ploidy),
            "DP": "." if dp is None else str(dp),
            "GQ": "0",
            "MIN_DP": "." if dp is None else str(dp),
            "PL": ",".join(["0"] * (ploidy + 1)),
        }
        cols.append(":".join(values[k] for k in keys))
    return vcflib.Variant(
        v.chrom,
        start,
        ".",
        ref,
        [alt],
        None,
        [],
        {"END": end},
        [{} for _ in v.samples],
        end,
        "\t".join(cols),
    )


class Writer:
    """Writes the records of one shard, or of the whole file, in order

    Kept records and block pieces are written by the shard they start in.
    A fill block waits for the records after it, which may cover part of
    it, and is cut at the shard boundaries, since the records that could
    cover it past a boundary belong to the next shard.
    """

    def __init__(self, outvcf, start, end):
        self.outvcf = outvcf
        self.start = start
        self.end = end
        self.chrom = None
        # The end of everything written or scheduled on the contig
        self.covered = 0
        # Scheduled fill blocks: (start, end, block builder)
        self.fills = []

    def record(self, v):
        """Write a kept record or a block piece"""
        self.contig(v.chrom)
        pending, self.fills = self.fills, []
        for start, end, make in pending:
            # The part before `v` is final; `v` covers the part it spans
            self.write_fill(start, min(end, v.pos), make)
            if max(start, v.end) < end:
                self.fills.append((max(start, v.end), end, make))
        self.covered = max(self.covered, v.end)
        if v.pos >= self.start:
            self.outvcf.emit(v)

    def fill(self, chrom, start, end, make):
        """Schedule a fill block over the uncovered part of [start, end)"""
        self.contig(chrom)
        start = max(start, self.covered)
        if start < end:
            self.fills.append((start, end, make))
            self.covered = end

    def contig(self, chrom):
        if chrom != self.chrom:
            self.finish()
            self.chrom, self.covered = chrom, 0

    def finish(self):
        """Write the scheduled fill blocks"""
        for start, end, make in self.fills:
            self.write_fill(start, end, make)
        self.fills = []

    def write_fill(self, start, end, make):
        start = max(start, self.start)
        if self.end is not None:
            end = min(end, self.end)
        if start < end:
            self.outvcf.emit(make(start, end))


def clip(invcf, outvcf, regions, reference, include, alt, keys):
    """Write the records, or parts of records, inside or outside regions

    When sharded, `invcf` yields every record overlapping the shard, so the
    records starting before the shard still count toward what is covered.
    """
    writer = Writer(
        outvcf, getattr(invcf, "start", -1), getattr(invcf, "end", None)
    )
    for v in invcf:
        if v.alt in REF_BLOCK_ALTS:
            for start, end in pieces(regions, v.chrom, v.pos, v.end, include):
                if start == v.pos and end == v.end:
                    writer.record(v)
                else:
                    writer.record(clip_block(v, start, end, reference))
        elif in_regions(regions, v.chrom, v.pos) == include:
            writer.record(v)
        else:
            make = functools.partial(
                fill_block, v, reference=reference, alt=alt, keys=keys
            )
            for start, end in pieces(regions, v.chrom, v.pos, v.end, include):
                writer.fill(v.chrom, start, end, make)
    writer.finish()


def main(args):
    regions = load_regions(args.bed)
    reference = Reference(args.reference)
    include = not args.exclude

    invcf = vcflib.VCF(args.vcf, "r")
    # The symbolic ALT of a fill block: the one the header declares
    alt = "<NON_REF>"
    if "*" in invcf.alts and "NON_REF" not in invcf.alts:
        alt = "<*>"
    keys = [k for k in FILL_KEYS if k in invcf.formats]
    outvcf = vcflib.VCF(args.output, "w")
    outvcf.copy_header(invcf)
    outvcf.emit_header()

    if args.threads < 2:
        clip(invcf, outvcf, regions, reference, include, alt, keys)
    else:
        sharder = vcflib.Sharder(args.threads)
        try:
            contig_lengths = [
                (contig, 0, int(d["length"]))
                for contig, d in invcf.contigs.items()
            ]
        except KeyError:
            return (
                "This script requires a VCF with contig lengths in"
                " the header when using multiple threads"
            )
        shards = sharder.cut(contig_lengths, args.step_size)
        sharder.run(
            shards,
            clip,
            [],
            invcf,
            outvcf,
            regions,
            reference,
            include,
            alt,
            keys,
        )

    outvcf.close()
    invcf.close()
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="gvcf_clip.py",
        usage="%(prog)s [options] -v GVCF -b BED -r REFERENCE output",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output", help="Output gVCF file name, or -")
    parser.add_argument("-v", "--vcf", required=True, help="Input gVCF")
    parser.add_argument("-b", "--bed", required=True, help="Region BED file")
    parser.add_argument(
        "-r", "--reference", required=True, help="Indexed reference FASTA"
    )
    parser.add_argument(
        "--exclude",
        action="store_true",
        help="Keep the positions outside the regions rather than inside",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=mp.cpu_count(),
        help="number of threads",
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=10 * 1000 * 1000,
        help=argparse.SUPPRESS,
    )
    sys.exit(main(parser.parse_args()))
