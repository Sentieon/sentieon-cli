#!/usr/bin/env python
"""Restrict a gVCF to the inside, or the outside, of a set of BED regions

The DNAscope hybrid pipeline splices its two calling passes together: the
first pass covers the positions outside the realigned regions and the
second pass the positions inside them. Variant records are kept or dropped
by their start position, as with `bcftools view -T`. Reference blocks are
clipped at the region boundaries, so the spliced gVCF has neither
overlapping blocks nor gaps.
"""

import argparse
import bisect
import collections
import io
import multiprocessing as mp
import sys

import vcflib

REF_BLOCK_ALTS = (["<NON_REF>"], ["<*>"])

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


def clip(invcf, outvcf, regions, reference, include):
    """Write the records, or parts of records, inside or outside regions

    When sharded, `invcf` yields every record overlapping the shard, so a
    record or a block piece is written only by the shard it starts in.
    """
    shard_start = getattr(invcf, "start", -1)
    for v in invcf:
        if v.alt not in REF_BLOCK_ALTS:
            if v.pos < shard_start:
                continue
            if in_regions(regions, v.chrom, v.pos) == include:
                outvcf.emit(v)
            continue

        for start, end in pieces(regions, v.chrom, v.pos, v.end, include):
            if start < shard_start:
                continue
            if start == v.pos and end == v.end:
                outvcf.emit(v)
            else:
                outvcf.emit(clip_block(v, start, end, reference))


def main(args):
    regions = load_regions(args.bed)
    reference = Reference(args.reference)
    include = not args.exclude

    invcf = vcflib.VCF(args.vcf, "r")
    outvcf = vcflib.VCF(args.output, "w")
    outvcf.copy_header(invcf)
    outvcf.emit_header()

    if args.threads < 2:
        clip(invcf, outvcf, regions, reference, include)
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
            shards, clip, [], invcf, outvcf, regions, reference, include
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
    parser.add_argument("output", help="Output gVCF file name")
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
