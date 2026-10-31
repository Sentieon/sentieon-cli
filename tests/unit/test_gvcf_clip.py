"""
Unit tests for the `gvcf_clip.py` helper script.

The script is not importable as part of the package (the `scripts`
directory has no `__init__.py`), so its helpers are loaded from the file
path, and the end-to-end tests run it the way the pipeline does, with
`sys.executable`.
"""

import gzip
import importlib.util
import pathlib
import subprocess
import sys

import pytest
import vcflib

SCRIPT = (
    pathlib.Path(__file__).resolve().parents[2]
    / "sentieon_cli"
    / "scripts"
    / "gvcf_clip.py"
)


def _load_script():
    spec = importlib.util.spec_from_file_location("gvcf_clip", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gvcf_clip = _load_script()

SEQ = "ACGT" * 25
HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100>",
    '##ALT=<ID=NON_REF,Description="Any other allele">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">',
    '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Min depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]
BLOCK_SAMPLE = "GT:PL\t0/0:0,30,300"
VARIANT_SAMPLE = "GT:PL\t0/1:500,0,600,500,600,900"
# Records in 1-based coordinates, around the BED regions 21-30 and 61-70
RECORDS = [
    f"chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=25;DP=10\t{BLOCK_SAMPLE}",
    f"chr1\t26\t.\tC\t<NON_REF>\t.\t.\t.\t{BLOCK_SAMPLE}",
    f"chr1\t27\t.\tG\tC,<NON_REF>\t50\t.\t.\t{VARIANT_SAMPLE}",
    f"chr1\t28\t.\tT\t<NON_REF>\t.\t.\tEND=40\t{BLOCK_SAMPLE}",
    f"chr1\t41\t.\tA\tT,<NON_REF>\t50\t.\t.\t{VARIANT_SAMPLE}",
    f"chr1\t42\t.\tC\t<NON_REF>\t.\t.\tEND=100\t{BLOCK_SAMPLE}",
]
REGIONS = "chr1\t20\t30\nchr1\t60\t70\n"


def _write_inputs(tmp_path, records=RECORDS, regions=REGIONS):
    """Write the reference, the region BED, and a bgzipped, indexed gVCF"""
    ref = tmp_path / "ref.fa"
    # Two FASTA lines, so lookups cross a line break
    ref.write_text(f">chr1\n{SEQ[:60]}\n{SEQ[60:]}\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t100\t6\t60\t61\n")

    bed = tmp_path / "regions.bed"
    bed.write_text(regions)

    gvcf = tmp_path / "in.g.vcf.gz"
    out = vcflib.VCF(str(gvcf), "w")
    out.headers = list(HEADER)
    out.parse_header()
    out.emit_header()
    for line in records:
        out.emit(out.parse(line))
    out.close()
    return gvcf, bed, ref


def _run(tmp_path, exclude, threads=1, step_size=None, **inputs):
    gvcf, bed, ref = _write_inputs(tmp_path, **inputs)
    name = "exclude" if exclude else "include"
    output = tmp_path / f"{name}-{threads}.g.vcf.gz"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "-v",
        str(gvcf),
        "-b",
        str(bed),
        "-r",
        str(ref),
        "-t",
        str(threads),
        str(output),
    ]
    if exclude:
        cmd.append("--exclude")
    if step_size:
        cmd.extend(["--step-size", str(step_size)])
    subprocess.run(cmd, check=True)
    assert pathlib.Path(str(output) + ".tbi").exists()
    with gzip.open(output, "rt") as fh:
        return [line.rstrip("\n") for line in fh if not line.startswith("#")]


def _columns(lines):
    """The CHROM, POS, REF, ALT and INFO columns of each record"""
    return [tuple(line.split("\t")[i] for i in (0, 1, 3, 4, 7)) for line in lines]


def _covered(lines):
    """Every 1-based position the records span, with repeats"""
    positions = []
    for line in lines:
        cols = line.split("\t")
        start = int(cols[1])
        end = start + len(cols[3]) - 1
        for fld in cols[7].split(";"):
            if fld.startswith("END="):
                end = int(fld[4:])
        positions.extend(range(start, end + 1))
    return positions


class TestRegions:
    def test_load_regions_sorts_and_merges(self, tmp_path):
        bed = tmp_path / "regions.bed"
        bed.write_text(
            "track name=test\n"
            "# comment\n"
            "chr1\t50\t60\textra\n"
            "chr1\t10\t20\n"
            "chr1\t15\t30\n"
            "chr1\t30\t35\n"
            "chr2\t0\t5\n"
        )
        regions = gvcf_clip.load_regions(str(bed))
        assert regions == {
            "chr1": ([10, 50], [35, 60]),
            "chr2": ([0], [5]),
        }

    def test_in_regions(self):
        regions = {"chr1": ([10, 50], [20, 60])}
        assert not gvcf_clip.in_regions(regions, "chr1", 9)
        assert gvcf_clip.in_regions(regions, "chr1", 10)
        assert gvcf_clip.in_regions(regions, "chr1", 19)
        assert not gvcf_clip.in_regions(regions, "chr1", 20)
        assert gvcf_clip.in_regions(regions, "chr1", 55)
        assert not gvcf_clip.in_regions(regions, "chr2", 15)

    @pytest.mark.parametrize(
        "start, end, include, expected",
        [
            # Outside every region
            (0, 5, True, []),
            (0, 5, False, [(0, 5)]),
            # Across the start of a region
            (5, 15, True, [(10, 15)]),
            (5, 15, False, [(5, 10)]),
            # Across the end of a region
            (15, 25, True, [(15, 20)]),
            (15, 25, False, [(20, 25)]),
            # Within a region
            (12, 18, True, [(12, 18)]),
            (12, 18, False, []),
            # Spanning two regions
            (0, 70, True, [(10, 20), (50, 60)]),
            (0, 70, False, [(0, 10), (20, 50), (60, 70)]),
            # Ending at a region start, starting at a region end
            (0, 10, True, []),
            (20, 50, False, [(20, 50)]),
        ],
    )
    def test_pieces(self, start, end, include, expected):
        regions = {"chr1": ([10, 50], [20, 60])}
        assert (
            gvcf_clip.pieces(regions, "chr1", start, end, include) == expected
        )

    def test_pieces_without_regions_on_contig(self):
        assert gvcf_clip.pieces({}, "chr1", 0, 10, False) == [(0, 10)]
        assert gvcf_clip.pieces({}, "chr1", 0, 10, True) == []

    @pytest.mark.parametrize(
        "info, expected",
        [
            ("END=40", "END=50"),
            ("DP=10;END=40;MQ=60", "DP=10;END=50;MQ=60"),
            ("BLOCKEND=3", "BLOCKEND=3;END=50"),
            (".", "END=50"),
        ],
    )
    def test_set_end(self, info, expected):
        assert gvcf_clip.set_end(info, 50) == expected


class TestClip:
    def test_exclude_clips_blocks_at_region_boundaries(self, tmp_path):
        records = _run(tmp_path, exclude=True)
        assert _columns(records) == [
            ("chr1", "1", "A", "<NON_REF>", "END=20;DP=10"),
            ("chr1", "31", "G", "<NON_REF>", "END=40"),
            ("chr1", "41", "A", "T,<NON_REF>", "."),
            ("chr1", "42", "C", "<NON_REF>", "END=60"),
            ("chr1", "71", "G", "<NON_REF>", "END=100"),
        ]

    def test_include_keeps_the_region_parts(self, tmp_path):
        records = _run(tmp_path, exclude=False)
        assert _columns(records) == [
            ("chr1", "21", "A", "<NON_REF>", "END=25;DP=10"),
            ("chr1", "26", "C", "<NON_REF>", "."),
            ("chr1", "27", "G", "C,<NON_REF>", "."),
            ("chr1", "28", "T", "<NON_REF>", "END=30"),
            ("chr1", "61", "A", "<NON_REF>", "END=70"),
        ]

    def test_clipped_blocks_keep_the_sample_columns(self, tmp_path):
        records = _run(tmp_path, exclude=True)
        for line in records:
            cols = line.split("\t")
            if cols[4] == "<NON_REF>":
                assert "\t".join(cols[8:]) == BLOCK_SAMPLE

    def test_the_two_sides_tile_the_input(self, tmp_path):
        outside = _run(tmp_path, exclude=True)
        inside = _run(tmp_path, exclude=False)
        positions = _covered(outside) + _covered(inside)
        assert sorted(positions) == list(range(1, 101))

    @pytest.mark.parametrize("exclude", [True, False])
    def test_sharded_matches_single_thread(self, tmp_path, exclude):
        single = _run(tmp_path, exclude)
        # Shard boundaries at 30 and 60 fall inside blocks
        sharded = _run(tmp_path, exclude, threads=2, step_size=30)
        assert sharded == single

    def test_empty_bed(self, tmp_path):
        assert _run(tmp_path, exclude=True, regions="") == RECORDS
        assert _run(tmp_path, exclude=False, regions="") == []


def _block(pos, end):
    return (
        f"chr1\t{pos}\t.\t{SEQ[pos - 1]}\t<NON_REF>\t.\t.\tEND={end}\t"
        "GT:DP:GQ:MIN_DP:PL\t0/0:10:30:8:0,30,300"
    )


def _deletion(pos, length, gt="0/1", dp=12):
    ref = SEQ[pos - 1 : pos - 1 + length]
    return (
        f"chr1\t{pos}\t.\t{ref}\t{ref[0]},<NON_REF>\t50\t.\t.\t"
        f"GT:DP:PL\t{gt}:{dp}:500,0,600,500,600,900"
    )


def _snv(pos):
    ref = SEQ[pos - 1]
    alt = "C" if ref != "C" else "T"
    return f"chr1\t{pos}\t.\t{ref}\t{alt},<NON_REF>\t50\t.\t.\t{VARIANT_SAMPLE}"


def _fill(pos, end, dp=12):
    return (
        "chr1",
        str(pos),
        SEQ[pos - 1],
        "<NON_REF>",
        f"END={end}",
        "GT:DP:GQ:MIN_DP:PL",
        f"0/0:{dp}:0:{dp}:0,0,0",
    )


def _summary(lines):
    """POS, END or REF span, and whether each record is a fill block"""
    out = []
    for line in lines:
        cols = line.split("\t")
        out.append(tuple(cols[i] for i in (0, 1, 3, 4, 7)) + tuple(cols[8:]))
    return out


OUTSIDE = list(range(1, 21)) + list(range(31, 61)) + list(range(71, 101))
INSIDE = list(range(21, 31)) + list(range(61, 71))


class TestFill:
    """Positions spanned only by a dropped variant record are filled"""

    def test_tail_of_a_dropped_deletion(self, tmp_path):
        # The deletion starts at 29, inside 21-30, and spans 29-33
        records = [_block(1, 28), _deletion(29, 5), _block(34, 60)]
        records.append(_block(61, 100))
        out = _run(tmp_path, exclude=True, records=records)

        summary = _summary(out)
        assert _fill(31, 33) in summary
        assert [x[1] for x in summary] == ["1", "31", "34", "71"]
        assert sorted(_covered(out)) == OUTSIDE

    def test_records_inside_the_tail_split_the_fill(self, tmp_path):
        # The deletion spans 29-36; a SNV at 33 is inside its span
        records = [_block(1, 28), _deletion(29, 8), _snv(33)]
        records += [_block(37, 60), _block(61, 100)]
        out = _run(tmp_path, exclude=True, records=records)

        summary = _summary(out)
        assert [x[1] for x in summary] == ["1", "31", "33", "34", "37", "71"]
        assert summary[1] == _fill(31, 32)
        assert summary[3] == _fill(34, 36)
        assert sorted(_covered(out)) == OUTSIDE

    def test_a_block_inside_the_tail_trims_the_fill(self, tmp_path):
        records = [_block(1, 28), _deletion(29, 8), _block(33, 60)]
        records.append(_block(61, 100))
        out = _run(tmp_path, exclude=True, records=records)

        summary = _summary(out)
        assert [x[1] for x in summary] == ["1", "31", "33", "71"]
        assert summary[1] == _fill(31, 32)
        assert sorted(_covered(out)) == OUTSIDE

    def test_no_fill_where_a_kept_record_covers_the_tail(self, tmp_path):
        # The deletion at 19 starts outside and spans 19-34
        records = [_block(1, 18), _deletion(19, 16), _deletion(29, 5)]
        records += [_block(35, 60), _block(61, 100)]
        out = _run(tmp_path, exclude=True, records=records)

        assert [x[1] for x in _summary(out)] == ["1", "19", "35", "71"]

    def test_include_side(self, tmp_path):
        # The deletion at 59 starts outside 61-70 and spans 59-63
        records = [_block(1, 58), _deletion(59, 5, gt="0/0", dp=7)]
        records.append(_block(64, 100))
        out = _run(tmp_path, exclude=False, records=records)

        summary = _summary(out)
        assert [x[1] for x in summary] == ["21", "61", "64"]
        assert summary[1] == _fill(61, 63, dp=7)
        assert sorted(_covered(out)) == INSIDE

    def test_the_two_sides_tile_the_contig(self, tmp_path):
        records = [_block(1, 28), _deletion(29, 8), _snv(33)]
        records += [_block(37, 58), _deletion(59, 5), _block(64, 100)]
        outside = _run(tmp_path, exclude=True, records=records)
        inside = _run(tmp_path, exclude=False, records=records)
        positions = _covered(outside) + _covered(inside)
        assert set(positions) == set(range(1, 101))
        # The deletion kept at 59 spans 61-63, inside 61-70, as in the input
        assert sorted(_covered(outside)) == sorted(OUTSIDE + [61, 62, 63])

    def test_sharded_matches_single_thread(self, tmp_path):
        records = [_block(1, 28), _deletion(29, 8), _snv(33)]
        records += [_block(37, 60), _block(61, 100)]
        single = _run(tmp_path, exclude=True, records=records)
        # A shard boundary between positions 32 and 33
        sharded = _run(
            tmp_path, exclude=True, threads=2, step_size=32, records=records
        )
        assert sharded == single

    def test_fill_is_cut_at_a_shard_boundary(self, tmp_path):
        records = [_block(1, 28), _deletion(29, 5), _block(34, 60)]
        records.append(_block(61, 100))
        out = _run(
            tmp_path, exclude=True, threads=2, step_size=32, records=records
        )

        summary = _summary(out)
        assert summary[1] == _fill(31, 32)
        assert summary[2] == _fill(33, 33)
        assert sorted(_covered(out)) == OUTSIDE

    def test_stdout_output(self, tmp_path):
        records = [_block(1, 28), _deletion(29, 5), _block(34, 60)]
        records.append(_block(61, 100))
        expected = _run(tmp_path, exclude=True, records=records)
        gvcf, bed, ref = _write_inputs(tmp_path, records=records)
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "-v",
                str(gvcf),
                "-b",
                str(bed),
                "-r",
                str(ref),
                "-t",
                "2",
                "--exclude",
                "-",
            ],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        lines = result.stdout.splitlines()
        assert lines[0] == "##fileformat=VCFv4.2"
        assert [x for x in lines if not x.startswith("#")] == expected
