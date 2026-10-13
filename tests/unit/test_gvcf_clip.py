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
