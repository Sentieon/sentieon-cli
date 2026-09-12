"""
Unit tests for the pure-Python VCF header reader
"""

import gzip
import pathlib
import struct
import zlib

import pytest

from sentieon_cli.shard import parse_vcf_contigs, vcf_contigs
from sentieon_cli.util import (
    VcfHeaderError,
    read_vcf_header,
    vcf_id,
    vcf_id_from_header,
)

HEADER = """##fileformat=VCFv4.2
##SentieonVcfID=population-test-20260101
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chrUn_KI270302v1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

BODY = "chr1\t1\t.\tA\tC\t50\tPASS\t.\n#not-a-header-line\n"


def bgzf(data: bytes) -> bytes:
    """A single BGZF block, the way a `bgzip`-compressed VCF starts.

    BGZF is gzip with an extra `BC` subfield, so the magic bytes the
    reader sniffs are the same and `gzip` decompresses it unchanged.
    """
    deflater = zlib.compressobj(6, zlib.DEFLATED, -zlib.MAX_WBITS)
    deflated = deflater.compress(data) + deflater.flush()
    bsize = 18 + len(deflated) + 8 - 1
    header = struct.pack(
        "<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize
    )
    trailer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + deflated + trailer


@pytest.fixture
def plain_vcf(tmp_path: pathlib.Path) -> pathlib.Path:
    """A plain-text VCF, named `.vcf.gz` to defeat the extension"""
    path = tmp_path / "plain.vcf.gz"
    path.write_text(HEADER + BODY)
    return path


@pytest.fixture
def gzip_vcf(tmp_path: pathlib.Path) -> pathlib.Path:
    """A gzip-compressed VCF, named `.vcf` to defeat the extension"""
    path = tmp_path / "compressed.vcf"
    with gzip.open(path, "wt") as fh:
        fh.write(HEADER + BODY)
    return path


@pytest.fixture
def bgzf_vcf(tmp_path: pathlib.Path) -> pathlib.Path:
    """A BGZF-compressed VCF, as `bgzip` would write it"""
    path = tmp_path / "bgzf.vcf.gz"
    path.write_bytes(bgzf((HEADER + BODY).encode()))
    return path


class TestReadVcfHeader:
    """`read_vcf_header` decompresses by magic bytes, not by name"""

    def test_a_plain_text_vcf_is_read(self, plain_vcf):
        assert read_vcf_header(plain_vcf) == HEADER.splitlines()

    def test_a_gzip_vcf_is_read(self, gzip_vcf):
        assert read_vcf_header(gzip_vcf) == HEADER.splitlines()

    def test_a_bgzf_vcf_is_read(self, bgzf_vcf):
        assert read_vcf_header(bgzf_vcf) == HEADER.splitlines()

    def test_the_body_is_never_read(self, plain_vcf):
        """Reading stops at `#CHROM`, so a `#` line below it is not a
        header line"""
        header = read_vcf_header(plain_vcf)
        assert header[-1].startswith("#CHROM")
        assert "#not-a-header-line" not in header

    def test_a_header_without_chrom_stops_at_the_body(self, tmp_path):
        path = tmp_path / "no_chrom.vcf"
        path.write_text("##fileformat=VCFv4.2\nchr1\t1\t.\tA\tC\n")
        assert read_vcf_header(path) == ["##fileformat=VCFv4.2"]

    def test_a_file_with_no_header_raises(self, tmp_path):
        path = tmp_path / "headerless.vcf"
        path.write_text("chr1\t1\t.\tA\tC\t50\tPASS\t.\n")
        with pytest.raises(VcfHeaderError):
            read_vcf_header(path)

    def test_an_empty_file_raises(self, tmp_path):
        path = tmp_path / "empty.vcf.gz"
        path.touch()
        with pytest.raises(VcfHeaderError):
            read_vcf_header(path)

    def test_a_missing_file_raises_oserror(self, tmp_path):
        with pytest.raises(OSError):
            read_vcf_header(tmp_path / "absent.vcf.gz")

    def test_truncated_gzip_raises(self, tmp_path, gzip_vcf):
        path = tmp_path / "truncated.vcf.gz"
        path.write_bytes(gzip_vcf.read_bytes()[:20])
        with pytest.raises((VcfHeaderError, OSError)):
            read_vcf_header(path)


class TestVcfContigs:
    """`vcf_contigs` keeps its return type and its `None` lengths"""

    def test_contigs_and_lengths_are_parsed(self, bgzf_vcf):
        assert vcf_contigs(bgzf_vcf) == {
            "chr1": 248956422,
            "chr2": 242193529,
            "chrUn_KI270302v1": None,
        }

    def test_a_contig_without_a_length_is_none(self):
        contigs = parse_vcf_contigs(["##contig=<ID=chr1>"])
        assert contigs == {"chr1": None}

    def test_a_quoted_assembly_does_not_confuse_the_parser(self):
        contigs = parse_vcf_contigs(
            ['##contig=<ID=chr1,length=10,assembly="b,37">']
        )
        assert contigs == {"chr1": 10}

    def test_a_missing_file_returns_an_empty_dict(self, tmp_path):
        assert vcf_contigs(tmp_path / "absent.vcf.gz") == {}

    def test_a_headerless_file_returns_an_empty_dict(self, tmp_path):
        path = tmp_path / "headerless.vcf"
        path.write_text("chr1\t1\t.\tA\tC\n")
        assert vcf_contigs(path) == {}


class TestVcfId:
    """`vcf_id` reports `None` the same way `vcf_contigs` reports `{}`"""

    def test_the_id_is_read(self, gzip_vcf):
        assert vcf_id(gzip_vcf) == "population-test-20260101"

    def test_a_header_without_an_id(self, tmp_path):
        path = tmp_path / "no_id.vcf"
        path.write_text("##fileformat=VCFv4.2\n")
        assert vcf_id(path) is None

    def test_a_missing_file_returns_none(self, tmp_path):
        assert vcf_id(tmp_path / "absent.vcf.gz") is None

    def test_an_empty_file_returns_none(self, tmp_path):
        path = tmp_path / "empty.vcf.gz"
        path.touch()
        assert vcf_id(path) is None

    def test_vcf_id_from_header(self):
        assert vcf_id_from_header(HEADER.splitlines()) == (
            "population-test-20260101"
        )
        assert vcf_id_from_header([]) is None


if __name__ == "__main__":
    pytest.main([__file__])
