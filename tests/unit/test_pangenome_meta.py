"""
Unit tests for the pangenome metadata reader and its run-time checks.

The binary fixtures in `tests/data/pangenome` are built by the
`generate.sh` script next to them; see that script for the `vg` commands
and for what each graph contains.
"""

import pathlib
import struct
import subprocess
import sys

import pytest

from sentieon_cli.pangenome_meta import (
    HAPL_MAGIC,
    PangenomeMetadataError,
    check_gfa,
    check_sample_gbz,
    detect_pangenome_reference,
    main,
    read_gbz_metadata,
    read_hapl_header,
    unpack_bits,
)

DATA = pathlib.Path(__file__).resolve().parents[1] / "data" / "pangenome"

# The tiny graphs hold two reference samples over two 30 bp contigs
TINY_CONTIGS = ["chr1", "chr2"]


class TestUnpackBits:
    """The packed-integer decoder behind every `simple-sds` vector"""

    def test_byte_wide_values(self):
        words = struct.pack("<Q", 0x0807060504030201)
        assert unpack_bits(words, 8, 8) == [1, 2, 3, 4, 5, 6, 7, 8]

    def test_values_straddle_a_word_boundary(self):
        # Five 13-bit values need 65 bits, so the fifth spans two words
        values = [1, 4095, 8191, 0, 1234]
        packed = 0
        for i, value in enumerate(values):
            packed |= value << (13 * i)
        words = packed.to_bytes(16, "little")
        assert unpack_bits(words, 5, 13) == values

    def test_no_values(self):
        assert unpack_bits(b"", 0, 8) == []


class TestReadGbzMetadata:
    """Parsing the tags and path names out of a GBZ"""

    def test_samples_contigs_and_tags(self):
        meta = read_gbz_metadata(DATA / "tiny.gbz")
        assert meta.sample_names == ["GRCh38", "CHM13", "HG01", "HG02"]
        assert meta.contig_names == TINY_CONTIGS
        assert meta.reference_samples == ["GRCh38", "CHM13"]
        assert meta.reference_samples_tag == "GRCh38 CHM13"
        assert meta.gbwt_tags["source"] == "jltsiren/gbwt"

    def test_path_names(self):
        meta = read_gbz_metadata(DATA / "tiny.gbz")
        names = sorted(meta.full_path_name(p) for p in meta.path_names)
        assert names == [
            "CHM13#0#chr1",
            "CHM13#0#chr1",
            "CHM13#0#chr2",
            "GRCh38#0#chr1",
            "GRCh38#0#chr2",
            "HG01#1#chr1",
            "HG01#1#chr2",
            "HG01#2#chr1",
            "HG01#2#chr2",
            "HG02#1#chr1",
            "HG02#2#chr2",
        ]
        # The second CHM13 chr1 path is a fragment starting at offset 20
        fragments = sorted(
            path[3]
            for path in meta.path_names
            if meta.sample_of(path) == "CHM13"
        )
        assert fragments == [0, 0, 20]

    def test_only_the_head_of_the_file_is_read(self):
        """The BWT and the graph sequences are skipped, not read"""
        meta = read_gbz_metadata(DATA / "tiny.gbz")
        assert 0 < meta.bytes_read <= (DATA / "tiny.gbz").stat().st_size

    def test_a_gfa_is_not_a_gbz(self):
        with pytest.raises(PangenomeMetadataError, match="not a GBZ"):
            read_gbz_metadata(DATA / "tiny.gfa")

    def test_a_missing_file_raises_oserror(self):
        with pytest.raises(OSError):
            read_gbz_metadata(DATA / "does-not-exist.gbz")

    def test_a_truncated_file_is_rejected(self, tmp_path):
        truncated = tmp_path / "truncated.gbz"
        truncated.write_bytes((DATA / "tiny.gbz").read_bytes()[:64])
        with pytest.raises(PangenomeMetadataError):
            read_gbz_metadata(truncated)


class TestDetectPangenomeReference:
    """Picking the unfragmented backbone out of the reference samples"""

    def test_grch38_backbone(self):
        reference = detect_pangenome_reference(DATA / "tiny.gbz")
        assert reference.ref_name == "GRCh38"
        assert reference.contig_prefix == "GRCh38#0#"
        assert reference.contigs == TINY_CONTIGS
        assert reference.reference_samples == ["GRCh38", "CHM13"]
        assert reference.n_paths == 2

    def test_chm13_backbone(self):
        """The same graph with the two reference samples swapped"""
        reference = detect_pangenome_reference(DATA / "tiny_chm13.gbz")
        assert reference.ref_name == "CHM13"
        assert reference.contig_prefix == "CHM13#0#"
        assert reference.contigs == TINY_CONTIGS
        assert reference.reference_samples == ["CHM13", "GRCh38"]

    def test_sampled_graph_with_one_reference_sample(self):
        reference = detect_pangenome_reference(DATA / "tiny_sampled.gbz")
        assert reference.ref_name == "GRCh38"
        assert reference.reference_samples == ["GRCh38"]

    def test_not_a_gbz(self):
        with pytest.raises(PangenomeMetadataError, match="not a GBZ"):
            detect_pangenome_reference(DATA / "tiny.gfa")

    def test_no_reference_samples_tag(self):
        with pytest.raises(
            PangenomeMetadataError, match="'reference_samples' tag is missing"
        ):
            detect_pangenome_reference(DATA / "tiny_no_ref_tag.gbz")

    def test_no_unfragmented_candidate(self):
        with pytest.raises(
            PangenomeMetadataError, match="split into fragments"
        ) as excinfo:
            detect_pangenome_reference(DATA / "tiny_no_backbone.gbz")
        # The candidates are named, so the message says what was rejected
        assert "GRCh38" in str(excinfo.value)
        assert "CHM13" in str(excinfo.value)


class TestReadHaplHeader:
    """The `vg haplotypes` index header"""

    def test_tiny_hapl(self):
        header = read_hapl_header(DATA / "tiny.hapl")
        assert header.version == 6
        # One top-level chain per backbone contig
        assert header.top_level_chains == 2
        assert header.k == 29
        # Version 6 added the header tags
        assert "pggname" in header.tags

    def test_a_gbz_is_not_a_hapl(self):
        with pytest.raises(
            PangenomeMetadataError, match="not a vg haplotypes file"
        ):
            read_hapl_header(DATA / "tiny.gbz")

    def test_a_version_4_header_has_no_tags(self, tmp_path):
        """A hand-written v4 header, as in the HPRC v2.0 GRCh38 index"""
        path = tmp_path / "v4.hapl"
        path.write_bytes(
            struct.pack("<II", HAPL_MAGIC, 4)
            + struct.pack("<QQQQQ", 195, 16, 100, 12345, 29)
        )
        header = read_hapl_header(path)
        assert header.version == 4
        assert header.top_level_chains == 195
        assert header.construction_jobs == 16
        assert header.total_subchains == 100
        assert header.total_kmers == 12345
        assert header.k == 29
        assert header.tags == {}

    def test_a_truncated_header_is_rejected(self, tmp_path):
        path = tmp_path / "short.hapl"
        path.write_bytes(struct.pack("<II", HAPL_MAGIC, 4) + b"\x00" * 8)
        with pytest.raises(PangenomeMetadataError):
            read_hapl_header(path)


class TestCheckSampleGbz:
    """The run-time check on the graph `vg haplotypes` writes"""

    def test_a_sampled_graph_passes(self):
        check_sample_gbz(DATA / "tiny_sampled.gbz", "GRCh38", "GRCh38#0#")

    def test_more_than_one_reference_sample_fails(self):
        """`--set-reference` was never applied, so both samples remain"""
        with pytest.raises(PangenomeMetadataError) as excinfo:
            check_sample_gbz(DATA / "tiny.gbz", "GRCh38", "GRCh38#0#")
        message = str(excinfo.value)
        assert "reference_samples='GRCh38 CHM13'" in message
        assert "reference name may be wrong" in message

    def test_an_empty_reference_samples_tag_fails(self):
        """What `vg haplotypes --set-reference <wrong name>` writes"""
        with pytest.raises(PangenomeMetadataError) as excinfo:
            check_sample_gbz(
                DATA / "tiny_no_ref_tag.gbz", "GRCh38", "GRCh38#0#"
            )
        message = str(excinfo.value)
        assert "reference_samples=''" in message
        assert "reference name may be wrong" in message

    def test_a_reference_sample_without_the_prefix_fails(self):
        with pytest.raises(
            PangenomeMetadataError, match="contig prefix 'GRCh38#1#'"
        ):
            check_sample_gbz(
                DATA / "tiny_sampled.gbz", "GRCh38", "GRCh38#1#"
            )


class TestCheckGfa:
    """The run-time check on the GFA `vg convert -Q` writes"""

    def test_an_rgfa_with_the_reference_tags_passes(self):
        check_gfa(DATA / "rgfa_ok.gfa", "GRCh38", "GRCh38#0#")

    def test_a_header_that_still_lists_the_reference_fails(self):
        with pytest.raises(PangenomeMetadataError) as excinfo:
            check_gfa(
                DATA / "rgfa_wrong_reference.gfa", "GRCh38", "GRCh38#0#"
            )
        message = str(excinfo.value)
        assert "RS:Z:CHM13 GRCh38" in message
        assert "reference name may be wrong" in message

    def test_a_file_without_sn_tags_fails(self):
        with pytest.raises(PangenomeMetadataError, match="SN:Z:GRCh38#0#"):
            check_gfa(DATA / "rgfa_no_tags.gfa", "GRCh38", "GRCh38#0#")

    def test_a_file_without_a_header_line_is_tolerated(self):
        check_gfa(DATA / "rgfa_no_header.gfa", "GRCh38", "GRCh38#0#")

    def test_the_wrong_contig_prefix_fails(self):
        """The header passes, but no segment carries the prefix"""
        with pytest.raises(PangenomeMetadataError, match="SN:Z:GRCh38#1#"):
            check_gfa(DATA / "rgfa_ok.gfa", "GRCh38", "GRCh38#1#")

    def test_a_reference_name_in_the_header_fails_first(self):
        """`RS:Z:` naming the reference is checked before the tags"""
        with pytest.raises(PangenomeMetadataError, match="RS:Z:CHM13"):
            check_gfa(DATA / "rgfa_ok.gfa", "CHM13", "CHM13#0#")

    def test_only_the_head_of_the_file_is_read(self):
        """The tags are past `max_lines`, so the check gives up"""
        with pytest.raises(
            PangenomeMetadataError, match="first 1 lines"
        ):
            check_gfa(
                DATA / "rgfa_ok.gfa", "GRCh38", "GRCh38#0#", max_lines=1
            )


class TestCli:
    """The `python -m sentieon_cli.pangenome_meta` entry point"""

    def test_detect_prints_the_reference(self, capsys):
        assert main(["detect", "--gbz", str(DATA / "tiny.gbz")]) == 0
        out = capsys.readouterr().out
        assert "pangenome_ref_name\tGRCh38" in out
        assert "pangenome_contig_prefix\tGRCh38#0#" in out
        assert "n_contigs\t2" in out

    def test_detect_fails_on_a_graph_without_a_backbone(self):
        argv = ["detect", "--gbz", str(DATA / "tiny_no_backbone.gbz")]
        assert main(argv) == 1

    def test_check_gbz_passes(self):
        argv = [
            "check-gbz",
            "--gbz",
            str(DATA / "tiny_sampled.gbz"),
            "--reference_name",
            "GRCh38",
            "--contig_prefix",
            "GRCh38#0#",
        ]
        assert main(argv) == 0

    def test_check_gbz_fails(self):
        argv = [
            "check-gbz",
            "--gbz",
            str(DATA / "tiny.gbz"),
            "--reference_name",
            "GRCh38",
            "--contig_prefix",
            "GRCh38#0#",
        ]
        assert main(argv) == 1

    def test_check_gbz_fails_on_a_missing_file(self):
        argv = [
            "check-gbz",
            "--gbz",
            str(DATA / "does-not-exist.gbz"),
            "--reference_name",
            "GRCh38",
            "--contig_prefix",
            "GRCh38#0#",
        ]
        assert main(argv) == 1

    def test_check_gfa_passes(self):
        argv = [
            "check-gfa",
            "--gfa",
            str(DATA / "rgfa_ok.gfa"),
            "--reference_name",
            "GRCh38",
            "--contig_prefix",
            "GRCh38#0#",
        ]
        assert main(argv) == 0

    def test_check_gfa_fails(self):
        argv = [
            "check-gfa",
            "--gfa",
            str(DATA / "rgfa_wrong_reference.gfa"),
            "--reference_name",
            "GRCh38",
            "--contig_prefix",
            "GRCh38#0#",
        ]
        assert main(argv) == 1

    def test_check_gfa_honors_max_lines(self):
        argv = [
            "check-gfa",
            "--gfa",
            str(DATA / "rgfa_ok.gfa"),
            "--reference_name",
            "GRCh38",
            "--contig_prefix",
            "GRCh38#0#",
            "--max_lines",
            "1",
        ]
        assert main(argv) == 1

    def test_the_module_runs_as_a_script(self):
        """The DAG jobs run the module with `python -m`"""
        failing = subprocess.run(
            [
                sys.executable,
                "-m",
                "sentieon_cli.pangenome_meta",
                "check-gbz",
                "--gbz",
                str(DATA / "tiny.gbz"),
                "--reference_name",
                "GRCh38",
                "--contig_prefix",
                "GRCh38#0#",
            ],
            capture_output=True,
            text=True,
        )
        assert failing.returncode == 1
        assert "reference name may be wrong" in failing.stderr

        passing = subprocess.run(
            [
                sys.executable,
                "-m",
                "sentieon_cli.pangenome_meta",
                "check-gbz",
                "--gbz",
                str(DATA / "tiny_sampled.gbz"),
                "--reference_name",
                "GRCh38",
                "--contig_prefix",
                "GRCh38#0#",
            ],
            capture_output=True,
            text=True,
        )
        assert passing.returncode == 0
