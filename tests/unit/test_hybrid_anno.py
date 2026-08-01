import gzip
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

from sentieon_cli.command_strings import cmd_pyexec_hybrid_anno
from sentieon_cli.scripts.hybrid_anno import cut_shards


REPOSITORY_ROOT = pathlib.Path(__file__).parents[2]
HYBRID_ANNO = REPOSITORY_ROOT / "sentieon_cli" / "scripts" / "hybrid_anno.py"

INPUT_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=40>
##contig=<ID=chr2,length=20>
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=LHC,Number=1,Type=Integer,Description="old description">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT\t0/0
chr1\t6\t.\tC\t<NON_REF>\t.\tPASS\tEND=15\tGT\t0/0
chr1\t11\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
chr1\t20\t.\tA\tT\t40\tPASS\tDP=3\tGT\t0/1
chr1\t21\t.\tA\tT\t40\tPASS\tLHC=9\tGT\t0/1
chr1\t30\t.\tG\tC\t30\tPASS\t.\tGT\t1/1
chr1\t31\t.\tT\tC\t30\tPASS\t.\tGT\t1/1
chr2\t1\t.\tA\tG\t20\tPASS\t.\tGT\t0/1
chr2\t10\t.\tC\tT\t20\tPASS\t.\tGT\t0/1
"""

BED = """\
chr1\t0\t10\t1
chr1\t5\t20\t2
chr1\t20\t30\t3
chr2\t0\t5\t4
"""

EXPECTED_ANNOTATED_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=40>
##contig=<ID=chr2,length=20>
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=LHC,Number=1,Type=Integer,Description="Longread hap count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5;LHC=1\tGT\t0/0
chr1\t6\t.\tC\t<NON_REF>\t.\tPASS\tEND=15;LHC=1\tGT\t0/0
chr1\t11\t.\tA\tG\t50\tPASS\t.;LHC=2\tGT\t0/1
chr1\t20\t.\tA\tT\t40\tPASS\tDP=3;LHC=2\tGT\t0/1
chr1\t21\t.\tA\tT\t40\tPASS\tLHC=9;LHC=3\tGT\t0/1
chr1\t30\t.\tG\tC\t30\tPASS\t.;LHC=3\tGT\t1/1
chr1\t31\t.\tT\tC\t30\tPASS\t.\tGT\t1/1
chr2\t1\t.\tA\tG\t20\tPASS\t.;LHC=4\tGT\t0/1
chr2\t10\t.\tC\tT\t20\tPASS\t.\tGT\t0/1
"""

EXPECTED_EMPTY_BED_VCF = EXPECTED_ANNOTATED_VCF.replace(
    "END=5;LHC=1", "END=5"
).replace(
    "END=15;LHC=1", "END=15"
).replace(
    ".;LHC=2", "."
).replace(
    "DP=3;LHC=2", "DP=3"
).replace(
    "LHC=9;LHC=3", "LHC=9"
).replace(
    ".;LHC=3", "."
).replace(
    ".;LHC=4", "."
)


def require_htslib_tools() -> tuple[str, str]:
    bgzip = shutil.which("bgzip")
    tabix = shutil.which("tabix")
    assert bgzip is not None, "bgzip is required for hybrid_anno tests"
    assert tabix is not None, "tabix is required for hybrid_anno tests"
    return bgzip, tabix


def make_indexed_vcf(tmp_path: pathlib.Path) -> pathlib.Path:
    bgzip, tabix = require_htslib_tools()
    source = tmp_path / "input.vcf"
    source.write_text(INPUT_VCF)
    compressed = tmp_path / "input.vcf.gz"
    with compressed.open("wb") as output:
        subprocess.run([bgzip, "-c", str(source)], check=True, stdout=output)
    subprocess.run([tabix, "-p", "vcf", str(compressed)], check=True)
    return compressed


def run_annotator(
    tmp_path: pathlib.Path,
    input_vcf: pathlib.Path,
    bed_text: str,
    threads: int,
    step_size: int | None = None,
) -> tuple[pathlib.Path, subprocess.CompletedProcess[str]]:
    bed = tmp_path / f"stage1.{threads}.bed"
    bed.write_text(bed_text)
    output = tmp_path / f"annotated.{threads}.vcf.gz"
    scratch = tmp_path / f"scratch.{threads}"
    scratch.mkdir()
    command = [
        sys.executable,
        str(HYBRID_ANNO),
        "-v",
        str(input_vcf),
        "-b",
        str(bed),
        "-t",
        str(threads),
    ]
    if step_size is not None:
        command.extend(["--step-size", str(step_size)])
    command.append(str(output))
    environment = os.environ.copy()
    environment["SENTIEON_TMPDIR"] = str(scratch)
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        env=environment,
    )
    return output, result


@pytest.mark.parametrize(
    ("threads", "step_size"),
    ((1, None), (4, 13)),
)
def test_decompressed_output_matches_v170_oracle(
    tmp_path: pathlib.Path,
    threads: int,
    step_size: int | None,
) -> None:
    input_vcf = make_indexed_vcf(tmp_path)
    output, result = run_annotator(
        tmp_path, input_vcf, BED, threads, step_size
    )
    assert result.returncode == 0, result.stderr
    expected_step_size = step_size or 10_000_000
    assert f"step_size={expected_step_size}" in result.stderr
    assert pathlib.Path(f"{output}.tbi").is_file()
    with gzip.open(output, "rt") as annotated:
        assert annotated.read() == EXPECTED_ANNOTATED_VCF

    _, tabix = require_htslib_tools()
    observed = subprocess.check_output(
        [tabix, str(output), "chr1:6-21"],
        text=True,
    )
    expected = "\n".join(
        line
        for line in EXPECTED_ANNOTATED_VCF.splitlines()
        if line.startswith("chr1\t")
        and 6 <= int(line.split("\t", 2)[1]) <= 21
    )
    assert observed == f"{expected}\n"


def test_empty_bed_preserves_records_and_adds_header(
    tmp_path: pathlib.Path,
) -> None:
    input_vcf = make_indexed_vcf(tmp_path)
    output, result = run_annotator(tmp_path, input_vcf, "", 2, 10)
    assert result.returncode == 0, result.stderr
    with gzip.open(output, "rt") as annotated:
        assert annotated.read() == EXPECTED_EMPTY_BED_VCF


def test_shard_boundaries_carry_across_contigs_like_vcflib() -> None:
    assert cut_shards([("chr1", 15), ("chr2", 10)], 10) == [
        (0, "chr1", 0, 10),
        (1, "chr1", 10, 15),
        (2, "chr2", 0, 5),
        (3, "chr2", 5, 10),
    ]


def test_missing_index_fails_without_touching_input(
    tmp_path: pathlib.Path,
) -> None:
    input_vcf = make_indexed_vcf(tmp_path)
    pathlib.Path(f"{input_vcf}.tbi").unlink()
    output, result = run_annotator(tmp_path, input_vcf, BED, 1)
    assert result.returncode == 1
    assert "input VCF index is missing" in result.stderr
    assert input_vcf.is_file()
    assert not output.exists()


def test_invalid_thread_count_fails(tmp_path: pathlib.Path) -> None:
    input_vcf = make_indexed_vcf(tmp_path)
    output, result = run_annotator(tmp_path, input_vcf, BED, 0)
    assert result.returncode == 1
    assert "--threads must be at least 1" in result.stderr
    assert not output.exists()


def test_missing_contig_length_fails(tmp_path: pathlib.Path) -> None:
    bgzip, tabix = require_htslib_tools()
    source = tmp_path / "missing-length.vcf"
    source.write_text(
        INPUT_VCF.replace(
            "##contig=<ID=chr1,length=40>",
            "##contig=<ID=chr1>",
        )
    )
    input_vcf = tmp_path / "missing-length.vcf.gz"
    with input_vcf.open("wb") as output_file:
        subprocess.run([bgzip, "-c", str(source)], check=True, stdout=output_file)
    subprocess.run([tabix, "-p", "vcf", str(input_vcf)], check=True)
    output, result = run_annotator(tmp_path, input_vcf, BED, 1)
    assert result.returncode == 1
    assert "contig chr1 has no length" in result.stderr
    assert not output.exists()


@pytest.mark.parametrize(
    "bed_text",
    (
        "chr1\t0\t10\n",
        "chr1\tbad\t10\t1\n",
        "chr1\t10\t5\t1\n",
        "chr1\t10\t20\t1\nchr1\t0\t5\t2\n",
    ),
)
def test_malformed_bed_fails(
    tmp_path: pathlib.Path,
    bed_text: str,
) -> None:
    input_vcf = make_indexed_vcf(tmp_path)
    output, result = run_annotator(tmp_path, input_vcf, bed_text, 1)
    assert result.returncode == 1
    assert "Error:" in result.stderr
    assert not output.exists()


def test_command_builder_interface_is_unchanged(tmp_path: pathlib.Path) -> None:
    output = tmp_path / "output.vcf.gz"
    input_vcf = tmp_path / "input.vcf.gz"
    bed = tmp_path / "stage1_hap.bed"
    script = tmp_path / "hybrid_anno.py"
    pipeline = cmd_pyexec_hybrid_anno(output, input_vcf, bed, script, 128)
    assert len(pipeline.nodes) == 1
    command = pipeline.nodes[0]
    assert command.executable == sys.executable
    assert command.args == [
        str(script),
        "-v",
        str(input_vcf),
        "-b",
        str(bed),
        "-t",
        "128",
        str(output),
    ]
