import gzip
import pathlib
import shutil
import subprocess
import sys

import pytest

from sentieon_cli.command_strings import cmd_pyexec_hybrid_transfer
from sentieon_cli.scripts.hybrid_transfer import (
    build_work_items,
    trim_record,
)

REPOSITORY_ROOT = pathlib.Path(__file__).parents[2]
HYBRID_TRANSFER = (
    REPOSITORY_ROOT / "sentieon_cli" / "scripts" / "hybrid_transfer.py"
)
TRIMALT = REPOSITORY_ROOT / "sentieon_cli" / "scripts" / "trimalt.py"

RAW_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=25>
##contig=<ID=chrUn,length=7>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=AF,Number=A,Type=Float,Description="Raw allele frequency">
##INFO=<ID=IA,Number=A,Type=Integer,Description="Number A fixture">
##INFO=<ID=IR,Number=R,Type=Integer,Description="Number R fixture">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=FA,Number=A,Type=Integer,Description="Number A fixture">
##FORMAT=<ID=FR,Number=R,Type=Integer,Description="Number R fixture">
##FORMAT=<ID=FG,Number=G,Type=Integer,Description="Number G fixture">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t1\t.\tA\tC\t50\tPASS\tAF=0.5;IA=11;IR=1,2\tGT:FA:FR:FG\t0/1:11:1,2:0,10,20
chr1\t2\t.\tA\tC,G\t50\tPASS\tAF=0.1,0.2;IA=12,13;IR=1,2,3\tGT:FA:FR:FG\t1/2:12,13:1,2,3:0,10,20,30,40,50
chr1\t3\t.\tG\t<NON_REF>\t.\tPASS\tEND=5\tGT:FA:FR:FG\t0/0:.:1,.:0,10,20
chr1\t9\t.\tT\tC\t40\tPASS\tAF=0.3;IA=19;IR=3,4\tGT:FA:FR:FG\t0/1:19:3,4:0,11,22
chr1\t10\t.\tA\tG\t40\tPASS\tAF=0.4;IA=20;IR=4,5\tGT:FA:FR:FG\t0/1:20:4,5:0,12,24
chr1\t11\t.\tA\tT\t40\tPASS\tAF=0.5;IA=21;IR=5,6\tGT:FA:FR:FG\t0/1:21:5,6:0,13,26
chr1\t12\ta\tC\tG\t40\tPASS\tAF=0.6;IA=22;IR=6,7\tGT:FA:FR:FG\t0/1:22:6,7:0,14,28
chr1\t12\tb\tC\tT\t40\tPASS\tAF=0.7;IA=23;IR=7,8\tGT:FA:FR:FG\t1:23:7,8:0,15
chr1\t15\t.\tAT\tCT\t40\tPASS\tAF=0.8;IA=24;IR=8,9\tGT:FA:FR:FG\t0/1:24:8,9:0,16,32
chr1\t20\t.\tG\tA\t40\tPASS\tAF=.;IA=25;IR=9,10\tGT:FA:FR:FG\t./.:25:9,10:0,17,34
chrUn\t1\t.\tA\tT\t30\tPASS\tAF=0.9;IA=31;IR=10,11\tGT:FA:FR:FG\t0/1:31:10,11:0,18,36
chrUn\t7\t.\tC\tG\t30\tPASS\tAF=0.4;IA=32;IR=11,12\tGT:FA:FR:FG\t0/1:32:11,12:0,19,38
"""

POPULATION_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=25>
##INFO=<ID=POPA,Number=A,Type=Integer,Description="Population count">
##INFO=<ID=POPR,Number=R,Type=Integer,Description="Population ref count">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t2\t.\tA\tC,T\t.\tPASS\tPOPA=5,7;POPR=20,5,7
chr1\t3\t.\tG\tA\t.\tPASS\tPOPA=8;POPR=20,8
chr1\t4\t.\tC\tT\t.\tPASS\tPOPA=9;POPR=20,9
chr1\t10\t.\tA\tG,C\t.\tPASS\tPOPA=10,11;POPR=20,10,11
chr1\t11\t.\tA\tC\t.\tPASS\tPOPA=12;POPR=20,12
chr1\t12\t.\tC\tA\t.\tPASS\tPOPA=13;POPR=20,13
chr1\t15\t.\tAT\tGT\t.\tPASS\tPOPA=14;POPR=20,14
chr1\t21\t.\tT\tA\t.\tPASS\tPOPA=15;POPR=20,15
"""


def require_tools() -> tuple[str, str, str]:
    bcftools = shutil.which("bcftools")
    bgzip = shutil.which("bgzip")
    tabix = shutil.which("tabix")
    assert bcftools is not None, "bcftools is required"
    assert bgzip is not None, "bgzip is required"
    assert tabix is not None, "tabix is required"
    return bcftools, bgzip, tabix


def make_indexed_vcf(
    tmp_path: pathlib.Path, name: str, contents: str
) -> pathlib.Path:
    _, bgzip, tabix = require_tools()
    source = tmp_path / f"{name}.vcf"
    source.write_text(contents)
    compressed = tmp_path / f"{name}.vcf.gz"
    with compressed.open("wb") as output:
        subprocess.run([bgzip, "-c", str(source)], check=True, stdout=output)
    subprocess.run([tabix, "-p", "vcf", str(compressed)], check=True)
    return compressed


def make_inputs(
    tmp_path: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path, pathlib.Path]:
    raw = make_indexed_vcf(tmp_path, "raw", RAW_VCF)
    population = make_indexed_vcf(tmp_path, "population", POPULATION_VCF)
    reference_index = tmp_path / "reference.fa.fai"
    reference_index.write_text("chr1\t25\t0\t25\t26\nchrUn\t7\t26\t7\t8\n")
    return raw, population, reference_index


def run_pipeline(
    commands: list[list[str]],
) -> tuple[bytes, bytes]:
    stdin: bytes | None = None
    stderr = bytearray()
    for command in commands:
        result = subprocess.run(
            command,
            input=stdin,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stderr.extend(result.stderr)
        assert result.returncode == 0, stderr.decode(errors="replace")
        stdin = result.stdout
    return stdin or b"", bytes(stderr)


def run_legacy_transfer(
    tmp_path: pathlib.Path,
    raw: pathlib.Path,
    population: pathlib.Path,
    reference_index: pathlib.Path,
    step_size: int,
) -> pathlib.Path:
    bcftools, _, _ = require_tools()
    population_header = subprocess.check_output(
        [bcftools, "view", "-h", str(population)], text=True
    )
    info_fields = []
    for line in population_header.splitlines():
        if line.startswith("##INFO") and ",Number=A" in line:
            info_fields.append(line.split("ID=", 1)[1].split(",", 1)[0])
    merge_rules = ",".join(f"{field}:sum" for field in info_fields)
    items = build_work_items([("chr1", 25), ("chrUn", 7)], {"chr1"}, step_size)
    shards: list[pathlib.Path] = []
    for item in items:
        bed = tmp_path / f"legacy.{item.number}.bed"
        bed.write_text(f"{item.contig}\t{item.start}\t{item.stop}\n")
        shard = tmp_path / f"legacy.{item.number}.vcf.gz"
        if item.merge_population:
            merged, _ = run_pipeline(
                [
                    [
                        bcftools,
                        "merge",
                        "--regions-file",
                        str(bed),
                        "--no-version",
                        "--regions-overlap",
                        "pos",
                        "-m",
                        "all",
                        "-i",
                        merge_rules,
                        str(raw),
                        str(population),
                    ],
                    [sys.executable, str(TRIMALT)],
                ]
            )
            result = subprocess.run(
                [
                    bcftools,
                    "view",
                    "--no-version",
                    "-W=tbi",
                    "-o",
                    str(shard),
                ],
                input=merged,
                capture_output=True,
            )
        else:
            result = subprocess.run(
                [
                    bcftools,
                    "view",
                    "--no-version",
                    "-W=tbi",
                    "-O",
                    "z",
                    "-o",
                    str(shard),
                    "--regions-file",
                    str(bed),
                    str(raw),
                ],
                capture_output=True,
            )
        assert result.returncode == 0, result.stderr.decode(errors="replace")
        shards.append(shard)

    output = tmp_path / "legacy.vcf.gz"
    result = subprocess.run(
        [
            bcftools,
            "concat",
            "-W=tbi",
            "--output",
            str(output),
            "--no-version",
            "--threads",
            "4",
            *(str(shard) for shard in shards),
        ],
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr.decode(errors="replace")
    return output


def run_fused_transfer(
    tmp_path: pathlib.Path,
    raw: pathlib.Path,
    population: pathlib.Path,
    reference_index: pathlib.Path,
) -> tuple[pathlib.Path, subprocess.CompletedProcess[str]]:
    output = tmp_path / "fused.vcf.gz"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = subprocess.run(
        [
            sys.executable,
            str(HYBRID_TRANSFER),
            "--raw-vcf",
            str(raw),
            "--population-vcf",
            str(population),
            "--reference-fai",
            str(reference_index),
            "--temp-dir",
            str(scratch),
            "--step-size",
            "10",
            "--threads",
            "4",
            "--workers",
            "2",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    return output, result


def decompressed(path: pathlib.Path) -> str:
    with gzip.open(path, "rt") as input_vcf:
        return input_vcf.read()


def record_body(path: pathlib.Path) -> list[str]:
    return [
        line
        for line in decompressed(path).splitlines()
        if not line.startswith("#")
    ]


def test_fused_transfer_matches_legacy_pipeline(
    tmp_path: pathlib.Path,
) -> None:
    raw, population, reference_index = make_inputs(tmp_path)
    legacy = run_legacy_transfer(
        tmp_path, raw, population, reference_index, step_size=10
    )
    fused, result = run_fused_transfer(
        tmp_path, raw, population, reference_index
    )
    assert result.returncode == 0, result.stderr
    assert pathlib.Path(f"{fused}.tbi").is_file()
    assert record_body(fused) == record_body(legacy)
    assert "workers=1 compression_threads=1" in result.stderr

    _, _, tabix = require_tools()
    for region in ("chr1:1-25", "chr1:10-12", "chrUn:1-7"):
        expected = subprocess.check_output([tabix, str(legacy), region])
        observed = subprocess.check_output([tabix, str(fused), region])
        assert observed == expected


def test_work_items_preserve_legacy_position_ownership() -> None:
    items = build_work_items([("chr1", 25), ("chrUn", 7)], {"chr1"}, 10)
    assert [
        (
            item.number,
            item.contig,
            item.start,
            item.stop,
            item.merge_population,
        )
        for item in items
    ] == [
        (0, "chr1", 1, 10, True),
        (1, "chr1", 11, 20, True),
        (2, "chr1", 21, 25, True),
        (3, "chrUn", 0, 7, False),
    ]


def test_trim_record_matches_number_arg_semantics() -> None:
    line = (
        b"chr1\t2\t.\tAT\tCT,GT\t50\tPASS\t"
        b"AF=0.2,.;IA=7,8;IR=1,2,3;FLAG\t"
        b"GT:FA:FR:FG\t0/1:7,8:1,2,3:0,10,20,30,40,50\n"
    )
    assert trim_record(
        line,
        {b"AF": b"A", b"IA": b"A", b"IR": b"R"},
        {b"FA": b"A", b"FR": b"R", b"FG": b"G"},
    ) == (
        b"chr1\t2\t.\tA\tC\t50\tPASS\t"
        b"AF=0.2;IA=7;IR=1,2;FLAG\t"
        b"GT:FA:FR:FG\t0/1:7:1,2:0,10,20\n"
    )


def test_missing_index_fails_without_replacing_output(
    tmp_path: pathlib.Path,
) -> None:
    raw, population, reference_index = make_inputs(tmp_path)
    pathlib.Path(f"{population}.tbi").unlink()
    output = tmp_path / "fused.vcf.gz"
    output.write_bytes(b"existing output")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = subprocess.run(
        [
            sys.executable,
            str(HYBRID_TRANSFER),
            "--raw-vcf",
            str(raw),
            "--population-vcf",
            str(population),
            "--reference-fai",
            str(reference_index),
            "--temp-dir",
            str(scratch),
            "--threads",
            "2",
            "--workers",
            "1",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "population VCF index is missing" in result.stderr
    assert output.read_bytes() == b"existing output"
    assert not pathlib.Path(f"{output}.tbi").exists()


@pytest.mark.parametrize(
    ("threads", "workers", "message"),
    (
        (0, 1, "--threads must be at least 1"),
        (2, 0, "--workers must be between 1 and 64"),
        (2, 3, "--workers must not exceed --threads"),
    ),
)
def test_invalid_budget_fails_atomically(
    tmp_path: pathlib.Path,
    threads: int,
    workers: int,
    message: str,
) -> None:
    raw, population, reference_index = make_inputs(tmp_path)
    output = tmp_path / "fused.vcf.gz"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = subprocess.run(
        [
            sys.executable,
            str(HYBRID_TRANSFER),
            "--raw-vcf",
            str(raw),
            "--population-vcf",
            str(population),
            "--reference-fai",
            str(reference_index),
            "--temp-dir",
            str(scratch),
            "--threads",
            str(threads),
            "--workers",
            str(workers),
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert message in result.stderr
    assert not output.exists()


def test_command_builder_uses_one_process(tmp_path: pathlib.Path) -> None:
    output = tmp_path / "output.vcf.gz"
    raw = tmp_path / "raw.vcf.gz"
    population = tmp_path / "population.vcf.gz"
    reference_index = tmp_path / "reference.fa.fai"
    scratch = tmp_path / "scratch"
    script = tmp_path / "hybrid_transfer.py"
    pipeline = cmd_pyexec_hybrid_transfer(
        output,
        raw,
        population,
        reference_index,
        scratch,
        script,
        128,
        32,
    )
    assert len(pipeline.nodes) == 1
    command = pipeline.nodes[0]
    assert command.executable == sys.executable
    assert command.args == [
        str(script),
        "--raw-vcf",
        str(raw),
        "--population-vcf",
        str(population),
        "--reference-fai",
        str(reference_index),
        "--temp-dir",
        str(scratch),
        "--threads",
        "128",
        "--workers",
        "32",
        str(output),
    ]
