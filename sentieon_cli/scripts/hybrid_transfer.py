#!/usr/bin/env python
"""Bounded, ordered population annotation transfer for Hybrid VCFs."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import BinaryIO, Iterator, Sequence, cast

DEFAULT_STEP_SIZE = 10_000_000
DEFAULT_WORKERS = 32
MAX_WORKERS = 64

_HEADER_FIELD_PATTERN = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')

_worker_raw_vcf: str | None = None
_worker_population_vcf: str | None = None
_worker_merge_rules: str | None = None
_worker_info_numbers: dict[bytes, bytes] | None = None
_worker_format_numbers: dict[bytes, bytes] | None = None
_worker_temp_dir: str | None = None


class HybridTransferError(RuntimeError):
    """A population-transfer contract or subprocess failed."""


@dataclass(frozen=True)
class WorkItem:
    """One legacy-compatible reference shard or raw-only contig."""

    number: int
    contig: str
    start: int
    stop: int
    merge_population: bool


@dataclass(frozen=True)
class WorkResult:
    """Header/body fragment produced by one ordered work item."""

    number: int
    contig: str
    header_path: str
    body_path: str
    input_records: int
    output_records: int


def require_index(path: pathlib.Path, label: str) -> None:
    """Require a tabix or CSI index without attempting discovery/fallback."""

    if not path.is_file():
        raise HybridTransferError(f"{label} {path} does not exist")
    tbi = pathlib.Path(f"{path}.tbi")
    csi = pathlib.Path(f"{path}.csi")
    if not tbi.is_file() and not csi.is_file():
        raise HybridTransferError(f"{label} index is missing for {path}")


def run_checked(command: Sequence[str], label: str) -> bytes:
    """Run a bounded metadata command and return stdout."""

    try:
        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError as error:
        raise HybridTransferError(
            f"required executable {command[0]} is missing"
        ) from error
    if result.returncode != 0:
        message = result.stderr.decode(errors="replace").strip()
        raise HybridTransferError(
            f"{label} failed with exit code {result.returncode}: {message}"
        )
    return result.stdout


def parse_structured_header(line: str) -> dict[str, str]:
    """Parse the simple structured-header grammar used by the legacy code."""

    try:
        start = line.index("<")
        end = line.index(">")
    except ValueError as error:
        raise HybridTransferError(
            f"malformed structured VCF header: {line}"
        ) from error
    return dict(_HEADER_FIELD_PATTERN.findall(line[start + 1 : end]))


def population_contract(
    raw_vcf: pathlib.Path, population_vcf: pathlib.Path
) -> tuple[str, set[str], bytes]:
    """Return legacy merge rules, population contigs, and merged header."""

    population_header = run_checked(
        ["bcftools", "view", "--no-version", "-h", str(population_vcf)],
        "population header read",
    )
    info_fields: list[str] = []
    population_contigs: set[str] = set()
    for raw_line in population_header.decode(errors="strict").splitlines():
        if raw_line.startswith("##INFO") and ",Number=A" in raw_line:
            fields = parse_structured_header(raw_line)
            if "ID" not in fields:
                raise HybridTransferError(
                    f"population INFO header has no ID: {raw_line}"
                )
            info_fields.append(fields["ID"])
        elif raw_line.startswith("##contig"):
            fields = parse_structured_header(raw_line)
            if "ID" not in fields:
                raise HybridTransferError(
                    f"population contig header has no ID: {raw_line}"
                )
            population_contigs.add(fields["ID"])

    merge_rules = ",".join(f"{field}:sum" for field in info_fields)
    merged_header = run_checked(
        [
            "bcftools",
            "merge",
            "--print-header",
            "--no-version",
            "--regions-overlap",
            "pos",
            "-m",
            "all",
            "-i",
            merge_rules,
            str(raw_vcf),
            str(population_vcf),
        ],
        "merged header construction",
    )
    if not merged_header.endswith(b"\n"):
        raise HybridTransferError(
            "merged VCF header is not newline terminated"
        )
    return merge_rules, population_contigs, merged_header


def parse_number_schemas(
    merged_header: bytes,
) -> tuple[dict[bytes, bytes], dict[bytes, bytes]]:
    """Collect Number=A/R/G INFO and FORMAT declarations for trimming."""

    infos: dict[bytes, bytes] = {}
    formats: dict[bytes, bytes] = {}
    for line in merged_header.decode(errors="strict").splitlines():
        target: dict[bytes, bytes] | None = None
        if line.startswith("##INFO"):
            target = infos
        elif line.startswith("##FORMAT"):
            target = formats
        if target is None:
            continue
        fields = parse_structured_header(line)
        number = fields.get("Number")
        identifier = fields.get("ID")
        if identifier is not None and number in ("A", "R", "G"):
            target[identifier.encode()] = number.encode()
    return infos, formats


def load_fai(path: pathlib.Path) -> list[tuple[str, int]]:
    """Load reference contigs in their declared order."""

    if not path.is_file():
        raise HybridTransferError(f"reference index {path} does not exist")
    contigs: list[tuple[str, int]] = []
    seen: set[str] = set()
    with path.open() as reference_index:
        for line_number, line in enumerate(reference_index, 1):
            columns = line.rstrip().split("\t")
            if len(columns) < 2:
                raise HybridTransferError(
                    f"reference index line {line_number} has fewer than "
                    "two fields"
                )
            contig = columns[0]
            try:
                length = int(columns[1])
            except ValueError as error:
                raise HybridTransferError(
                    f"reference index line {line_number} has an invalid length"
                ) from error
            if contig in seen:
                raise HybridTransferError(
                    f"reference index contains duplicate contig {contig}"
                )
            if length <= 0:
                raise HybridTransferError(
                    f"reference index line {line_number} has a "
                    "non-positive length"
                )
            seen.add(contig)
            contigs.append((contig, length))
    if not contigs:
        raise HybridTransferError("reference index declares no contigs")
    return contigs


def build_work_items(
    contigs: Sequence[tuple[str, int]],
    population_contigs: set[str],
    step_size: int,
) -> list[WorkItem]:
    """Reproduce the legacy 10 Mb shard and unusual-contig ownership."""

    items: list[WorkItem] = []
    for contig, length in contigs:
        if contig not in population_contigs:
            items.append(
                WorkItem(
                    number=len(items),
                    contig=contig,
                    start=0,
                    stop=length,
                    merge_population=False,
                )
            )
            continue

        start = 1
        while start <= length:
            stop = min(start + step_size - 1, length)
            items.append(
                WorkItem(
                    number=len(items),
                    contig=contig,
                    start=start,
                    stop=stop,
                    merge_population=True,
                )
            )
            start = stop + 1
    return items


def parse_key_value(value: bytes) -> tuple[bytes, bytes | None]:
    """Parse one INFO key/value, representing flags with ``None``."""

    fields = value.split(b"=", 1)
    if len(fields) == 2:
        return fields[0], fields[1]
    return fields[0], None


def filter_values(values: bytes, keep: Sequence[bool]) -> bytes:
    """Filter a comma-delimited VCF vector using legacy zip truncation."""

    return b",".join(
        value for value, retained in zip(values.split(b","), keep) if retained
    )


def trim_record(
    line: bytes,
    info_numbers: dict[bytes, bytes],
    format_numbers: dict[bytes, bytes],
) -> bytes | None:
    """Port ``trimalt.py`` record semantics without object construction."""

    values = line.rstrip().split(b"\t")
    if len(values) < 9 or values[8] == b".":
        return None

    alts = values[4].split(b",")
    keep_alt: list[bool] | None = None
    info: list[tuple[bytes, bytes | None]] | None = None

    if b"<NON_REF>" in alts:
        nonref_index = alts.index(b"<NON_REF>")
        keep_alt = [True] * (nonref_index + 1)
        keep_alt.extend([False] * (len(alts) - nonref_index - 1))
    elif len(values) >= 8 and values[7] != b".":
        info = [parse_key_value(item) for item in values[7].split(b";")]
        for key, value in info:
            if key == b"AF":
                if value is None:
                    raise HybridTransferError("INFO/AF is declared as a flag")
                keep_alt = [item != b"." for item in value.split(b",")]
                break

    if keep_alt is None:
        return None
    if all(keep_alt):
        return line

    keep_ref = [True, *keep_alt]
    if info is None and len(values) >= 8 and values[7] != b".":
        info = [parse_key_value(item) for item in values[7].split(b";")]

    if values[4] != b".":
        retained_alts = [
            alt for alt, retained in zip(alts, keep_alt) if retained
        ]
        reference = values[3]
        try:
            suffix_length = (
                min(
                    len(reference),
                    min(
                        len(alt)
                        for alt in retained_alts
                        if alt != b"<NON_REF>"
                    ),
                )
                - 1
            )
        except ValueError:
            suffix_length = 0
        if suffix_length > 0:
            reference = reference[:-suffix_length]
            retained_alts = [
                (alt if alt == b"<NON_REF>" else alt[:-suffix_length])
                for alt in retained_alts
            ]
        values[3] = reference
        values[4] = b",".join(retained_alts)

    if info is not None:
        rewritten_info: list[bytes] = []
        for key, value in info:
            number = info_numbers.get(key) if value != b"." else None
            if number == b"A":
                if value is None:
                    raise HybridTransferError(
                        "Number=A INFO/"
                        f"{key.decode(errors='replace')} is a flag"
                    )
                rewritten_info.append(
                    key + b"=" + filter_values(value, keep_alt)
                )
            elif number == b"R":
                if value is None:
                    raise HybridTransferError(
                        "Number=R INFO/"
                        f"{key.decode(errors='replace')} is a flag"
                    )
                rewritten_info.append(
                    key + b"=" + filter_values(value, keep_ref)
                )
            elif value is None:
                rewritten_info.append(key)
            else:
                rewritten_info.append(key + b"=" + value)
        values[7] = b";".join(rewritten_info)

    if len(values) >= 10 and values[9] != b".":
        format_keys = values[8].split(b":")
        diploid_genotypes = [
            first and second
            for index, first in enumerate(keep_ref)
            for second in keep_ref[: index + 1]
        ]
        for sample_index, sample in enumerate(values[9:], 9):
            sample_values = sample.split(b":")
            for field_index, value in enumerate(sample_values):
                if field_index >= len(format_keys):
                    break
                number = (
                    format_numbers.get(format_keys[field_index])
                    if value != b"."
                    else None
                )
                if number == b"A":
                    sample_values[field_index] = filter_values(value, keep_alt)
                elif number == b"R":
                    sample_values[field_index] = filter_values(value, keep_ref)
                elif number == b"G":
                    vector_length = len(value.split(b","))
                    if vector_length == len(keep_ref):
                        keep_genotype = keep_ref
                    elif vector_length == len(diploid_genotypes):
                        keep_genotype = diploid_genotypes
                    else:
                        name = format_keys[field_index].decode(
                            errors="replace"
                        )
                        raise HybridTransferError(
                            f"FORMAT/{name} Number=G vector has "
                            "unexpected length"
                        )
                    sample_values[field_index] = filter_values(
                        value, keep_genotype
                    )
            values[sample_index] = b":".join(sample_values)

    return b"\t".join(values) + b"\n"


def init_worker(
    raw_vcf: str,
    population_vcf: str,
    merge_rules: str,
    info_numbers: dict[bytes, bytes],
    format_numbers: dict[bytes, bytes],
    temp_dir: str,
) -> None:
    """Initialize immutable process-local transfer state."""

    global _worker_raw_vcf
    global _worker_population_vcf
    global _worker_merge_rules
    global _worker_info_numbers
    global _worker_format_numbers
    global _worker_temp_dir

    _worker_raw_vcf = raw_vcf
    _worker_population_vcf = population_vcf
    _worker_merge_rules = merge_rules
    _worker_info_numbers = info_numbers
    _worker_format_numbers = format_numbers
    _worker_temp_dir = temp_dir


def worker_command(item: WorkItem, region_path: pathlib.Path) -> list[str]:
    """Build the exact legacy merge or unusual-contig view primitive."""

    if _worker_raw_vcf is None:
        raise HybridTransferError("transfer worker was not initialized")
    if item.merge_population:
        if _worker_population_vcf is None or _worker_merge_rules is None:
            raise HybridTransferError("transfer worker was not initialized")
        return [
            "bcftools",
            "merge",
            "--regions-file",
            str(region_path),
            "--no-version",
            "--regions-overlap",
            "pos",
            "-m",
            "all",
            "-i",
            _worker_merge_rules,
            _worker_raw_vcf,
            _worker_population_vcf,
        ]
    return [
        "bcftools",
        "view",
        "--no-version",
        "--regions-file",
        str(region_path),
        _worker_raw_vcf,
    ]


def consume_records(
    stream: BinaryIO,
    header: BinaryIO,
    body: BinaryIO,
    trim: bool,
) -> tuple[int, int]:
    """Split one bcftools stream and optionally apply exact trimming."""

    if _worker_info_numbers is None or _worker_format_numbers is None:
        raise HybridTransferError("transfer worker was not initialized")
    input_records = 0
    output_records = 0
    saw_column_header = False
    for line in stream:
        if line.startswith(b"#"):
            if input_records:
                raise HybridTransferError(
                    "VCF header appeared after record data"
                )
            header.write(line)
            if line.startswith(b"#CHROM\t"):
                saw_column_header = True
            continue
        if not saw_column_header:
            raise HybridTransferError("bcftools output has no #CHROM header")
        input_records += 1
        output_line = (
            trim_record(line, _worker_info_numbers, _worker_format_numbers)
            if trim
            else line
        )
        if output_line is not None:
            body.write(output_line)
            output_records += 1
    if not saw_column_header:
        raise HybridTransferError("bcftools output has no #CHROM header")
    return input_records, output_records


def process_work_item(item: WorkItem) -> WorkResult:
    """Run one merge/view shard and emit an uncompressed ordered fragment."""

    if _worker_temp_dir is None:
        raise HybridTransferError("transfer worker was not initialized")
    base = pathlib.Path(_worker_temp_dir)
    prefix = f"fragment.{item.number:08d}"
    region_path = base / f"{prefix}.bed"
    header_path = base / f"{prefix}.header.vcf"
    body_path = base / f"{prefix}.body.vcf"
    error_path = base / f"{prefix}.stderr"
    region_path.write_text(f"{item.contig}\t{item.start}\t{item.stop}\n")

    command = worker_command(item, region_path)
    try:
        with error_path.open("wb") as error_file:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=error_file,
            )
            if process.stdout is None:
                process.terminate()
                raise HybridTransferError(
                    "bcftools stdout pipe was not created"
                )
            process_stdout = cast(BinaryIO, process.stdout)
            try:
                with (
                    header_path.open("wb") as header,
                    body_path.open("wb") as body,
                ):
                    input_records, output_records = consume_records(
                        process_stdout,
                        header,
                        body,
                        trim=item.merge_population,
                    )
            except BaseException:
                process.terminate()
                process.wait()
                raise
            finally:
                process_stdout.close()
            return_code = process.wait()
    except FileNotFoundError as error:
        raise HybridTransferError(
            "required executable bcftools is missing"
        ) from error

    if return_code != 0:
        message = error_path.read_text(errors="replace").strip()
        raise HybridTransferError(
            f"work item {item.number} ({item.contig}) failed with exit code "
            f"{return_code}: {message}"
        )
    return WorkResult(
        number=item.number,
        contig=item.contig,
        header_path=str(header_path),
        body_path=str(body_path),
        input_records=input_records,
        output_records=output_records,
    )


def column_header(path: pathlib.Path) -> bytes:
    """Return the #CHROM line for a fragment header."""

    with path.open("rb") as header:
        for line in header:
            if line.startswith(b"#CHROM\t"):
                return line
    raise HybridTransferError(f"fragment header {path} has no #CHROM line")


def copy_path(path: pathlib.Path, output: BinaryIO) -> None:
    """Copy a fragment with a large bounded userspace buffer."""

    with path.open("rb") as fragment:
        shutil.copyfileobj(fragment, output, length=16 * 1024 * 1024)


def partial_paths(
    output_path: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path]:
    """Allocate a unique sibling path accepted by bcftools write-index."""

    descriptor, partial_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.partial.",
        suffix=".vcf.gz",
        dir=output_path.parent,
    )
    os.close(descriptor)
    partial_path = pathlib.Path(partial_name)
    partial_path.unlink()
    return partial_path, pathlib.Path(f"{partial_path}.tbi")


def publish_pair(
    partial_path: pathlib.Path,
    partial_index: pathlib.Path,
    output_path: pathlib.Path,
) -> None:
    """Replace the completed VCF/index pair only after both validate."""

    if not partial_path.is_file() or not partial_index.is_file():
        raise HybridTransferError("bcftools did not create the VCF/index pair")
    run_checked(
        ["bcftools", "index", "--nrecords", str(partial_path)],
        "output index validation",
    )
    output_index = pathlib.Path(f"{output_path}.tbi")
    backup_vcf = pathlib.Path(f"{partial_path}.previous-vcf")
    backup_index = pathlib.Path(f"{partial_path}.previous-index")
    had_vcf = output_path.exists()
    had_index = output_index.exists()
    try:
        if had_vcf:
            os.replace(output_path, backup_vcf)
        if had_index:
            os.replace(output_index, backup_index)
        os.replace(partial_path, output_path)
        os.replace(partial_index, output_index)
    except BaseException:
        if output_path.exists():
            output_path.unlink()
        if output_index.exists():
            output_index.unlink()
        if had_vcf and backup_vcf.exists():
            os.replace(backup_vcf, output_path)
        if had_index and backup_index.exists():
            os.replace(backup_index, output_index)
        raise
    finally:
        if backup_vcf.exists():
            backup_vcf.unlink()
        if backup_index.exists():
            backup_index.unlink()


def stream_results(
    results: Iterator[WorkResult],
    output_path: pathlib.Path,
    compression_threads: int,
    temp_dir: pathlib.Path,
) -> tuple[int, int]:
    """Publish ordered fragment bodies through one bcftools BGZF writer."""

    partial_path, partial_index = partial_paths(output_path)
    publisher_error_path = temp_dir / "publisher.stderr"
    publisher: subprocess.Popen[bytes] | None = None
    expected_number = 0
    expected_column_header: bytes | None = None
    input_records = 0
    output_records = 0
    try:
        with publisher_error_path.open("wb") as publisher_error:
            publisher = subprocess.Popen(
                [
                    "bcftools",
                    "view",
                    "--no-version",
                    "-O",
                    "z",
                    "--threads",
                    str(compression_threads),
                    "-o",
                    str(partial_path),
                    "-W=tbi",
                    "-",
                ],
                stdin=subprocess.PIPE,
                stderr=publisher_error,
            )
            if publisher.stdin is None:
                publisher.terminate()
                raise HybridTransferError(
                    "bcftools stdin pipe was not created"
                )
            publisher_stdin = cast(BinaryIO, publisher.stdin)
            try:
                for result in results:
                    if result.number != expected_number:
                        raise HybridTransferError(
                            "transfer work results are not complete and "
                            "ordered"
                        )
                    header_path = pathlib.Path(result.header_path)
                    body_path = pathlib.Path(result.body_path)
                    observed_column_header = column_header(header_path)
                    if expected_column_header is None:
                        expected_column_header = observed_column_header
                        copy_path(header_path, publisher_stdin)
                    elif observed_column_header != expected_column_header:
                        raise HybridTransferError(
                            f"fragment {result.number} sample header differs"
                        )
                    copy_path(body_path, publisher_stdin)
                    input_records += result.input_records
                    output_records += result.output_records
                    header_path.unlink()
                    body_path.unlink()
                    expected_number += 1
            except BaseException:
                publisher_stdin.close()
                publisher.terminate()
                publisher.wait()
                raise
            publisher_stdin.close()
            return_code = publisher.wait()

        if expected_number == 0:
            raise HybridTransferError("transfer constructed no work items")
        if return_code != 0:
            message = publisher_error_path.read_text(errors="replace").strip()
            raise HybridTransferError(
                f"final VCF publication failed with exit code {return_code}: "
                f"{message}"
            )
        publish_pair(partial_path, partial_index, output_path)
        return input_records, output_records
    except FileNotFoundError as error:
        raise HybridTransferError(
            "required executable bcftools is missing"
        ) from error
    finally:
        if publisher is not None and publisher.poll() is None:
            publisher.terminate()
            publisher.wait()
        if partial_path.exists():
            partial_path.unlink()
        if partial_index.exists():
            partial_index.unlink()


def run_transfer(args: argparse.Namespace) -> int:
    """Validate, execute, and atomically publish the fused transfer."""

    if args.threads < 1:
        raise HybridTransferError("--threads must be at least 1")
    if args.workers < 1 or args.workers > MAX_WORKERS:
        raise HybridTransferError(
            f"--workers must be between 1 and {MAX_WORKERS}"
        )
    if args.workers > args.threads:
        raise HybridTransferError("--workers must not exceed --threads")
    if args.step_size < 1:
        raise HybridTransferError("--step-size must be at least 1")

    raw_vcf = pathlib.Path(args.raw_vcf).resolve()
    population_vcf = pathlib.Path(args.population_vcf).resolve()
    reference_fai = pathlib.Path(args.reference_fai).resolve()
    output_path = pathlib.Path(args.output).resolve()
    base_temp_dir = pathlib.Path(args.temp_dir).resolve()
    require_index(raw_vcf, "raw VCF")
    require_index(population_vcf, "population VCF")
    if not output_path.parent.is_dir():
        raise HybridTransferError(
            f"output directory {output_path.parent} does not exist"
        )
    if not base_temp_dir.is_dir():
        raise HybridTransferError(
            f"temporary directory {base_temp_dir} does not exist"
        )

    merge_rules, population_contigs, merged_header = population_contract(
        raw_vcf, population_vcf
    )
    info_numbers, format_numbers = parse_number_schemas(merged_header)
    contigs = load_fai(reference_fai)
    items = build_work_items(contigs, population_contigs, args.step_size)
    worker_budget = max(1, (args.threads - 1) // 2)
    workers = min(args.workers, len(items), worker_budget)
    compression_threads = max(0, args.threads - (2 * workers) - 1)
    temp_dir = pathlib.Path(
        tempfile.mkdtemp(prefix="hybrid-transfer.", dir=base_temp_dir)
    )
    success = False
    pool: mp.pool.Pool | None = None
    try:
        print(
            "hybrid_transfer: "
            f"workers={workers} compression_threads={compression_threads} "
            f"step_size={args.step_size} shards={len(items)} "
            f"tmpdir={temp_dir}",
            file=sys.stderr,
        )
        context = mp.get_context()
        pool = context.Pool(
            processes=workers,
            initializer=init_worker,
            initargs=(
                str(raw_vcf),
                str(population_vcf),
                merge_rules,
                info_numbers,
                format_numbers,
                str(temp_dir),
            ),
        )
        results = pool.imap(process_work_item, items, chunksize=1)
        input_records, output_records = stream_results(
            results,
            output_path,
            compression_threads,
            temp_dir,
        )
        pool.close()
        pool.join()
        pool = None
        print(
            "hybrid_transfer: "
            f"input_records={input_records} output_records={output_records}",
            file=sys.stderr,
        )
        success = True
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if success:
            shutil.rmtree(temp_dir)
        else:
            print(
                f"hybrid_transfer: preserved failed work directory {temp_dir}",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="sentieon pyexec hybrid_transfer.py",
        usage=(
            "%(prog)s --raw-vcf RAW --population-vcf POP "
            "--reference-fai REF.fai --temp-dir DIR [options] output.vcf.gz"
        ),
    )
    parser.add_argument("output", help="output BGZF VCF")
    parser.add_argument(
        "--raw-vcf", required=True, help="annotated Hybrid VCF"
    )
    parser.add_argument(
        "--population-vcf", required=True, help="population annotation VCF"
    )
    parser.add_argument(
        "--reference-fai", required=True, help="reference FASTA index"
    )
    parser.add_argument(
        "--temp-dir", required=True, help="existing scratch directory"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=mp.cpu_count(),
        help="total thread budget",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=DEFAULT_STEP_SIZE,
        help=argparse.SUPPRESS,
    )
    try:
        sys.exit(run_transfer(parser.parse_args()))
    except HybridTransferError as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
