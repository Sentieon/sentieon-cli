#!/usr/bin/env python
from __future__ import annotations

import argparse
import bisect
import collections
import multiprocessing as mp
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
from typing import BinaryIO, Iterable, Iterator, Sequence

import vcflib
from vcflib import bgzf, tabix


EXTRA_HEADERS = (
    '##INFO=<ID=LHC,Number=1,Type=Integer,Description="Longread hap count">',
)

DEFAULT_STEP_SIZE = 10_000_000
MAX_ANNOTATION_WORKERS = 96
MAX_COMPRESSION_THREADS = 32

_worker_vcf_path: str | None = None
_worker_bed_intervals: dict[str, list[tuple[int, int, int]]] | None = None
_worker_temp_dir: str | None = None
_worker_vcf: BinaryIO | None = None
_worker_index: tabix.Tabix | None = None


class HybridAnnoError(RuntimeError):
    pass


def load_bed(path: pathlib.Path) -> dict[str, list[tuple[int, int, int]]]:
    if not path.exists():
        raise HybridAnnoError(f"input bed file {path} does not exist")

    intervals: dict[str, list[tuple[int, int, int]]] = {}
    with path.open() as bed_file:
        for line_number, line in enumerate(bed_file, 1):
            columns = line.rstrip().split("\t")
            if len(columns) < 4:
                raise HybridAnnoError(
                    f"wrong format in {path} line {line_number}: {line.rstrip()}"
                )
            try:
                start = int(columns[1])
                end = int(columns[2])
                count = int(columns[3])
            except ValueError as error:
                raise HybridAnnoError(
                    f"non-integer BED value in {path} line {line_number}"
                ) from error
            if start < 0 or end < start:
                raise HybridAnnoError(
                    f"invalid BED interval in {path} line {line_number}"
                )
            interval = (start, end, count)
            contig_intervals = intervals.setdefault(columns[0], [])
            if contig_intervals and interval < contig_intervals[-1]:
                raise HybridAnnoError(
                    f"BED intervals are not sorted in {path} line {line_number}"
                )
            contig_intervals.append(interval)
    return intervals


def copy_header_lines(input_vcf: vcflib.VCF) -> list[str]:
    headers: collections.OrderedDict[
        str, collections.OrderedDict[str | None, str]
    ] = collections.OrderedDict()
    pattern = re.compile(r"^##([^=]+)=(<ID=([^,]+).*>)?")

    for line in input_vcf.headers:
        match = pattern.match(line)
        if match is None:
            field, identifier = line, None
        else:
            field, identifier = match.group(1), match.group(3)
        headers.setdefault(field, collections.OrderedDict())[identifier] = line

    for line in EXTRA_HEADERS:
        match = pattern.match(line)
        if match is None:
            continue
        field, identifier = match.group(1), match.group(3)
        headers.setdefault(field, collections.OrderedDict())[identifier] = line

    output = [
        line
        for field_headers in headers.values()
        for line in field_headers.values()
        if not line.startswith("#CHROM")
    ]
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if input_vcf.samples:
        columns.append("FORMAT")
        columns.extend(input_vcf.samples)
    output.append("\t".join(columns))
    return output


def load_vcf_metadata(
    path: pathlib.Path,
) -> tuple[bytes, list[tuple[str, int]]]:
    if not path.exists():
        raise HybridAnnoError(f"input file {path} does not exist")
    if not str(path).endswith(".gz"):
        raise HybridAnnoError("input VCF must be BGZF-compressed")
    if not pathlib.Path(f"{path}.tbi").exists() and not pathlib.Path(
        f"{path}.csi"
    ).exists():
        raise HybridAnnoError(f"input VCF index is missing for {path}")

    input_vcf = vcflib.VCF(str(path), "r")
    try:
        contigs: list[tuple[str, int]] = []
        for contig, description in input_vcf.contigs.items():
            if "length" not in description:
                raise HybridAnnoError(
                    f"contig {contig} has no length in the VCF header"
                )
            try:
                length = int(description["length"])
            except ValueError as error:
                raise HybridAnnoError(
                    f"contig {contig} has an invalid length"
                ) from error
            if length <= 0:
                raise HybridAnnoError(
                    f"contig {contig} has a non-positive length"
                )
            contigs.append((contig, length))
        if not contigs:
            raise HybridAnnoError("input VCF declares no contigs")
        header = ("\n".join(copy_header_lines(input_vcf)) + "\n").encode()
    finally:
        input_vcf.close()
    return header, contigs


def cut_shards(
    contigs: Sequence[tuple[str, int]], step_size: int
) -> list[tuple[int, str, int, int]]:
    shards: list[tuple[int, str, int, int]] = []
    filled = 0
    for contig, length in contigs:
        start = 0
        while start < length:
            chunk_size = min(length - start, step_size - filled)
            end = start + chunk_size
            shards.append((len(shards), contig, start, end))
            start = end
            filled += chunk_size
            if filled == step_size:
                filled = 0
    return shards


def init_worker(
    vcf_path: str,
    bed_intervals: dict[str, list[tuple[int, int, int]]],
    temp_dir: str,
) -> None:
    global _worker_vcf_path
    global _worker_bed_intervals
    global _worker_temp_dir
    global _worker_vcf
    global _worker_index

    _worker_vcf_path = vcf_path
    _worker_bed_intervals = bed_intervals
    _worker_temp_dir = temp_dir
    _worker_vcf = bgzf.open(vcf_path, "rb")
    _worker_index = tabix.Tabix(vcf_path, "r")


def close_worker() -> None:
    global _worker_vcf
    if _worker_vcf is not None:
        _worker_vcf.close()
    _worker_vcf = None


def parse_record(line: bytes) -> tuple[list[bytes], str, int, int]:
    stripped = line.rstrip()
    columns = stripped.split(b"\t", 8)
    if len(columns) < 8:
        raise HybridAnnoError(
            f"VCF record has fewer than eight columns: {stripped!r}"
        )
    try:
        contig = columns[0].decode()
        position = int(columns[1]) - 1
    except (UnicodeDecodeError, ValueError) as error:
        raise HybridAnnoError(f"invalid VCF record: {stripped!r}") from error
    if position < 0:
        raise HybridAnnoError(f"invalid VCF position: {stripped!r}")

    record_end = position + len(columns[3])
    for item in columns[7].split(b";"):
        if item.startswith(b"END="):
            try:
                record_end = int(item[4:])
            except ValueError as error:
                raise HybridAnnoError(
                    f"invalid INFO/END in VCF record: {stripped!r}"
                ) from error
    return columns, contig, position, record_end


def iter_shard_records(
    contig: str, start: int, end: int
) -> Iterator[tuple[list[bytes], int]]:
    if _worker_vcf is None or _worker_index is None:
        raise HybridAnnoError("annotation worker was not initialized")

    ranges = list(_worker_index.query(contig, start, end))
    sequential = False
    for range_start, range_end in ranges:
        _worker_vcf.seek(range_start)
        while True:
            line = _worker_vcf.readline()
            if not line:
                return
            if line.startswith(b"#"):
                continue
            columns, record_contig, position, record_end = parse_record(line)
            if record_contig != contig or position >= end:
                return
            if not sequential and record_end <= start:
                if _worker_vcf.tell() >= range_end:
                    break
                continue
            sequential = True
            yield columns, position
        if sequential:
            return


def annotate_shard(
    shard: tuple[int, str, int, int],
) -> tuple[int, str, int, int]:
    shard_number, contig, start, end = shard
    if _worker_bed_intervals is None or _worker_temp_dir is None:
        raise HybridAnnoError("annotation worker was not initialized")

    output_path = os.path.join(
        _worker_temp_dir, f"fragment.{shard_number:08d}.vcf"
    )
    contig_intervals = _worker_bed_intervals.get(contig)
    bed_index = 0
    bed_initialized = False
    record_count = 0
    annotated_count = 0

    with open(output_path, "wb") as output:
        for columns, position in iter_shard_records(contig, start, end):
            hap_count = -1
            if contig_intervals:
                if not bed_initialized:
                    bed_index = bisect.bisect_left(
                        contig_intervals, (position, position + 1, 0)
                    )
                    if bed_index > 0:
                        bed_index -= 1
                    bed_initialized = True
                while (
                    bed_index < len(contig_intervals)
                    and contig_intervals[bed_index][1] <= position
                ):
                    bed_index += 1
                if (
                    bed_index < len(contig_intervals)
                    and contig_intervals[bed_index][0] <= position
                ):
                    hap_count = contig_intervals[bed_index][2]

            if position < start:
                continue
            if hap_count != -1:
                columns[7] += f";LHC={hap_count}".encode()
                annotated_count += 1
            output.write(b"\t".join(columns))
            output.write(b"\n")
            record_count += 1

    return shard_number, output_path, record_count, annotated_count


def generate_fragments(
    shards: Sequence[tuple[int, str, int, int]],
    workers: int,
    vcf_path: pathlib.Path,
    bed_intervals: dict[str, list[tuple[int, int, int]]],
    temp_dir: pathlib.Path,
) -> tuple[list[pathlib.Path], int, int]:
    results: Iterable[tuple[int, str, int, int]]
    initializer_arguments = (str(vcf_path), bed_intervals, str(temp_dir))

    if workers == 1:
        init_worker(*initializer_arguments)
        try:
            results = [annotate_shard(shard) for shard in shards]
        finally:
            close_worker()
    else:
        context = mp.get_context()
        with context.Pool(
            processes=workers,
            initializer=init_worker,
            initargs=initializer_arguments,
        ) as pool:
            results = pool.imap(annotate_shard, shards, chunksize=1)
            results = list(results)

    ordered = sorted(results)
    expected_numbers = list(range(len(shards)))
    observed_numbers = [result[0] for result in ordered]
    if observed_numbers != expected_numbers:
        raise HybridAnnoError("annotation shards were not returned exactly once")
    paths = [pathlib.Path(result[1]) for result in ordered]
    record_count = sum(result[2] for result in ordered)
    annotated_count = sum(result[3] for result in ordered)
    return paths, record_count, annotated_count


def compress_fragments(
    bgzip_path: str,
    header: bytes,
    fragments: Sequence[pathlib.Path],
    threads: int,
    output_path: pathlib.Path,
) -> None:
    with output_path.open("wb") as compressed_output:
        process = subprocess.Popen(
            [bgzip_path, "-@", str(threads), "-c"],
            stdin=subprocess.PIPE,
            stdout=compressed_output,
            stderr=subprocess.PIPE,
        )
        if process.stdin is None or process.stderr is None:
            process.kill()
            raise HybridAnnoError("failed to open bgzip pipes")
        try:
            process.stdin.write(header)
            for fragment in fragments:
                with fragment.open("rb") as fragment_file:
                    shutil.copyfileobj(
                        fragment_file, process.stdin, length=16 * 1024 * 1024
                    )
            process.stdin.close()
        except BaseException:
            process.kill()
            process.wait()
            raise
        error_output = process.stderr.read().decode(errors="replace")
        return_code = process.wait()
    if return_code != 0:
        raise HybridAnnoError(
            f"bgzip failed with exit code {return_code}: {error_output.strip()}"
        )


def publish_output(
    output_path: pathlib.Path,
    header: bytes,
    fragments: Sequence[pathlib.Path],
    compression_threads: int,
    bgzip_path: str,
    tabix_path: str,
) -> None:
    output_parent = output_path.parent
    if not output_parent.is_dir():
        raise HybridAnnoError(
            f"output directory {output_parent} does not exist"
        )

    descriptor, partial_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.partial.",
        suffix=".vcf.gz",
        dir=output_parent,
    )
    os.close(descriptor)
    partial_path = pathlib.Path(partial_name)
    partial_index = pathlib.Path(f"{partial_path}.tbi")
    try:
        compress_fragments(
            bgzip_path,
            header,
            fragments,
            compression_threads,
            partial_path,
        )
        result = subprocess.run(
            [tabix_path, "-f", "-p", "vcf", str(partial_path)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise HybridAnnoError(
                "tabix failed with exit code "
                f"{result.returncode}: {result.stderr.strip()}"
            )
        os.replace(partial_path, output_path)
        os.replace(partial_index, pathlib.Path(f"{output_path}.tbi"))
        stale_csi = pathlib.Path(f"{output_path}.csi")
        if stale_csi.exists():
            stale_csi.unlink()
    finally:
        for temporary_path in (partial_path, partial_index):
            if temporary_path.exists():
                temporary_path.unlink()


def main(args: argparse.Namespace) -> int:
    if args.threads < 1:
        raise HybridAnnoError("--threads must be at least 1")
    if args.step_size is not None and args.step_size < 1:
        raise HybridAnnoError("--step-size must be at least 1")

    input_path = pathlib.Path(args.vcf).resolve()
    output_path = pathlib.Path(args.output).resolve()
    bed_path = pathlib.Path(args.bed).resolve()
    if input_path == output_path:
        raise HybridAnnoError("input and output VCF paths must differ")
    if not str(output_path).endswith(".gz"):
        raise HybridAnnoError("output VCF must end in .gz")

    bgzip_path = shutil.which("bgzip")
    tabix_path = shutil.which("tabix")
    if bgzip_path is None:
        raise HybridAnnoError("bgzip is not available on PATH")
    if tabix_path is None:
        raise HybridAnnoError("tabix is not available on PATH")

    bed_intervals = load_bed(bed_path)
    header, contigs = load_vcf_metadata(input_path)
    workers = min(args.threads, MAX_ANNOTATION_WORKERS)
    compression_threads = min(args.threads, MAX_COMPRESSION_THREADS)
    step_size = args.step_size or DEFAULT_STEP_SIZE
    shards = cut_shards(contigs, step_size)

    temp_dir = pathlib.Path(
        tempfile.mkdtemp(
            prefix="hybrid-anno.",
            dir=os.getenv("SENTIEON_TMPDIR"),
        )
    )
    success = False
    try:
        print(
            "hybrid_anno: "
            f"workers={workers} compression_threads={compression_threads} "
            f"step_size={step_size} shards={len(shards)} "
            f"tmpdir={temp_dir}",
            file=sys.stderr,
        )
        fragments, record_count, annotated_count = generate_fragments(
            shards,
            workers,
            input_path,
            bed_intervals,
            temp_dir,
        )
        publish_output(
            output_path,
            header,
            fragments,
            compression_threads,
            bgzip_path,
            tabix_path,
        )
        print(
            "hybrid_anno: "
            f"records={record_count} annotated={annotated_count}",
            file=sys.stderr,
        )
        success = True
    finally:
        if success:
            shutil.rmtree(temp_dir)
        else:
            print(
                f"hybrid_anno: preserved failed work directory {temp_dir}",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="sentieon pyexec hybrid_anno.py",
        usage="%(prog)s [options] -v VCF -b BED output",
    )
    parser.add_argument("output", help="Output vcf file name")
    parser.add_argument(
        "-v", "--vcf", required=True, help="Input vcf file name"
    )
    parser.add_argument(
        "-b",
        "--bed",
        required=True,
        help="region haplotype count from long reads",
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
        default=None,
        help=argparse.SUPPRESS,
    )
    try:
        sys.exit(main(parser.parse_args()))
    except HybridAnnoError as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
