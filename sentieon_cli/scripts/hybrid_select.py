#!/usr/bin/env python
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import pathlib
import shutil
import sys
import tempfile
from typing import BinaryIO, Iterable, Iterator, Sequence

import vcflib
from vcflib import bgzf, tabix

DEFAULT_STEP_SIZE = 10_000_000
DEFAULT_SLOP_SIZE = 1_000
MAX_SELECT_WORKERS = 96

_worker_vcf: BinaryIO | None = None
_worker_index: tabix.Tabix | None = None
_worker_fai_lengths: dict[str, int] | None = None
_worker_temp_dir: str | None = None
_worker_slop_size: int | None = None


class HybridSelectError(RuntimeError):
    pass


def load_fai(path: pathlib.Path) -> dict[str, int]:
    if not path.is_file():
        raise HybridSelectError(f"reference index {path} does not exist")

    lengths: dict[str, int] = {}
    with path.open() as reference_index:
        for line_number, line in enumerate(reference_index, 1):
            columns = line.rstrip().split("\t")
            if len(columns) < 2:
                raise HybridSelectError(
                    f"reference index line {line_number} has fewer than two fields"
                )
            try:
                length = int(columns[1])
            except ValueError as error:
                raise HybridSelectError(
                    f"reference index line {line_number} has an invalid length"
                ) from error
            if length <= 0:
                raise HybridSelectError(
                    f"reference index line {line_number} has a non-positive length"
                )
            if columns[0] in lengths:
                raise HybridSelectError(
                    f"reference index contains duplicate contig {columns[0]}"
                )
            lengths[columns[0]] = length
    if not lengths:
        raise HybridSelectError("reference index declares no contigs")
    return lengths


def load_vcf_contigs(path: pathlib.Path) -> list[tuple[str, int]]:
    if not path.is_file():
        raise HybridSelectError(f"input file {path} does not exist")
    if not str(path).endswith(".gz"):
        raise HybridSelectError("input VCF must be BGZF-compressed")
    if (
        not pathlib.Path(f"{path}.tbi").is_file()
        and not pathlib.Path(f"{path}.csi").is_file()
    ):
        raise HybridSelectError(f"input VCF index is missing for {path}")

    input_vcf = vcflib.VCF(str(path), "r")
    try:
        genotype = input_vcf.formats.get("GT")
        observed_type = None
        if genotype is not None:
            observed_type = {
                key: value
                for key, value in genotype.items()
                if key in ("Number", "Type")
            }
        expected_type = {"Number": "1", "Type": "String"}
        if observed_type != expected_type:
            raise HybridSelectError(
                "VCF FORMAT/GT is not Number=1,Type=String"
            )
        if not input_vcf.samples:
            raise HybridSelectError("input VCF contains no samples")

        contigs: list[tuple[str, int]] = []
        for contig, description in input_vcf.contigs.items():
            if "length" not in description:
                raise HybridSelectError(
                    f"contig {contig} has no length in the VCF header"
                )
            try:
                length = int(description["length"])
            except ValueError as error:
                raise HybridSelectError(
                    f"contig {contig} has an invalid length"
                ) from error
            if length <= 0:
                raise HybridSelectError(
                    f"contig {contig} has a non-positive length"
                )
            contigs.append((contig, length))
        if not contigs:
            raise HybridSelectError("input VCF declares no contigs")
        return contigs
    finally:
        input_vcf.close()


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


def parse_integer_list(value: bytes | None, field: str) -> list[int] | None:
    if value is None or value in (b"", b"."):
        return None
    try:
        return [int(item) for item in value.split(b",")]
    except ValueError as error:
        raise HybridSelectError(
            f"invalid FORMAT/{field} value {value!r}"
        ) from error


def sample_fields(columns: list[bytes]) -> dict[bytes, bytes]:
    if len(columns) < 10:
        raise HybridSelectError("VCF record contains no sample column")
    keys = columns[8].split(b":")
    values = columns[9].split(b":")
    return dict(zip(keys, values))


def has_info_field(info: bytes, name: bytes) -> bool:
    return any(item.split(b"=", 1)[0] == name for item in info.split(b";"))


def record_passes(columns: list[bytes]) -> bool:
    sample = sample_fields(columns)
    lad = parse_integer_list(sample.get(b"LAD"), "LAD") or [0, 0, 0]
    lpl = parse_integer_list(sample.get(b"LPL"), "LPL") or [0, 0, 0]
    spl = parse_integer_list(sample.get(b"SPL"), "SPL") or [0, 0, 0]

    try:
        lminpos = lpl.index(0)
        sminpos = spl.index(0)
        lpl_ref = lpl[0]
        lpl = lpl.copy()
        spl = spl.copy()
        lpl.remove(0)
        spl.remove(0)
        lpl_conf = min(lpl)
        spl_conf = min(spl)
    except (ValueError, IndexError) as error:
        raise HybridSelectError(
            "invalid LPL/SPL genotype likelihood vector"
        ) from error

    if sum(lad) < 2:
        return False

    is_str = has_info_field(columns[7], b"STR")
    if not is_str and lminpos != 0:
        lpl_conf = lpl_ref
    if lpl_conf < 30.0:
        return False

    if spl_conf >= 30.0:
        if lminpos == sminpos:
            return False
        if lminpos == 0 and is_str:
            return False
    return True


def parse_record(line: bytes) -> tuple[list[bytes], str, int, int]:
    stripped = line.rstrip()
    columns = stripped.split(b"\t")
    if len(columns) < 10:
        raise HybridSelectError(
            f"VCF record has fewer than ten columns: {stripped!r}"
        )
    try:
        contig = columns[0].decode()
        position = int(columns[1]) - 1
    except (UnicodeDecodeError, ValueError) as error:
        raise HybridSelectError(f"invalid VCF record: {stripped!r}") from error
    if position < 0:
        raise HybridSelectError(f"invalid VCF position: {stripped!r}")

    record_end = position + len(columns[3])
    for item in columns[7].split(b";"):
        if item.startswith(b"END="):
            try:
                record_end = int(item[4:])
            except ValueError as error:
                raise HybridSelectError(
                    f"invalid INFO/END in VCF record: {stripped!r}"
                ) from error
    if record_end <= position:
        raise HybridSelectError(f"invalid VCF record end: {stripped!r}")
    return columns, contig, position, record_end


def init_worker(
    vcf_path: str,
    fai_lengths: dict[str, int],
    temp_dir: str,
    slop_size: int,
) -> None:
    global _worker_vcf
    global _worker_index
    global _worker_fai_lengths
    global _worker_temp_dir
    global _worker_slop_size

    _worker_vcf = bgzf.open(vcf_path, "rb")
    _worker_index = tabix.Tabix(vcf_path, "r")
    _worker_fai_lengths = fai_lengths
    _worker_temp_dir = temp_dir
    _worker_slop_size = slop_size


def close_worker() -> None:
    global _worker_vcf
    if _worker_vcf is not None:
        _worker_vcf.close()
    _worker_vcf = None


def iter_shard_records(
    contig: str, start: int, end: int
) -> Iterator[tuple[list[bytes], int, int]]:
    if _worker_vcf is None or _worker_index is None:
        raise HybridSelectError("selection worker was not initialized")

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
            if position >= start:
                yield columns, position, record_end
        if sequential:
            return


def select_shard(
    shard: tuple[int, str, int, int],
) -> tuple[int, str, int, int]:
    shard_number, contig, start, end = shard
    if (
        _worker_fai_lengths is None
        or _worker_temp_dir is None
        or _worker_slop_size is None
    ):
        raise HybridSelectError("selection worker was not initialized")
    if contig not in _worker_fai_lengths:
        raise HybridSelectError(
            f"contig {contig} is absent from the reference index"
        )

    output_path = os.path.join(
        _worker_temp_dir, f"fragment.{shard_number:08d}.bed"
    )
    record_count = 0
    selected_count = 0
    contig_length = _worker_fai_lengths[contig]
    with open(output_path, "wb") as output:
        for columns, position, record_end in iter_shard_records(
            contig, start, end
        ):
            record_count += 1
            if not record_passes(columns):
                continue
            bed_start = max(0, position - _worker_slop_size)
            bed_end = min(contig_length, record_end + _worker_slop_size)
            output.write(f"{contig}\t{bed_start}\t{bed_end}\n".encode())
            selected_count += 1
    return shard_number, output_path, record_count, selected_count


def generate_fragments(
    shards: Sequence[tuple[int, str, int, int]],
    workers: int,
    vcf_path: pathlib.Path,
    fai_lengths: dict[str, int],
    temp_dir: pathlib.Path,
    slop_size: int,
) -> tuple[list[pathlib.Path], int, int]:
    initializer_arguments = (
        str(vcf_path),
        fai_lengths,
        str(temp_dir),
        slop_size,
    )
    results: Iterable[tuple[int, str, int, int]]
    if workers == 1:
        init_worker(*initializer_arguments)
        try:
            results = [select_shard(shard) for shard in shards]
        finally:
            close_worker()
    else:
        context = mp.get_context()
        with context.Pool(
            processes=workers,
            initializer=init_worker,
            initargs=initializer_arguments,
        ) as pool:
            results = list(pool.imap(select_shard, shards, chunksize=1))

    ordered = sorted(results)
    if [item[0] for item in ordered] != list(range(len(shards))):
        raise HybridSelectError(
            "selection shards were not returned exactly once"
        )
    return (
        [pathlib.Path(item[1]) for item in ordered],
        sum(item[2] for item in ordered),
        sum(item[3] for item in ordered),
    )


def publish_output(
    output_path: pathlib.Path, fragments: Sequence[pathlib.Path]
) -> None:
    if not output_path.parent.is_dir():
        raise HybridSelectError(
            f"output directory {output_path.parent} does not exist"
        )
    descriptor, partial_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.partial.",
        suffix=".bed",
        dir=output_path.parent,
    )
    os.close(descriptor)
    partial_path = pathlib.Path(partial_name)
    try:
        with partial_path.open("wb") as output:
            for fragment in fragments:
                with fragment.open("rb") as input_fragment:
                    shutil.copyfileobj(
                        input_fragment, output, length=16 * 1024 * 1024
                    )
        os.replace(partial_path, output_path)
    finally:
        if partial_path.exists():
            partial_path.unlink()


def main(args: argparse.Namespace) -> int:
    if args.threads < 1:
        raise HybridSelectError("--threads must be at least 1")
    if args.step_size < 1:
        raise HybridSelectError("--step-size must be at least 1")
    if args.slop_size < 0:
        raise HybridSelectError("--slop-size must not be negative")

    input_path = pathlib.Path(args.vcf).resolve()
    output_path = pathlib.Path(args.output).resolve()
    reference_index = pathlib.Path(args.reference_fai).resolve()
    contigs = load_vcf_contigs(input_path)
    fai_lengths = load_fai(reference_index)
    for contig, _length in contigs:
        if contig not in fai_lengths:
            raise HybridSelectError(
                f"contig {contig} is absent from the reference index"
            )

    workers = min(args.threads, MAX_SELECT_WORKERS)
    shards = cut_shards(contigs, args.step_size)
    temp_dir = pathlib.Path(
        tempfile.mkdtemp(
            prefix="hybrid-select.",
            dir=os.getenv("SENTIEON_TMPDIR"),
        )
    )
    success = False
    try:
        print(
            "hybrid_select: "
            f"workers={workers} step_size={args.step_size} "
            f"shards={len(shards)} slop_size={args.slop_size} "
            f"tmpdir={temp_dir}",
            file=sys.stderr,
        )
        fragments, record_count, selected_count = generate_fragments(
            shards,
            workers,
            input_path,
            fai_lengths,
            temp_dir,
            args.slop_size,
        )
        publish_output(output_path, fragments)
        print(
            f"hybrid_select: records={record_count} selected={selected_count}",
            file=sys.stderr,
        )
        success = True
    finally:
        if success:
            shutil.rmtree(temp_dir)
        else:
            print(
                f"hybrid_select: preserved failed work directory {temp_dir}",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="sentieon pyexec hybrid_select.py",
        usage="%(prog)s [options] -v VCF --reference-fai REF.fai output.bed",
    )
    parser.add_argument("output", help="output BED file name")
    parser.add_argument(
        "-v", "--vcf", required=True, help="input VCF file name"
    )
    parser.add_argument(
        "--reference-fai", required=True, help="reference FASTA index"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=mp.cpu_count(),
        help="number of worker processes",
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=DEFAULT_STEP_SIZE,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--slop-size",
        type=int,
        default=DEFAULT_SLOP_SIZE,
        help=argparse.SUPPRESS,
    )
    try:
        sys.exit(main(parser.parse_args()))
    except HybridSelectError as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
