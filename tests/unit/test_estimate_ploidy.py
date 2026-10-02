"""
Unit tests for the `estimate_ploidy.py` helper script.

The script is not importable as part of the package (the `scripts`
directory has no `__init__.py`), so it is loaded from its file path.
"""

import importlib.util
import io
import json
import pathlib
import sys
from types import SimpleNamespace

import pytest

SCRIPT = (
    pathlib.Path(__file__).resolve().parents[2]
    / "sentieon_cli"
    / "scripts"
    / "estimate_ploidy.py"
)


def _load_script():
    spec = importlib.util.spec_from_file_location("estimate_ploidy", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


estimate_ploidy = _load_script()


def _idxstats(contigs):
    """Format idxstats output for the supplied {contig: (length, reads)}"""
    lines = [
        f"{ctg}\t{length}\t{reads}\t0" for ctg, (length, reads) in contigs
    ]
    lines.append("*\t0\t0\t0")
    return "\n".join(lines) + "\n"


def _chr_contigs(x_reads, y_reads, prefix="chr"):
    """Autosomes at ploidy two plus the requested sex chromosome reads"""
    contigs = [(f"{prefix}{i}", (1000, 100)) for i in range(1, 23)]
    contigs.append((f"{prefix}X", (1000, x_reads)))
    if y_reads is not None:
        contigs.append((f"{prefix}Y", (1000, y_reads)))
    return contigs


def _run(monkeypatch, idxstats_contigs, **kwargs):
    """Run the script's `main` with a mocked `samtools idxstats`"""
    monkeypatch.setattr(
        estimate_ploidy.subprocess,
        "run",
        lambda *a, **kw: SimpleNamespace(stdout=_idxstats(idxstats_contigs)),
    )
    outfile = io.StringIO()
    args = SimpleNamespace(
        input_bam=["sample.bam"],
        contigs=kwargs.get("contigs", estimate_ploidy.DEFAULT_CONTIGS),
        autosomes=kwargs.get("autosomes", estimate_ploidy.DEFAULT_AUTOSOMES),
        x_contig=kwargs.get("x_contig", estimate_ploidy.DEFAULT_X_CONTIG),
        y_contig=kwargs.get("y_contig", estimate_ploidy.DEFAULT_Y_CONTIG),
        outfile=outfile,
    )
    with pytest.raises(SystemExit) as excinfo:
        estimate_ploidy.main(args)
    assert excinfo.value.code == 0
    return json.loads(outfile.getvalue())


def test_female_sample_with_chr_prefixed_contigs(monkeypatch):
    results = _run(monkeypatch, _chr_contigs(x_reads=100, y_reads=0))

    assert results["sex"] == "female"
    assert results["contigs"]["chrX"]["ploidy"] == "2"
    assert results["contigs"]["chrY"]["ploidy"] == "0"


def test_male_sample_with_chr_prefixed_contigs(monkeypatch):
    results = _run(monkeypatch, _chr_contigs(x_reads=50, y_reads=50))

    assert results["sex"] == "male"
    assert results["contigs"]["chrX"]["ploidy"] == "1"
    assert results["contigs"]["chrY"]["ploidy"] == "1"


def test_male_sample_with_bare_contig_names(monkeypatch):
    # A b37-style reference, using the contig name arguments
    results = _run(
        monkeypatch,
        _chr_contigs(x_reads=50, y_reads=50, prefix=""),
        contigs=[str(i) for i in range(1, 23)] + ["X", "Y"],
        autosomes=[str(i) for i in range(1, 23)],
        x_contig="X",
        y_contig="Y",
    )

    assert results["sex"] == "male"
    assert results["contigs"]["X"]["ploidy"] == "1"


def test_a_missing_y_contig_reports_an_unknown_sex(monkeypatch):
    # An absent sex chromosome used to raise an uncaught KeyError
    results = _run(monkeypatch, _chr_contigs(x_reads=100, y_reads=None))

    assert results["sex"] == "Unknown"
    assert "chrY" not in results["contigs"]


def test_no_matching_contigs_reports_an_unknown_sex(monkeypatch):
    # A reference whose contig names do not match the supplied ones
    results = _run(monkeypatch, [("scaffold1", (1000, 100))])

    assert results == {"contigs": {}, "sex": "Unknown"}


def test_the_contig_arguments_are_parsed(monkeypatch):
    argv = [
        "estimate_ploidy.py",
        "-i",
        "sample.bam",
        "--contigs",
        "1",
        "2",
        "--autosomes",
        "1",
        "2",
        "--x_contig",
        "X",
        "--y_contig",
        "Y",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    args = estimate_ploidy.process_args()

    assert args.contigs == ["1", "2"]
    assert args.autosomes == ["1", "2"]
    assert args.x_contig == "X"
    assert args.y_contig == "Y"


def test_the_contig_arguments_default_to_chr_prefixed_names(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["estimate_ploidy.py", "-i", "s.bam"])
    args = estimate_ploidy.process_args()

    assert args.contigs == estimate_ploidy.DEFAULT_CONTIGS
    assert args.x_contig == "chrX"
    assert args.y_contig == "chrY"
