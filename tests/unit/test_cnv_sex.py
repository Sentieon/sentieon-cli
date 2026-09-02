"""
Unit tests for sex-aware CNV calling: reference build detection, PAR BED
selection, the CNVscope algo arguments, and the `--sample_sex` argument.
"""

import argparse
import logging
import multiprocessing as mp
import pathlib
from typing import List

import pytest

from sentieon_cli import shard
from sentieon_cli.driver import CNVscope
from sentieon_cli.pipeline import BasePipeline
from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext
from sentieon_cli.stages.cnv import CNVscopeStage
from sentieon_cli.shard import (
    BUILD_SIGNATURES,
    detect_reference_build,
    par_bed_for_build,
    ploidy_contigs_for_build,
)
from sentieon_cli.util import (
    SampleSex,
    cnvscope_sex_args,
    sample_sex_arg,
)

PACKAGE_LOGGER = "sentieon_cli"


class _RecordingHandler(logging.Handler):
    """Collects formatted messages from the package logger."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


@pytest.fixture
def messages():
    """Capture the package logger, which does not propagate to the root."""
    handler = _RecordingHandler()
    package_logger = logging.getLogger(PACKAGE_LOGGER)
    package_logger.addHandler(handler)
    yield handler.messages
    package_logger.removeHandler(handler)


def _fai(contigs) -> dict:
    """A parsed fasta index with the supplied contig lengths"""
    return {
        ctg: {"length": length, "offset": 0, "linebases": 60, "linewidth": 61}
        for ctg, length in contigs.items()
    }


class _MinimalPipeline(BasePipeline):
    """A concrete pipeline for exercising the shared sex machinery"""

    def validate(self) -> None:
        pass

    def configure(self) -> None:
        pass

    def build_dag(self) -> DAG:
        return DAG()


def _pipeline(**kwargs) -> _MinimalPipeline:
    pipeline = _MinimalPipeline()
    pipeline.setup_logging(argparse.Namespace(loglevel="DEBUG"))
    for key, value in kwargs.items():
        setattr(pipeline, key, value)
    return pipeline


# Reference build detection


@pytest.mark.parametrize("build", sorted(BUILD_SIGNATURES))
def test_detect_reference_build_identifies_each_build(build):
    assert detect_reference_build(_fai(BUILD_SIGNATURES[build])) == build


def test_detect_reference_build_distinguishes_hg19_from_b37():
    # hg19 and b37 share contig lengths and differ only by the 'chr' prefix
    hg19 = detect_reference_build(_fai(BUILD_SIGNATURES["hg19"]))
    b37 = detect_reference_build(_fai(BUILD_SIGNATURES["b37"]))
    assert hg19 == "hg19"
    assert b37 == "b37"


def test_detect_reference_build_requires_all_signature_contigs():
    signature = dict(BUILD_SIGNATURES["hg38"])
    del signature["chrX"]
    assert detect_reference_build(_fai(signature)) is None


def test_detect_reference_build_requires_exact_lengths():
    signature = dict(BUILD_SIGNATURES["hg38"])
    signature["chr2"] = signature["chr2"] - 1
    assert detect_reference_build(_fai(signature)) is None


def test_detect_reference_build_with_an_empty_index():
    assert detect_reference_build({}) is None


def test_detect_reference_build_with_an_unknown_reference():
    assert detect_reference_build(_fai({"contig1": 1000})) is None


# PAR BED selection


def test_par_bed_for_build_without_a_build():
    assert par_bed_for_build(None) is None


def test_par_bed_for_build_with_an_absent_file(messages):
    assert par_bed_for_build("not_a_build") is None
    assert any("PAR" in msg for msg in messages)


def test_par_bed_for_build_treats_an_empty_file_as_absent(
    tmp_path, monkeypatch, messages
):
    # The packaged BED files ship empty until their contents are added
    empty_bed = tmp_path.joinpath("par_hg38.bed")
    empty_bed.touch()
    monkeypatch.setattr(shard, "files", lambda pkg: tmp_path)

    assert par_bed_for_build("hg38") is None
    assert any("PAR" in msg for msg in messages)


def test_par_bed_for_build_returns_a_packaged_file(tmp_path, monkeypatch):
    par_bed = tmp_path.joinpath("par_hg38.bed")
    par_bed.write_text("chrX\t10000\t2781479\n")
    monkeypatch.setattr(shard, "files", lambda pkg: tmp_path)

    assert par_bed_for_build("hg38") == par_bed


# Ploidy estimation contigs


def test_ploidy_contigs_for_b37_use_bare_contig_names():
    contigs = ploidy_contigs_for_build("b37")
    assert contigs.contigs is not None
    assert contigs.contigs[0] == "1"
    assert contigs.autosomes == [str(i) for i in range(1, 23)]
    assert contigs.x_contig == "X"
    assert contigs.y_contig == "Y"


@pytest.mark.parametrize("build", ["hg38", "hg19", "chm13", None])
def test_ploidy_contigs_default_to_the_script_defaults(build):
    contigs = ploidy_contigs_for_build(build)
    assert contigs.contigs is None
    assert contigs.autosomes is None
    assert contigs.x_contig is None
    assert contigs.y_contig is None


# The `--sample_sex` argument


@pytest.mark.parametrize("value", ["male", "MALE", " Male "])
def test_sample_sex_arg_parses_male(value):
    assert sample_sex_arg(value) is SampleSex.MALE


@pytest.mark.parametrize("value", ["female", "Female", "FEMALE"])
def test_sample_sex_arg_parses_female(value):
    assert sample_sex_arg(value) is SampleSex.FEMALE


@pytest.mark.parametrize("value", ["m", "f", "unknown", ""])
def test_sample_sex_arg_rejects_other_values(value):
    with pytest.raises(argparse.ArgumentTypeError):
        sample_sex_arg(value)


# The CNVscope `--sex` and `--par` arguments


def test_cnvscope_sex_args_for_a_male_sample():
    par_bed = pathlib.Path("/par.bed")
    assert cnvscope_sex_args(SampleSex.MALE, par_bed) == ("M", par_bed)


def test_cnvscope_sex_args_for_a_female_sample():
    # A female sample does not need the PAR regions
    par_bed = pathlib.Path("/par.bed")
    assert cnvscope_sex_args(SampleSex.FEMALE, par_bed) == ("F", None)


@pytest.mark.parametrize("sex", [None, SampleSex.UNKNOWN])
def test_cnvscope_sex_args_without_a_known_sex(sex, messages):
    assert cnvscope_sex_args(sex, None) == (None, None)
    assert any("diploid" in msg for msg in messages)


def test_cnvscope_algo_emits_the_sex_arguments():
    cmd = CNVscope(
        pathlib.Path("out.vcf.gz"),
        pathlib.Path("bundle/cnv.model"),
        sex="M",
        par=pathlib.Path("par.bed"),
    ).build_cmd()

    assert cmd == [
        "--algo",
        "CNVscope",
        "--model",
        "bundle/cnv.model",
        "--sex",
        "M",
        "--par",
        "par.bed",
        "out.vcf.gz",
    ]


def test_cnvscope_algo_omits_the_sex_arguments_by_default():
    cmd = CNVscope(
        pathlib.Path("out.vcf.gz"),
        pathlib.Path("bundle/cnv.model"),
    ).build_cmd()

    assert "--sex" not in cmd
    assert "--par" not in cmd
    assert cmd[-1] == "out.vcf.gz"


# The shared fail-fast helpers


@pytest.mark.parametrize(
    "sample_sex", [None, SampleSex.MALE, SampleSex.FEMALE, SampleSex.UNKNOWN]
)
def test_validate_cnv_par_requires_a_par_bed_for_every_sex(
    sample_sex, messages
):
    # CNV calling needs the PAR BED file whatever the sample sex is, and
    # the run stops during validation without one
    pipeline = _pipeline(sample_sex=sample_sex, cnv_par_bed=None)

    with pytest.raises(SystemExit) as excinfo:
        pipeline.validate_cnv_par(True)
    assert excinfo.value.code == 2
    assert any("--par_bed" in msg for msg in messages)


@pytest.mark.parametrize(
    "sample_sex", [None, SampleSex.MALE, SampleSex.FEMALE, SampleSex.UNKNOWN]
)
def test_validate_cnv_par_accepts_a_par_bed(sample_sex):
    pipeline = _pipeline(
        sample_sex=sample_sex, cnv_par_bed=pathlib.Path("/par.bed")
    )

    pipeline.validate_cnv_par(True)  # no SystemExit


def test_validate_cnv_par_is_skipped_when_cnvs_are_not_called():
    # A PAR BED file is only needed by CNV calling
    pipeline = _pipeline(sample_sex=SampleSex.MALE, cnv_par_bed=None)

    pipeline.validate_cnv_par(False)  # no SystemExit


def test_validate_cnv_par_errors_in_a_dry_run_too():
    # There is no dry-run exception: a dry-run of a real command needs
    # the same arguments the real run would use
    pipeline = _pipeline(cnv_par_bed=None, dry_run=True)

    with pytest.raises(SystemExit) as excinfo:
        pipeline.validate_cnv_par(True)
    assert excinfo.value.code == 2


def test_cnv_sex_args_for_a_male_sample_with_a_par_bed():
    par_bed = pathlib.Path("/par.bed")
    pipeline = _pipeline(sample_sex=SampleSex.MALE, cnv_par_bed=par_bed)

    assert pipeline.cnv_sex_args() == ("M", par_bed)


def test_resolve_cnv_par_bed_prefers_the_supplied_file(tmp_path):
    par_bed = tmp_path.joinpath("custom_par.bed")
    par_bed.write_text("chrX\t10000\t2781479\n")
    pipeline = _pipeline()

    pipeline.resolve_cnv_par_bed(_fai(BUILD_SIGNATURES["hg38"]), par_bed)

    assert pipeline.cnv_par_bed == par_bed
    # The build is still identified, for the ploidy estimation contigs
    assert pipeline.reference_build == "hg38"


def test_resolve_cnv_par_bed_falls_back_to_the_packaged_file(
    tmp_path, monkeypatch
):
    packaged = tmp_path.joinpath("par_b37.bed")
    packaged.write_text("X\t60001\t2699520\n")
    monkeypatch.setattr(shard, "files", lambda pkg: tmp_path)
    pipeline = _pipeline()

    pipeline.resolve_cnv_par_bed(_fai(BUILD_SIGNATURES["b37"]))

    assert pipeline.reference_build == "b37"
    assert pipeline.cnv_par_bed == packaged


def test_resolve_cnv_par_bed_skips_the_lookup_without_cnv_calling(messages):
    # The build is still identified, for the ploidy estimation contigs,
    # but a run without CNV calling needs no PAR BED file
    pipeline = _pipeline()

    pipeline.resolve_cnv_par_bed(
        _fai(BUILD_SIGNATURES["hg38"]), None, cnv_will_run=False
    )

    assert pipeline.reference_build == "hg38"
    assert pipeline.cnv_par_bed is None
    assert not any("PAR" in msg for msg in messages)


# `CNVscopeStage` passes the sample sex through to CNVscope


def _cnvscope_cmd(**kwargs) -> str:
    tmp_dir = pathlib.Path("/tmp/cnv")
    output_vcf = pathlib.Path("/out/sample.vcf.gz")
    ctx = StageContext(
        reference=pathlib.Path("/ref/reference.fa"),
        output_vcf=output_vcf,
        tmp_dir=tmp_dir,
        cores=mp.cpu_count(),
        dry_run=True,
        skip_version_check=True,
    )
    result = CNVscopeStage(
        ctx=ctx,
        inputs=[pathlib.Path("/in/sample.cram")],
        model=pathlib.Path("/bundle/model.bundle/cnv.model"),
        cnvscope_vcf=tmp_dir.joinpath("cnvscope.vcf.gz"),
        cnv_vcf=pathlib.Path(
            str(output_vcf).replace(".vcf.gz", ".cnv.vcf.gz")
        ),
        **kwargs,
    ).add_to(DAG())
    return str(result.cnvscope_job.shell)


def test_call_cnvs_passes_the_sex_of_a_male_sample():
    cmd = _cnvscope_cmd(
        sample_sex=SampleSex.MALE, par_bed=pathlib.Path("/par.bed")
    )
    assert "--sex M" in cmd
    assert "--par /par.bed" in cmd


def test_call_cnvs_passes_the_sex_of_a_female_sample():
    cmd = _cnvscope_cmd(
        sample_sex=SampleSex.FEMALE, par_bed=pathlib.Path("/par.bed")
    )
    assert "--sex F" in cmd
    assert "--par" not in cmd


def test_call_cnvs_without_a_sample_sex():
    cmd = _cnvscope_cmd()
    assert "--sex" not in cmd
    assert "--par" not in cmd
