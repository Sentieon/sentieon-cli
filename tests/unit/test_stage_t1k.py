"""
Unit tests for the T1K HLA/KIR genotyping stage
"""

import pathlib

from sentieon_cli.dag import DAG
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.t1k import (
    DEFAULT_T1K_HLA_LOCUS,
    DEFAULT_T1K_KIR_LOCUS,
    T1K_MIN_VERSIONS,
    T1KStage,
)


def make_ctx(tmp_path: pathlib.Path, cores: int = 4) -> StageContext:
    """A StageContext over a temporary directory"""
    return StageContext(
        reference=tmp_path / "ref.fa",
        output_vcf=tmp_path / "output.vcf.gz",
        tmp_dir=tmp_path,
        cores=cores,
        dry_run=True,
        skip_version_check=True,
    )


def make_stage(tmp_path: pathlib.Path, **kwargs) -> T1KStage:
    """A T1KStage genotyping HLA, unless the caller says otherwise"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        inputs=[tmp_path / "sample.cram"],
        gene_seq=tmp_path / "hla_seq.fa",
        gene_coord=tmp_path / "hla_coord.fa",
        locus=DEFAULT_T1K_HLA_LOCUS,
        preset="hla-wgs",
        tag="hla",
        output_dir=tmp_path / "output_hla",
    )
    defaults.update(kwargs)
    return T1KStage(**defaults)  # type: ignore[arg-type]


def make_kir_stage(tmp_path: pathlib.Path, **kwargs) -> T1KStage:
    """A T1KStage genotyping KIR"""
    return make_stage(
        tmp_path,
        gene_seq=tmp_path / "kir_seq.fa",
        gene_coord=tmp_path / "kir_coord.fa",
        locus=DEFAULT_T1K_KIR_LOCUS,
        preset="kir-wgs",
        tag="kir",
        output_dir=tmp_path / "output_kir",
        **kwargs,
    )


class TestConstants:
    """The version constant and the default loci live with the stage"""

    def test_no_minimum_version(self):
        assert T1K_MIN_VERSIONS == {"run-t1k": None}

    def test_default_loci(self):
        assert DEFAULT_T1K_HLA_LOCUS == "chr6:28510020-33480577"
        assert DEFAULT_T1K_KIR_LOCUS == "chr19:53100000-55800000"


class TestCommands:
    """The extraction and `run-t1k` commands the stage builds"""

    def test_hla_extract_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert str(result.extract_job.shell) == (
            f"sentieon driver --input {tmp_path}/sample.cram "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            "--interval chr6:28510020-33480577 "
            f"--algo ReadWriter {tmp_path}/sample_hla.bam"
        )

    def test_hla_t1k_command(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert str(result.t1k_job.shell) == (
            "run-t1k --abnormalUnmapFlag -t 4 --preset hla-wgs "
            f"-f {tmp_path}/hla_seq.fa -c {tmp_path}/hla_coord.fa "
            f"--od {tmp_path}/output_hla -b {tmp_path}/sample_hla.bam"
        )

    def test_kir_extract_command(self, tmp_path):
        result = make_kir_stage(tmp_path).add_to(DAG())

        assert str(result.extract_job.shell) == (
            f"sentieon driver --input {tmp_path}/sample.cram "
            f"--reference {tmp_path}/ref.fa --thread_count 4 "
            "--interval chr19:53100000-55800000 "
            f"--algo ReadWriter {tmp_path}/sample_kir.bam"
        )

    def test_kir_t1k_command(self, tmp_path):
        result = make_kir_stage(tmp_path).add_to(DAG())

        assert str(result.t1k_job.shell) == (
            "run-t1k --abnormalUnmapFlag -t 4 --preset kir-wgs "
            f"-f {tmp_path}/kir_seq.fa -c {tmp_path}/kir_coord.fa "
            f"--od {tmp_path}/output_kir -b {tmp_path}/sample_kir.bam"
        )

    def test_multiple_inputs_are_merged_by_the_extract(self, tmp_path):
        result = make_stage(
            tmp_path,
            inputs=[tmp_path / "one.cram", tmp_path / "two.cram"],
        ).add_to(DAG())

        shell = str(result.extract_job.shell)
        assert f"--input {tmp_path}/one.cram" in shell
        assert f"--input {tmp_path}/two.cram" in shell

    def test_a_custom_locus(self, tmp_path):
        result = make_stage(tmp_path, locus="chr6:1-100").add_to(DAG())

        assert "--interval chr6:1-100 " in str(result.extract_job.shell)


class TestJobMetadata:
    """Job names, task name and thread counts"""

    def test_hla_job_names(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.extract_job.name == "t1k-hla-extract"
        assert result.t1k_job.name == "t1k-hla"

    def test_kir_job_names(self, tmp_path):
        result = make_kir_stage(tmp_path).add_to(DAG())

        assert result.extract_job.name == "t1k-kir-extract"
        assert result.t1k_job.name == "t1k-kir"

    def test_task_name_and_threads(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.jobs:
            assert job.task_name == "t1k"
            assert job.threads == 4

    def test_overridden_task_name(self, tmp_path):
        result = make_stage(tmp_path, task_name="hla").add_to(DAG())

        for job in result.jobs:
            assert job.task_name == "hla"


class TestDagWiring:
    """The edges the stage inserts"""

    def test_upstream_lands_on_the_extract_only(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.extract_job] == {upstream}
        assert dag.waiting_jobs[result.t1k_job] == {result.extract_job}

    def test_the_extract_is_a_root_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        assert result.extract_job in dag.ready_jobs
        assert dag.waiting_jobs[result.t1k_job] == {result.extract_job}

    def test_result_fields(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.jobs == [result.extract_job, result.t1k_job]
        assert result.terminal == {result.t1k_job}
        assert result.extracted_bam == tmp_path / "sample_hla.bam"
        assert result.output_dir == tmp_path / "output_hla"
