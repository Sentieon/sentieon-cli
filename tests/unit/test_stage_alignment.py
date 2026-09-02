"""
Unit tests for the alignment and read-extraction stages
"""

import logging
import pathlib
from unittest.mock import patch

import pytest

from sentieon_cli.dag import DAG
from sentieon_cli.stages.alignment import (
    BwaExtractStage,
    BwaFastqStage,
    BwaRealignStage,
    Minimap2FastqStage,
    Minimap2RealignStage,
    ReadWriterStage,
    aln_suffix,
    find_unzip,
)
from sentieon_cli.stages.base import StageContext, rm_job

RG_LINES = ["@RG\tID:rg1\tSM:sample"]


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


@pytest.fixture
def rg_lines():
    """Stub out the `samtools view -H` call that reads the @RG lines"""
    with patch(
        "sentieon_cli.command_strings.get_rg_lines", return_value=RG_LINES
    ) as mock:
        yield mock


class TestHelpers:
    """The shared helpers the stages and pipelines use"""

    def test_aln_suffix(self):
        assert aln_suffix(True) == "bam"
        assert aln_suffix(False) == "cram"

    def test_find_unzip_prefers_igzip(self):
        with patch(
            "sentieon_cli.stages.alignment.shutil.which",
            return_value="/usr/bin/igzip",
        ):
            assert find_unzip(logging.getLogger("test-unzip")) == "igzip"

    def test_find_unzip_falls_back_to_gzip(self, caplog):
        logger = logging.getLogger("test-unzip-fallback")
        with patch(
            "sentieon_cli.stages.alignment.shutil.which", return_value=None
        ):
            with caplog.at_level(logging.WARNING, logger=logger.name):
                assert find_unzip(logger) == "gzip"
        assert "igzip is recommended" in caplog.text
        assert caplog.records[0].levelno == logging.WARNING

    def test_find_unzip_takes_the_log_level(self, caplog):
        logger = logging.getLogger("test-unzip-level")
        with patch(
            "sentieon_cli.stages.alignment.shutil.which", return_value=None
        ):
            with caplog.at_level(logging.INFO, logger=logger.name):
                assert find_unzip(logger, logging.INFO) == "gzip"
        assert caplog.records[0].levelno == logging.INFO


class TestBwaRealignStage:
    """Re-aligning BAM/CRAM input with bwa"""

    def make(self, tmp_path, **kwargs) -> BwaRealignStage:
        defaults = dict(
            ctx=make_ctx(tmp_path),
            inputs=[tmp_path / "in0.bam", tmp_path / "in1.bam"],
            model_bundle=tmp_path / "model.bundle",
        )
        defaults.update(kwargs)
        return BwaRealignStage(**defaults)  # type: ignore[arg-type]

    def test_one_job_per_input(self, tmp_path, rg_lines):
        result = self.make(tmp_path).build()

        assert [job.name for job in result.jobs] == [
            "bam-align-0",
            "bam-align-1",
        ]
        for job in result.jobs:
            assert job.threads == 4
            assert job.task_name == "alignment"
        assert result.terminal == set(result.jobs)

    def test_header_files_are_written(self, tmp_path, rg_lines):
        self.make(tmp_path).build()

        for i in range(2):
            header = tmp_path / f"input_{i}.hdr"
            assert header.read_text() == RG_LINES[0] + "\n"
        assert rg_lines.call_count == 2

    def test_duplicate_marking_writes_bam_to_the_tmp_dir(
        self, tmp_path, rg_lines
    ):
        result = self.make(tmp_path).build()

        assert result.outputs == [
            tmp_path / "bwa_sorted_0.bam",
            tmp_path / "bwa_sorted_1.bam",
        ]
        assert "--bam_compression 1" in str(result.jobs[0].shell)

    def test_without_duplicate_marking_writes_cram_output(
        self, tmp_path, rg_lines
    ):
        result = self.make(tmp_path, duplicate_marking="none").build()

        assert result.outputs == [
            tmp_path / "output_bwa_sorted_0.cram",
            tmp_path / "output_bwa_sorted_1.cram",
        ]
        assert "--bam_compression 1" not in str(result.jobs[0].shell)

    def test_bam_format_without_duplicate_marking(self, tmp_path, rg_lines):
        result = self.make(
            tmp_path, duplicate_marking="none", bam_format=True
        ).build()

        assert result.outputs == [
            tmp_path / "output_bwa_sorted_0.bam",
            tmp_path / "output_bwa_sorted_1.bam",
        ]

    def test_bam_cleanup_paths(self, tmp_path, rg_lines):
        result = self.make(tmp_path).build()

        assert result.cleanup_paths == [
            tmp_path / "bwa_sorted_0.bam",
            pathlib.Path(str(tmp_path / "bwa_sorted_0.bam") + ".bai"),
            tmp_path / "bwa_sorted_1.bam",
            pathlib.Path(str(tmp_path / "bwa_sorted_1.bam") + ".bai"),
        ]

    def test_cram_cleanup_paths_add_the_crai(self, tmp_path, rg_lines):
        result = self.make(tmp_path, duplicate_marking="none").build()

        out = str(tmp_path / "output_bwa_sorted_0.cram")
        assert result.cleanup_paths[:3] == [
            pathlib.Path(out),
            pathlib.Path(out + ".bai"),
            pathlib.Path(out + ".crai"),
        ]

    def test_collate_and_input_ref(self, tmp_path, rg_lines):
        result = self.make(
            tmp_path,
            inputs=[tmp_path / "in0.cram"],
            collate=True,
            input_ref=tmp_path / "decode.fa",
        ).build()

        shell = str(result.jobs[0].shell)
        assert "samtools collate" in shell
        assert f"--reference {tmp_path}/decode.fa" in shell

    def test_build_inserts_nothing(self, tmp_path, rg_lines):
        dag = DAG()
        self.make(tmp_path).build()

        assert not dag.ready_jobs and not dag.waiting_jobs


class TestBwaFastqStage:
    """Aligning short-read fastq with bwa"""

    def make(self, tmp_path, **kwargs) -> BwaFastqStage:
        defaults = dict(
            ctx=make_ctx(tmp_path),
            r1_fastq=[tmp_path / "r1.fq.gz"],
            r2_fastq=[tmp_path / "r2.fq.gz"],
            readgroups=[r"@RG\tID:rg1\tSM:sample"],
            model_bundle=tmp_path / "model.bundle",
        )
        defaults.update(kwargs)
        return BwaFastqStage(**defaults)  # type: ignore[arg-type]

    def test_without_numa_one_job_takes_every_core(self, tmp_path):
        result = self.make(tmp_path).build()

        assert [job.name for job in result.jobs] == ["bam-align-0-0"]
        job = result.jobs[0]
        assert job.threads == 4
        assert job.resources == {"node0": 1}
        assert job.task_name == "alignment"
        assert "taskset" not in str(job.shell)
        assert "--split" not in str(job.shell)

    def test_numa_splits_the_pair_across_the_nodes(self, tmp_path):
        result = self.make(
            tmp_path, ctx=make_ctx(tmp_path, cores=8), numa_nodes=["0-3", "4-7"]
        ).build()

        assert [job.name for job in result.jobs] == [
            "bam-align-0-0",
            "bam-align-0-1",
        ]
        assert [job.threads for job in result.jobs] == [4, 4]
        assert [job.resources for job in result.jobs] == [
            {"node0": 1},
            {"node1": 1},
        ]
        first, second = (str(job.shell) for job in result.jobs)
        assert "taskset -c 0-3" in first and "0/2" in first
        assert "taskset -c 4-7" in second and "1/2" in second

    def test_one_job_per_pair_and_node(self, tmp_path):
        result = self.make(
            tmp_path,
            r1_fastq=[tmp_path / "a1.fq.gz", tmp_path / "b1.fq.gz"],
            r2_fastq=[tmp_path / "a2.fq.gz", tmp_path / "b2.fq.gz"],
            readgroups=[r"@RG\tID:a", r"@RG\tID:b"],
            numa_nodes=["0-3", "4-7"],
        ).build()

        assert [job.name for job in result.jobs] == [
            "bam-align-0-0",
            "bam-align-0-1",
            "bam-align-1-0",
            "bam-align-1-1",
        ]

    def test_no_fastq_builds_nothing(self, tmp_path):
        result = self.make(tmp_path, r1_fastq=[], r2_fastq=[], readgroups=[])

        built = result.build()
        assert built.jobs == []
        assert built.outputs == []
        assert built.cleanup_paths == []

    def test_unzip_reaches_the_command(self, tmp_path):
        result = self.make(tmp_path, unzip="igzip").build()

        assert "igzip -dc" in str(result.jobs[0].shell)

    def test_outputs_and_cleanup_paths(self, tmp_path):
        result = self.make(tmp_path, duplicate_marking="none").build()

        out = tmp_path / "output_bwa_sorted_fq_0_0.cram"
        assert result.outputs == [out]
        assert result.cleanup_paths == [
            out,
            pathlib.Path(str(out) + ".bai"),
            pathlib.Path(str(out) + ".crai"),
        ]


class TestMinimap2RealignStage:
    """Re-aligning long-read input with minimap2"""

    def make(self, tmp_path, **kwargs) -> Minimap2RealignStage:
        defaults = dict(
            ctx=make_ctx(tmp_path),
            inputs=[tmp_path / "lr.bam"],
            model_bundle=tmp_path / "model.bundle",
            sample_name="sample",
        )
        defaults.update(kwargs)
        return Minimap2RealignStage(**defaults)  # type: ignore[arg-type]

    def test_job_metadata_and_outputs(self, tmp_path, rg_lines):
        result = self.make(tmp_path).build()

        assert [job.name for job in result.jobs] == ["bam-realign-0"]
        assert result.jobs[0].threads == 4
        assert result.jobs[0].task_name == "alignment"
        assert result.outputs == [tmp_path / "output_mm2_sorted_0.cram"]
        assert result.cleanup_paths == []

    def test_default_model_comes_from_the_bundle(self, tmp_path, rg_lines):
        result = self.make(tmp_path).build()

        assert f"-x {tmp_path}/model.bundle/minimap2.model" in str(
            result.jobs[0].shell
        )

    def test_minimap2_model_overrides_the_bundle_default(
        self, tmp_path, rg_lines
    ):
        model = tmp_path / "model.bundle" / "minimap2_lr.model"
        result = self.make(tmp_path, minimap2_model=model).build()

        assert f"-x {model}" in str(result.jobs[0].shell)

    def test_bam_format_and_input_ref(self, tmp_path, rg_lines):
        result = self.make(
            tmp_path, bam_format=True, input_ref=tmp_path / "lr_ref.fa"
        ).build()

        assert result.outputs == [tmp_path / "output_mm2_sorted_0.bam"]
        assert f"--reference {tmp_path}/lr_ref.fa" in str(result.jobs[0].shell)


class TestMinimap2FastqStage:
    """Aligning long-read fastq with minimap2"""

    def test_one_job_per_fastq(self, tmp_path):
        result = Minimap2FastqStage(
            ctx=make_ctx(tmp_path),
            fastq=[tmp_path / "a.fq.gz", tmp_path / "b.fq.gz"],
            readgroups=[r"@RG\tID:a", r"@RG\tID:b"],
            model_bundle=tmp_path / "model.bundle",
            unzip="igzip",
        ).build()

        assert [job.name for job in result.jobs] == ["align-0", "align-1"]
        assert result.outputs == [
            tmp_path / "output_mm2_sorted_fq_0.cram",
            tmp_path / "output_mm2_sorted_fq_1.cram",
        ]
        assert "igzip -dc" in str(result.jobs[0].shell)
        assert result.jobs[0].task_name == "alignment"

    def test_no_fastq_builds_nothing(self, tmp_path):
        result = Minimap2FastqStage(
            ctx=make_ctx(tmp_path),
            fastq=[],
            readgroups=[],
            model_bundle=tmp_path / "model.bundle",
        ).build()

        assert result.jobs == []
        assert result.outputs == []


class TestBwaExtractStage:
    """The pangenome bwa alignment and read-extraction pass"""

    def make(self, tmp_path, **kwargs) -> BwaExtractStage:
        defaults = dict(
            ctx=make_ctx(tmp_path),
            output_bam=tmp_path / "sample-bwa.bam",
            output_fastq=tmp_path / "sample-extract.fq.gz",
            r1_fastq=[tmp_path / "r1.fq.gz"],
            r2_fastq=[tmp_path / "r2.fq.gz"],
            readgroup=r"@RG\tID:rg1-bwa\tSM:sample",
            extract_model=tmp_path / "model.bundle" / "extract.model",
            bwa_model=tmp_path / "model.bundle" / "bwa.model",
        )
        defaults.update(kwargs)
        return BwaExtractStage(**defaults)  # type: ignore[arg-type]

    def test_single_job(self, tmp_path):
        result = self.make(tmp_path).build()

        assert len(result.jobs) == 1
        job = result.jobs[0]
        assert job.name == "bwa-extract"
        assert job.threads == 4
        assert job.task_name == "alignment"
        assert result.outputs == [
            tmp_path / "sample-bwa.bam",
            tmp_path / "sample-extract.fq.gz",
        ]
        assert result.cleanup_paths == []

    def test_readgroup_is_used_verbatim(self, tmp_path):
        result = self.make(
            tmp_path, readgroup=r"@RG\tID:rg1\tSM:sample\tLR:0"
        ).build()

        assert r"@RG\tID:rg1\tSM:sample\tLR:0" in str(result.jobs[0].shell)

    def test_unzip_reaches_the_command(self, tmp_path):
        result = self.make(tmp_path, unzip="igzip").build()

        assert "igzip -dc" in str(result.jobs[0].shell)


class TestReadWriterStage:
    """`sentieon driver --algo ReadWriter` merges and extracts"""

    def test_merge(self, tmp_path):
        result = ReadWriterStage(
            ctx=make_ctx(tmp_path),
            inputs=[tmp_path / "one.bam", tmp_path / "two.bam"],
            output=tmp_path / "merged.bam",
            name="merge-bam",
            task_name="alignment",
            threads=0,
        ).build()

        job = result.jobs[0]
        assert job.name == "merge-bam"
        assert job.threads == 0
        assert job.task_name == "alignment"
        assert str(job.shell) == (
            f"sentieon driver --input {tmp_path}/one.bam "
            f"--input {tmp_path}/two.bam --reference {tmp_path}/ref.fa "
            f"--thread_count 4 --algo ReadWriter {tmp_path}/merged.bam"
        )
        assert result.outputs == [tmp_path / "merged.bam"]

    def test_string_interval_extracts_a_locus(self, tmp_path):
        result = ReadWriterStage(
            ctx=make_ctx(tmp_path),
            inputs=[tmp_path / "sample.bam"],
            output=tmp_path / "sample_hla.bam",
            name="t1k-hla-extract",
            task_name="t1k",
            interval="chr6:28510020-33480577",
        ).build()

        job = result.jobs[0]
        assert job.name == "t1k-hla-extract"
        assert job.threads == 4  # defaults to the run's core count
        assert job.task_name == "t1k"
        assert "--interval chr6:28510020-33480577" in str(job.shell)


class TestAddTo:
    """`add_to` inserts what `build` returned"""

    def test_jobs_are_roots_without_upstream(self, tmp_path):
        dag = DAG()
        result = BwaFastqStage(
            ctx=make_ctx(tmp_path, cores=8),
            r1_fastq=[tmp_path / "r1.fq.gz"],
            r2_fastq=[tmp_path / "r2.fq.gz"],
            readgroups=[r"@RG\tID:rg1"],
            model_bundle=tmp_path / "model.bundle",
            numa_nodes=["0-3", "4-7"],
        ).add_to(dag)

        assert len(result.jobs) == 2
        for job in result.jobs:
            assert job in dag.ready_jobs
        assert not dag.waiting_jobs

    def test_every_job_waits_on_the_upstream_jobs(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = Minimap2FastqStage(
            ctx=make_ctx(tmp_path),
            fastq=[tmp_path / "a.fq.gz", tmp_path / "b.fq.gz"],
            readgroups=[r"@RG\tID:a", r"@RG\tID:b"],
            model_bundle=tmp_path / "model.bundle",
        ).add_to(dag, [upstream])

        for job in result.jobs:
            assert dag.waiting_jobs[job] == {upstream}
