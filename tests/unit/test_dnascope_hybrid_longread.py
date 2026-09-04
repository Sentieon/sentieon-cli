"""
Unit tests for the DNAscope hybrid pipeline's long-read jobs

The hybrid pipeline used to inherit its long-read realignment, mosdepth and
SV calling from `DNAscopeLRPipeline`; it now composes the same stages
itself, with its own argument defaults.
"""

import pathlib
import tempfile
from abc import ABC
from unittest.mock import MagicMock, patch

import pytest

from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.pipeline import BasePipeline

RG_LINES = ["@RG\tID:lr_rg1\tSM:sample1"]


class TestDNAscopeHybridLongRead:
    """The long-read half of the hybrid DAG"""

    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_input_ref = self.mock_dir / "input_reference.fa"
        self.mock_sr_aln = [self.mock_dir / "short.bam"]
        self.mock_lr_aln = [self.mock_dir / "long.bam"]
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_bed = self.mock_dir / "interval.bed"

        for file_path in [
            self.mock_ref,
            self.mock_input_ref,
            self.mock_sr_aln[0],
            self.mock_lr_aln[0],
            self.mock_bundle,
            self.mock_bed,
        ]:
            file_path.touch()

    def create_pipeline(self):
        """Create a DNAscopeHybridPipeline that realigns its long reads"""
        with patch("sys.exit"):
            pipeline = DNAscopeHybridPipeline()

        pipeline.logger = MagicMock()

        pipeline.output_vcf = self.mock_vcf
        pipeline.reference = self.mock_ref
        pipeline.model_bundle = self.mock_bundle
        pipeline.bed = self.mock_bed
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = self.mock_dir
        pipeline.pop_vcf = None

        # State normally set by validate()
        pipeline.fai_data = {"chr1": {"length": 1000}}
        pipeline.shards = [MagicMock()]
        pipeline.shards[0].contig = "chr1"
        pipeline.shards[0].start = 1
        pipeline.shards[0].stop = 1000
        pipeline.pop_vcf_contigs = {"chr1": 1000}
        pipeline.lr_aln_readgroups = [[{"ID": "lr_rg1", "SM": "sample1"}]]
        pipeline.sr_aln_readgroups = [[{"ID": "sr_rg1", "SM": "sample1"}]]
        pipeline.hybrid_rg_sm = "sample1"
        pipeline.hybrid_set_rg = False
        pipeline.shortread_tech = "Illumina"
        pipeline.longread_tech = "ONT"
        pipeline.sr_aln = self.mock_sr_aln
        pipeline.lr_aln = self.mock_lr_aln

        # Realign the long reads
        pipeline.lr_align_input = True
        pipeline.lr_input_ref = self.mock_input_ref

        return pipeline

    def build_dag(self, pipeline):
        """Build the DAG and index the jobs by name"""
        with patch(
            "sentieon_cli.command_strings.get_rg_lines", return_value=RG_LINES
        ):
            dag = pipeline.build_dag()
        jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        return dag, {job.name: job for job in jobs}

    def dep_names(self, dag, job):
        """The names of a job's dependencies"""
        return {dep.name for dep in dag.waiting_jobs.get(job, set())}

    def test_minimap2_args_default(self):
        """The hybrid pipeline aligns with `-Y`, not the long-read `-YL`"""
        pipeline = self.create_pipeline()
        _dag, jobs = self.build_dag(pipeline)

        assert "bam-realign-0" in jobs
        cmd_str = str(jobs["bam-realign-0"].shell)
        assert " -Y " in cmd_str
        assert "-YL" not in cmd_str

    def test_lr_fastq_taglist_is_used(self):
        """`--lr_fastq_taglist` reaches the minimap2 command"""
        pipeline = self.create_pipeline()
        pipeline.lr_fastq_taglist = "RG,MM"
        _dag, jobs = self.build_dag(pipeline)

        assert "-T RG,MM" in str(jobs["bam-realign-0"].shell)

    def test_long_read_jobs_wait_for_realignment(self):
        """mosdepth and LongReadSV read the realigned long reads"""
        pipeline = self.create_pipeline()
        dag, jobs = self.build_dag(pipeline)

        assert self.dep_names(dag, jobs["mosdepth-0"]) == {"bam-realign-0"}
        assert self.dep_names(dag, jobs["LongReadSV"]) == {"bam-realign-0"}

    def test_realigned_long_reads_are_called(self):
        """The realigned alignment, not the input, is passed downstream"""
        pipeline = self.create_pipeline()
        _dag, jobs = self.build_dag(pipeline)

        realigned = str(self.mock_vcf).replace(".vcf.gz", "_mm2_sorted_0.cram")
        for name in ("mosdepth-0", "LongReadSV", "dnascope-1"):
            assert realigned in str(jobs[name].shell)


class TestDNAscopeHybridBases:
    """The hybrid pipeline no longer inherits from the other pipelines"""

    def test_mro(self):
        assert DNAscopeHybridPipeline.__mro__ == (
            DNAscopeHybridPipeline,
            BasePipeline,
            ABC,
            object,
        )

    def test_not_a_dnascope_or_longread_pipeline(self):
        with patch("sys.exit"):
            pipeline = DNAscopeHybridPipeline()

        assert not isinstance(pipeline, DNAscopePipeline)
        assert not isinstance(pipeline, DNAscopeLRPipeline)

    @pytest.mark.parametrize(
        "attribute",
        [
            "sample_input",
            "r1_fastq",
            "r2_fastq",
            "readgroups",
            "duplicate_marking",
            "align",
            "collate_align",
            "consensus",
            "fastq",
            "fastq_taglist",
            "haploid_bed",
            "input_ref",
            "interval_padding",
            "no_ramdisk",
            "pcr_free",
            "repeat_model",
            "skip_small_variants",
            "tech",
            "use_pbsv",
            "cnv_excluded_regions",
            "assay",
        ],
    )
    def test_parent_only_attributes_are_gone(self, attribute):
        """The attributes that came from the parent `__init__`s"""
        with patch("sys.exit"):
            pipeline = DNAscopeHybridPipeline()

        assert not hasattr(pipeline, attribute)
