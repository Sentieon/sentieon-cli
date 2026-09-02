"""
Unit tests for the DNAscope hybrid first-stage realignment jobs
"""

import os
import pathlib
import sys
import tempfile
from unittest.mock import MagicMock, patch

import packaging.version

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.dnascope_hybrid import (
    CALLING_MIN_VERSIONS,
    DNAscopeHybridPipeline,
)


class TestDNAscopeHybridStage1:
    """Test the first-stage haplotype alignment jobs"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_sr_aln = [self.mock_dir / "short.bam"]
        self.mock_lr_aln = [self.mock_dir / "long.bam"]
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_bed = self.mock_dir / "interval.bed"

        for file_path in [
            self.mock_ref,
            self.mock_sr_aln[0],
            self.mock_lr_aln[0],
            self.mock_bundle,
            self.mock_bed,
        ]:
            file_path.touch()

        self.stage1_fifo = self.mock_dir / "stage1_hap.fq"
        self.stage1_hap_bam = self.mock_dir / "stage1_hap.bam"

    def create_pipeline(self):
        """Create a DNAscopeHybridPipeline for testing"""
        with patch("sys.exit"):
            pipeline = DNAscopeHybridPipeline()

        pipeline.logger = MagicMock()

        # Configure arguments
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

        return pipeline

    def build_dag(self, pipeline):
        """Build the DAG and index the jobs by name"""
        with patch(
            "sentieon_cli.dnascope_hybrid.check_version", return_value=True
        ):
            dag = pipeline.build_dag()
        jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        return dag, {job.name: job for job in jobs}

    def dep_names(self, dag, job):
        """The names of a job's dependencies"""
        return {dep.name for dep in dag.waiting_jobs.get(job, set())}

    def test_fifo_job(self):
        """The haplotype fastq fifo is created before the stage1 jobs"""
        pipeline = self.create_pipeline()
        dag, jobs = self.build_dag(pipeline)

        assert "stage1-fifo" in jobs
        cmd_str = str(jobs["stage1-fifo"].shell)
        assert cmd_str == f"mkfifo {self.stage1_fifo}"
        assert self.dep_names(dag, jobs["stage1-fifo"]) == set()

    def test_stage1_hap_command(self):
        """The haplotype alignments are sorted into stage1_hap.bam"""
        pipeline = self.create_pipeline()
        _dag, jobs = self.build_dag(pipeline)

        assert "first-stage-hap" in jobs
        hap_job = jobs["first-stage-hap"]
        cmd_str = str(hap_job.shell)
        assert "--algo HybridStage1" in cmd_str
        assert "HybridStage1.model" in cmd_str
        # The unsorted haplotype BAM is written to stdout and sorted
        assert "--hap_bam -" in cmd_str
        assert "sentieon util sort" in cmd_str
        assert f"-o {self.stage1_hap_bam}" in cmd_str
        # The driver writes BAM, not SAM, so no conversion is needed
        assert "--sam2bam" not in cmd_str
        # The fastq output goes to the fifo
        assert str(self.stage1_fifo) in cmd_str

        # The job must run concurrently with the fifo reader
        assert hap_job.threads == 0

    def test_stage1_reads_the_fifo(self):
        """The bwa job reads the haplotype fastq from the fifo"""
        pipeline = self.create_pipeline()
        _dag, jobs = self.build_dag(pipeline)

        cmd_str = str(jobs["first-stage"].shell)
        # The fifo is a plain `cat` argument, read before the proc sub
        assert cmd_str.startswith(f"cat {self.stage1_fifo} ")
        assert f"<({self.stage1_fifo}" not in cmd_str
        # The haplotype driver moved to its own job; only the insertion
        # driver remains in the bwa pipeline
        assert "--hap_bam" not in cmd_str
        assert "HybridStage1.model" not in cmd_str
        assert "HybridStage1_ins.model" in cmd_str
        assert "sentieon bwa mem" in cmd_str
        assert "--sam2bam" in cmd_str

    def test_stage1_dependencies(self):
        """The stage1 jobs wait for the fifo and the merged BED"""
        pipeline = self.create_pipeline()
        dag, jobs = self.build_dag(pipeline)

        assert self.dep_names(dag, jobs["first-stage-hap"]) == {
            "stage1-fifo",
            "concat-merge-bed",
        }
        assert self.dep_names(dag, jobs["first-stage"]) == {
            "stage1-fifo",
            "concat-merge-bed",
        }
        # The second stage reads both stage1 outputs
        assert self.dep_names(dag, jobs["second-stage"]) == {
            "first-stage",
            "first-stage-hap",
        }
        # The haplotype VCF removed here is written by the haplotype job
        assert self.dep_names(dag, jobs["rm-tmp2"]) == {
            "first-stage",
            "first-stage-hap",
        }

    def test_cleanup_jobs(self):
        """Each cleanup job waits for the jobs writing what it removes"""
        pipeline = self.create_pipeline()
        dag, jobs = self.build_dag(pipeline)

        assert self.dep_names(dag, jobs["rm-tmp1"]) == {"concat-merge-bed"}
        assert self.dep_names(dag, jobs["rm-tmp2"]) == {
            "first-stage",
            "first-stage-hap",
        }
        assert self.dep_names(dag, jobs["rm-tmp3"]) == {"second-stage"}
        assert self.dep_names(dag, jobs["rm-tmp4"]) == {"third-stage"}
        assert self.dep_names(dag, jobs["rm-tmp5"]) == {"concat-calls"}

    def test_retain_tmpdir_skips_cleanup(self):
        """`--retain_tmpdir` keeps the intermediate files"""
        pipeline = self.create_pipeline()
        pipeline.retain_tmpdir = True
        _dag, jobs = self.build_dag(pipeline)

        assert [name for name in jobs if name.startswith("rm-tmp")] == []

    def test_second_stage_inputs(self):
        """The second stage consumes the sorted haplotype BAM"""
        pipeline = self.create_pipeline()
        _dag, jobs = self.build_dag(pipeline)

        cmd_str = str(jobs["second-stage"].shell)
        assert "--algo HybridStage2" in cmd_str
        assert f"--input {self.stage1_hap_bam}" in cmd_str
        assert f"--input {self.mock_dir / 'hybrid_stage1.bam'}" in cmd_str

    def test_driver_min_version(self):
        """The unsorted `--hap_bam` output requires 202503.04"""
        assert CALLING_MIN_VERSIONS[
            "sentieon driver"
        ] >= packaging.version.Version("202503.04")
