"""
Unit tests for DNAscopeHybridPipeline pop_vcf logic
"""

import pathlib
import pytest
import tempfile
import sys
import os
import json
from unittest.mock import patch, MagicMock

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline
from sentieon_cli.dag import DAG


class TestDNAscopeHybridPopVcf:
    """Test DNAscope hybrid pipeline pop_vcf logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        # Create mock files
        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_sr_aln = [self.mock_dir / "short.bam"]
        self.mock_lr_aln = [self.mock_dir / "long.bam"]
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_pop_vcf = self.mock_dir / "population.vcf.gz"
        self.mock_bed = self.mock_dir / "interval.bed"

        # Create empty files
        for file_path in [
            self.mock_ref,
            self.mock_sr_aln[0],
            self.mock_lr_aln[0],
            self.mock_bundle,
            self.mock_pop_vcf,
            self.mock_bed,
        ]:
            file_path.touch()

    def create_pipeline(self):
        """Create a DNAscopeHybridPipeline for testing"""
        # We need to mock sys.exit to avoid exiting during initialization if validation fails
        with patch("sys.exit"):
            pipeline = DNAscopeHybridPipeline()
        
        # Setup mocks
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
        pipeline.skip_pop_vcf_id_check = False
        
        # Mock fai related
        pipeline.fai_data = {"chr1": {"length": 1000}}
        pipeline.shards = [MagicMock()]
        pipeline.shards[0].contig = "chr1"
        pipeline.shards[0].start = 1
        pipeline.shards[0].stop = 1000
        
        # Mock pop vcf contigs
        pipeline.pop_vcf_contigs = {"chr1": 1000}

        # Mock readgroup info (needed for build_dag)
        pipeline.lr_aln_readgroups = [[{"ID": "lr_rg1", "SM": "sample1"}]]
        pipeline.sr_aln_readgroups = [[{"ID": "sr_rg1", "SM": "sample1"}]]
        pipeline.hybrid_rg_sm = "sample1"
        pipeline.hybrid_set_rg = False
        pipeline.shortread_tech = "Illumina"
        pipeline.longread_tech = "ONT"
        pipeline.sr_aln = self.mock_sr_aln
        pipeline.lr_aln = self.mock_lr_aln
        
        return pipeline

    @patch("sentieon_cli.dnascope_hybrid.ar_load")
    @patch("sentieon_cli.dnascope_hybrid.vcf_id")
    def test_validation_requires_pop_vcf(self, mock_vcf_id, mock_ar_load):
        """Test that validation fails if model bundle requires pop_vcf but none provided"""
        pipeline = self.create_pipeline()
        
        # Mock bundle info requiring SentieonVcfID
        mock_bundle_info = {
            "longReadPlatform": "ONT",
            "shortReadPlatform": "Illumina",
            "SentieonVcfID": "some_id",
            "minScriptVersion": "2.0"
        }
        mock_ar_load.return_value = json.dumps(mock_bundle_info).encode()
        pipeline.dry_run = False
        
        with patch("sys.exit") as mock_exit:
             pipeline.validate_bundle()
             mock_exit.assert_called_with(2)

    @patch("sentieon_cli.dnascope_hybrid.ar_load")
    def test_validation_pop_vcf_mismatch(self, mock_ar_load):
        """Test validation fails if pop_vcf ID mismatches"""
        pipeline = self.create_pipeline()
        pipeline.pop_vcf = self.mock_pop_vcf
        
        mock_bundle_info = {
            "longReadPlatform": "ONT",
            "shortReadPlatform": "Illumina",
            "SentieonVcfID": "expected_id",
            "minScriptVersion": "2.0"
        }
        mock_ar_load.return_value = json.dumps(mock_bundle_info).encode()
        pipeline.dry_run = False

        with patch("sentieon_cli.dnascope_hybrid.vcf_id", return_value="wrong_id"), \
             patch("sys.exit") as mock_exit:
                 
            pipeline.validate_bundle()
            mock_exit.assert_called_with(2)

    def test_transfer_jobs_creation(self):
        """Test that transfer jobs are added to the DAG when pop_vcf is present"""
        pipeline = self.create_pipeline()
        pipeline.pop_vcf = self.mock_pop_vcf
        
        # Build DAG
        with patch("sentieon_cli.dnascope_hybrid.check_version", return_value=True):
             dag = pipeline.build_dag()
        
        # Verify transfer jobs exist
        job_names = [job.name for job in dag.waiting_jobs]
        
        # Expect merge-trim jobs
        assert any("merge-trim" in name for name in job_names)
        assert "merge-trim-concat" in job_names

    def test_model_apply_input_with_pop_vcf(self):
        """Test that model apply uses the transfer output when pop_vcf is present"""
        pipeline = self.create_pipeline()
        pipeline.pop_vcf = self.mock_pop_vcf
        
        with patch("sentieon_cli.dnascope_hybrid.check_version", return_value=True):
             dag = pipeline.build_dag()
        
        # Find model-apply job
        apply_job = None
        for job in dag.waiting_jobs:
            if job.name == "model-apply":
                apply_job = job
                break
        
        assert apply_job is not None
        
        # Check dependency: apply_job should depend on merge-trim-concat
        deps = dag.waiting_jobs[apply_job]
        dep_names = [j.name for j in deps]
        assert "merge-trim-concat" in dep_names
