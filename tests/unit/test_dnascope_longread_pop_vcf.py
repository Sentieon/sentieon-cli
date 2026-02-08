"""
Unit tests for DNAscopeLRPipeline pop_vcf logic
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

from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.dag import DAG


class TestDNAscopeLRPopVcf:
    """Test DNAscope longread pipeline pop_vcf logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        # Create mock files
        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_aln = [self.mock_dir / "input.bam"]
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_pop_vcf = self.mock_dir / "population.vcf.gz"
        self.mock_bed = self.mock_dir / "interval.bed"

        # Create empty files
        for file_path in [
            self.mock_ref,
            self.mock_aln[0],
            self.mock_bundle,
            self.mock_pop_vcf,
            self.mock_bed,
        ]:
            file_path.touch()

    def create_pipeline(self):
        """Create a DNAscopeLRPipeline for testing"""
        with patch("sys.exit"):
            pipeline = DNAscopeLRPipeline()
        
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

        # Mock sample inputs
        pipeline.sample_input = self.mock_aln
        pipeline.lr_aln = self.mock_aln
        
        return pipeline

    @patch("sentieon_cli.dnascope_longread.vcf_id")
    def test_validation_requires_pop_vcf(self, mock_vcf_id):
        """Test that validation fails if model bundle requires pop_vcf but none provided"""
        pipeline = self.create_pipeline()
        
        # Mock archive loading
        with patch('sentieon_cli.dnascope_longread.ar_load') as mock_ar_load:
            # Mock bundle info
            bundle_info = {
                "platform": "HiFi",
                "minScriptVersion": "1.5.2",
                "pipeline": "DNAscope LongRead",
                "SentieonVcfID": "some_value",
            }
            bundle_members = [
                "diploid_hp_model",
                "diploid_model",
                "diploid_model_unphased",
                "gvcf_model",
                "haploid_hp_model",
                "haploid_model",
                "longreadsv.model",
                "minimap2.model",
            ]
            mock_ar_load.side_effect = [
                bundle_members,  # First call for bundle members
                json.dumps(bundle_info).encode(),  # Second call for bundle_info.json
            ]
            
            pipeline.dry_run = False
            
            with patch("sys.exit") as mock_exit:
                pipeline.validate_bundle()
                mock_exit.assert_called_with(2)

    def test_validation_pop_vcf_mismatch(self):
        """Test validation fails if pop_vcf ID mismatches"""
        pipeline = self.create_pipeline()
        pipeline.pop_vcf = self.mock_pop_vcf
            
        pipeline.dry_run = False

        with patch("sentieon_cli.dnascope_longread.vcf_id", return_value="wrong_id"), \
            patch("sys.exit") as mock_exit, \
            patch('sentieon_cli.dnascope_longread.ar_load') as mock_ar_load:

            bundle_info = {
                "platform": "HiFi",
                "minScriptVersion": "1.5.2",
                "pipeline": "DNAscope LongRead"
            }
            bundle_members = [
                "diploid_hp_model",
                "diploid_model",
                "diploid_model_unphased",
                "gvcf_model",
                "haploid_hp_model",
                "haploid_model",
                "longreadsv.model",
                "minimap2.model",
            ]
            mock_ar_load.side_effect = [
                bundle_members,  # First call for bundle members
                json.dumps(bundle_info).encode(),  # Second call for bundle_info.json
            ]

            pipeline.validate_bundle()
            mock_exit.assert_called_with(2)

    def test_transfer_jobs_creation(self):
        """Test that transfer jobs are added to the DAG when pop_vcf is present"""
        pipeline = self.create_pipeline()
        pipeline.pop_vcf = self.mock_pop_vcf
        
        # Mock methods to avoid real execution/checking
        with patch("sentieon_cli.dnascope_longread.check_version", return_value=True), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.lr_align_inputs", return_value=([], set())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.lr_align_fastq", return_value=([], set())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.mosdepth", return_value=set()), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.merge_input_files", return_value=(pathlib.Path("merged.bam"), MagicMock())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.pbsv", return_value=(MagicMock(), MagicMock())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.hificnv", return_value=MagicMock()), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.call_svs", return_value=MagicMock()):
             
             dag = pipeline.build_dag()
        
        # Verify transfer jobs exist
        job_names = [job.name for job in dag.waiting_jobs]
        
        # Expect merge-trim jobs
        assert any("merge-trim" in name for name in job_names)
        assert "merge-trim-concat" in job_names

    def test_model_apply_dependency(self):
        """Test that model apply uses the transfer output when pop_vcf is present"""
        pipeline = self.create_pipeline()
        pipeline.pop_vcf = self.mock_pop_vcf
        
        with patch("sentieon_cli.dnascope_longread.check_version", return_value=True), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.lr_align_inputs", return_value=([], set())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.lr_align_fastq", return_value=([], set())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.mosdepth", return_value=set()), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.merge_input_files", return_value=(pathlib.Path("merged.bam"), MagicMock())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.pbsv", return_value=(MagicMock(), MagicMock())), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.hificnv", return_value=MagicMock()), \
             patch("sentieon_cli.dnascope_longread.DNAscopeLRPipeline.call_svs", return_value=MagicMock()):
             
             dag = pipeline.build_dag()
        
        # Find first-modelapply job
        apply_job = None
        for job in dag.waiting_jobs:
            if job.name == "first-modelapply":
                apply_job = job
                break
        
        assert apply_job is not None
        
        # Check dependency: apply_job should depend on merge-trim-concat (transfer_vcf_job)
        deps = dag.waiting_jobs[apply_job]
        dep_names = [j.name for j in deps]
        assert "merge-trim-concat" in dep_names
