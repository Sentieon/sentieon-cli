"""
Unit tests for SentieonPangenome pipeline logic
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

from sentieon_cli.sentieon_pangenome import SentieonPangenome
from sentieon_cli.dag import DAG


class TestSentieonPangenome:
    """Test Sentieon pangenome pipeline logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        # Create mock files
        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_bam = self.mock_dir / "sample.bam"
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_gbz = self.mock_dir / "pangenome.grch38.gbz"
        self.mock_hapl = self.mock_dir / "pangenome.grch38.hapl"
        self.mock_pop_vcf = self.mock_dir / "population.vcf.gz"

        # Create empty files
        for file_path in [
            self.mock_ref,
            self.mock_bam,
            self.mock_bundle,
            self.mock_gbz,
            self.mock_hapl,
            self.mock_pop_vcf,
        ]:
            file_path.touch()

        # Create reference index
        with open(str(self.mock_ref) + ".fai", "w") as f:
            f.write("chr1\t1000\t0\t80\t81\n")

    def create_pipeline(self):
        """Create a SentieonPangenome pipeline for testing"""
        pipeline = SentieonPangenome()
        
        # Setup mocks
        pipeline.logger = MagicMock()
        
        # Configure arguments
        pipeline.output_vcf = self.mock_vcf
        pipeline.reference = self.mock_ref
        pipeline.model_bundle = self.mock_bundle
        pipeline.sample_input = [self.mock_bam]
        pipeline.gbz = self.mock_gbz
        pipeline.hapl = self.mock_hapl
        pipeline.pop_vcf = self.mock_pop_vcf
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.skip_contig_checks = True
        pipeline.skip_pangenome_name_checks = True
        pipeline.skip_pop_vcf_id_check = True
        pipeline.tmp_dir = self.mock_dir
        
        # Mock parsing fai
        pipeline.fai_data = {"chr1": {"length": 1000}}
        pipeline.shards = [MagicMock()]
        pipeline.shards[0].contig = "chr1"
        pipeline.shards[0].start = 1
        pipeline.shards[0].stop = 1000
        
        # Mock pop vcf contigs
        pipeline.pop_vcf_contigs = {"chr1": 1000}

        # Mock collection readgroups
        pipeline.bam_readgroups = [{"ID": "rg1", "SM": "sample1"}]
        pipeline.fastq_readgroup = {}
        pipeline.tech = "Illumina"
        pipeline.pcr_free = False

        return pipeline

    def test_model_apply_default(self):
        """Test that model apply job is created by default"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_first_dag()

        # Verify DAG is created
        assert isinstance(dag, DAG)

        # Check for model-apply job
        job_names = [job.name for job in dag.waiting_jobs.keys()]
        assert "model-apply" in job_names

    def test_skip_model_apply(self):
        """Test that model apply job is skipped when requested"""
        pipeline = self.create_pipeline()
        pipeline.skip_model_apply = True
        dag = pipeline.build_first_dag()

        # Verify DAG is created
        assert isinstance(dag, DAG)

        # Check that model-apply job is NOT present
        job_names = [job.name for job in dag.waiting_jobs]
        assert "model-apply" not in job_names
        
        # Verify that the concat job writes to the final output
        # Find the merge-trim-concat job
        concat_job = None
        for job in dag.waiting_jobs:
            if job.name == "merge-trim-concat":
                concat_job = job
                break
        
        assert concat_job is not None
        # Check that the first argument (output file) is set to the final output VCF
        # Check that the output file is in the arguments
        args_str = " ".join([str(arg) for arg in concat_job.shell.nodes[0].args])
        assert str(pipeline.output_vcf) in args_str

