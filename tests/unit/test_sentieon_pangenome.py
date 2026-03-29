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

    def test_call_svs(self):
        """Test that PangenomeSV is added when --call_svs is enabled"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        assert isinstance(dag, DAG)

        # The dnascope-raw job should exist
        all_jobs = list(dag.waiting_jobs.keys()) + list(
            dag.ready_jobs.keys()
        )
        job_names = [job.name for job in all_jobs]
        assert "dnascope-raw" in job_names

        # Find the dnascope-raw job and check its command includes
        # both DNAscope and PangenomeSV
        dnascope_job = None
        for job in all_jobs:
            if job.name == "dnascope-raw":
                dnascope_job = job
                break
        assert dnascope_job is not None
        cmd_str = str(dnascope_job.shell)
        assert "--algo DNAscope" in cmd_str
        assert "--algo PangenomeSV" in cmd_str
        assert "--gfa_file" in cmd_str

        # SV output should use _sv.vcf.gz suffix
        sv_vcf = str(pipeline.output_vcf).replace(
            ".vcf.gz", "_sv.vcf.gz"
        )
        assert sv_vcf in cmd_str

    def test_call_svs_disabled_by_default(self):
        """Test that PangenomeSV is not added by default"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_first_dag()

        all_jobs = list(dag.waiting_jobs.keys()) + list(
            dag.ready_jobs.keys()
        )
        dnascope_job = None
        for job in all_jobs:
            if job.name == "dnascope-raw":
                dnascope_job = job
                break
        assert dnascope_job is not None
        cmd_str = str(dnascope_job.shell)
        assert "--algo DNAscope" in cmd_str
        assert "--algo PangenomeSV" not in cmd_str

    def test_skip_small_variants(self):
        """Test that DNAscope, transfer, and model-apply are skipped"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        dag = pipeline.build_first_dag()

        assert isinstance(dag, DAG)

        all_jobs = list(dag.waiting_jobs.keys()) + list(
            dag.ready_jobs.keys()
        )
        job_names = [job.name for job in all_jobs]
        assert "dnascope-raw" not in job_names
        assert "model-apply" not in job_names
        assert "merge-trim-concat" not in job_names

    def test_skip_small_variants_with_call_svs(self):
        """Test SV-only mode: driver runs PangenomeSV without DNAscope"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        assert isinstance(dag, DAG)

        all_jobs = list(dag.waiting_jobs.keys()) + list(
            dag.ready_jobs.keys()
        )
        job_names = [job.name for job in all_jobs]

        # The driver job should still run for SV calling
        assert "dnascope-raw" in job_names

        # Transfer and model-apply should be skipped
        assert "model-apply" not in job_names
        assert "merge-trim-concat" not in job_names

        # The driver command should have PangenomeSV but NOT DNAscope
        dnascope_job = None
        for job in all_jobs:
            if job.name == "dnascope-raw":
                dnascope_job = job
                break
        assert dnascope_job is not None
        cmd_str = str(dnascope_job.shell)
        assert "--algo PangenomeSV" in cmd_str
        assert "--algo DNAscope" not in cmd_str

    def _get_all_job_names(self, dag):
        """Helper to get all job names from a DAG"""
        all_jobs = list(dag.waiting_jobs.keys()) + list(
            dag.ready_jobs.keys()
        )
        return [job.name for job in all_jobs], all_jobs

    def test_cnv_jobs_with_call_svs(self):
        """Test that CNV jobs are added when --call_svs is enabled"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        job_names, all_jobs = self._get_all_job_names(dag)
        assert "cnvscope" in job_names
        assert "cnv-model-apply" in job_names
        assert "indel2cnv" in job_names
        assert "combine-cnv" in job_names

    def test_cnv_jobs_not_present_without_call_svs(self):
        """Test that CNV jobs are not added by default"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_first_dag()

        job_names, _ = self._get_all_job_names(dag)
        assert "cnvscope" not in job_names
        assert "cnv-model-apply" not in job_names
        assert "indel2cnv" not in job_names
        assert "combine-cnv" not in job_names

    def test_cnv_jobs_with_skip_small_variants(self):
        """Test CNV jobs are added in SV-only mode"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        job_names, _ = self._get_all_job_names(dag)
        assert "cnvscope" in job_names
        assert "cnv-model-apply" in job_names
        assert "indel2cnv" in job_names
        assert "combine-cnv" in job_names

    def test_cnvscope_command(self):
        """Test CNVscope driver command has correct algo and model"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        _, all_jobs = self._get_all_job_names(dag)
        cnvscope_job = None
        for job in all_jobs:
            if job.name == "cnvscope":
                cnvscope_job = job
                break
        assert cnvscope_job is not None
        cmd_str = str(cnvscope_job.shell)
        assert "--algo CNVscope" in cmd_str
        assert "cnvscope.model" in cmd_str
        # Input should be the sample BAM (BAM input mode)
        assert str(pipeline.sample_input[0]) in cmd_str

    def test_cnv_model_apply_command(self):
        """Test CNVModelApply driver command"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        _, all_jobs = self._get_all_job_names(dag)
        job = None
        for j in all_jobs:
            if j.name == "cnv-model-apply":
                job = j
                break
        assert job is not None
        cmd_str = str(job.shell)
        assert "--algo CNVModelApply" in cmd_str
        assert "cnvscope.model" in cmd_str

    def test_indel2cnv_command(self):
        """Test indel2cnv script command"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        _, all_jobs = self._get_all_job_names(dag)
        job = None
        for j in all_jobs:
            if j.name == "indel2cnv":
                job = j
                break
        assert job is not None
        cmd_str = str(job.shell)
        assert "indel2cnv.py" in cmd_str
        assert str(pipeline.reference) in cmd_str

    def test_combine_cnv_command(self):
        """Test combine_cnv script command"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        dag = pipeline.build_first_dag()

        _, all_jobs = self._get_all_job_names(dag)
        job = None
        for j in all_jobs:
            if j.name == "combine-cnv":
                job = j
                break
        assert job is not None
        cmd_str = str(job.shell)
        assert "combine_cnv.py" in cmd_str
        assert "--cnv" in cmd_str
        assert "--converted" in cmd_str
        # Output should use _cnv.vcf.gz suffix
        cnv_vcf = str(pipeline.output_vcf).replace(
            ".vcf.gz", "_cnv.vcf.gz"
        )
        assert cnv_vcf in cmd_str

    def test_cnv_with_skip_model_apply(self):
        """Test CNV jobs are added even with skip_model_apply"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        pipeline.skip_model_apply = True
        dag = pipeline.build_first_dag()

        job_names, _ = self._get_all_job_names(dag)
        assert "model-apply" not in job_names
        assert "cnvscope" in job_names
        assert "cnv-model-apply" in job_names
        assert "indel2cnv" in job_names
        assert "combine-cnv" in job_names

