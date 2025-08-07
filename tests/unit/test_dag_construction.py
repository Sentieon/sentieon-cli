"""
Unit tests for DAG construction logic
"""

import pathlib
import pytest
import tempfile
import sys
import os
from unittest.mock import patch

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.dag import DAG
from sentieon_cli.job import Job


class TestDAGConstruction:
    """Test DAG construction for pipelines"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()

        # Create mock files
        self.mock_vcf = pathlib.Path(self.temp_dir) / "output.vcf.gz"
        self.mock_ref = pathlib.Path(self.temp_dir) / "reference.fa"
        self.mock_bam = pathlib.Path(self.temp_dir) / "sample.bam"
        self.mock_bundle = pathlib.Path(self.temp_dir) / "model.bundle"

        # Create empty files
        for file_path in [self.mock_ref, self.mock_bam, self.mock_bundle]:
            file_path.touch()

        # Create reference index
        (pathlib.Path(str(self.mock_ref) + ".fai")).touch()

    def create_basic_dnascope_pipeline(self):
        """Create a basic DNAscope pipeline for testing"""
        pipeline = DNAscopePipeline()

        # Setup logging first
        from tests.utils.test_helpers import create_mock_args
        args = create_mock_args()
        pipeline.setup_logging(args)

        pipeline.output_vcf = self.mock_vcf
        pipeline.reference = self.mock_ref
        pipeline.model_bundle = self.mock_bundle
        pipeline.sample_input = [self.mock_bam]
        pipeline.sr_r1_fastq = []
        pipeline.sr_r2_fastq = []
        pipeline.sr_readgroups = []
        pipeline.sr_duplicate_marking = "markdup"
        pipeline.cores = 2
        pipeline.dry_run = True  # Important: use dry run to avoid actual execution
        pipeline.skip_version_check = True
        pipeline.tmp_dir = pathlib.Path(self.temp_dir)

        return pipeline

    @patch('sentieon_cli.util.library_preloaded')
    def test_basic_dag_structure(self, mock_lib_preloaded):
        """Test basic DAG structure is created correctly"""
        mock_lib_preloaded.return_value = True

        pipeline = self.create_basic_dnascope_pipeline()
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Verify DAG is created
        assert isinstance(dag, DAG)

        # Should have jobs for dedup, variant calling, etc.
        total_jobs = (len(dag.waiting_jobs) +
                     len(dag.ready_jobs) +
                     len(dag.finished_jobs))
        assert total_jobs > 0

    @patch('sentieon_cli.util.library_preloaded')
    def test_dag_job_dependencies(self, mock_lib_preloaded):
        """Test that jobs have proper dependencies"""
        mock_lib_preloaded.return_value = True

        pipeline = self.create_basic_dnascope_pipeline()
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Extract all jobs
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())

        # Should have some ready jobs (no dependencies)
        assert len(dag.ready_jobs) > 0

        # Check that waiting jobs have valid dependencies
        for job, dependencies in dag.waiting_jobs.items():
            for dep in dependencies:
                assert dep in all_jobs, f"Dependency {dep} not found in job list"

    @patch('sentieon_cli.util.library_preloaded')
    def test_conditional_job_creation(self, mock_lib_preloaded):
        """Test that conditional jobs are created based on parameters"""
        mock_lib_preloaded.return_value = True

        # Test with metrics enabled
        pipeline1 = self.create_basic_dnascope_pipeline()
        pipeline1.skip_metrics = False
        pipeline1.validate()
        pipeline1.configure()
        dag1 = pipeline1.build_dag()

        # Test with metrics disabled
        pipeline2 = self.create_basic_dnascope_pipeline()
        pipeline2.skip_metrics = True
        pipeline2.validate()
        pipeline2.configure()
        dag2 = pipeline2.build_dag()

        # Should have different numbers of jobs
        total_jobs1 = (len(dag1.waiting_jobs) +
                      len(dag1.ready_jobs) +
                      len(dag1.finished_jobs))
        total_jobs2 = (len(dag2.waiting_jobs) +
                      len(dag2.ready_jobs) +
                      len(dag2.finished_jobs))

        # Pipeline with metrics should have more jobs
        assert total_jobs1 > total_jobs2

    @patch('sentieon_cli.util.library_preloaded')
    def test_skip_small_variants_dag(self, mock_lib_preloaded):
        """Test DAG when small variant calling is skipped"""
        mock_lib_preloaded.return_value = True

        pipeline = self.create_basic_dnascope_pipeline()
        pipeline.skip_small_variants = True
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Extract job names
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        job_names = [job.name for job in all_jobs]

        # Should not have variant calling jobs
        variant_job_keywords = ["variant-calling", "model-apply", "gvcftyper"]
        for keyword in variant_job_keywords:
            assert not any(keyword in name for name in job_names)

    @patch('sentieon_cli.util.library_preloaded')
    def test_gvcf_mode_dag_differences(self, mock_lib_preloaded):
        """Test DAG differences when gVCF mode is enabled"""
        mock_lib_preloaded.return_value = True

        # Test without gVCF
        pipeline1 = self.create_basic_dnascope_pipeline()
        pipeline1.gvcf = False
        pipeline1.validate()
        pipeline1.configure()
        dag1 = pipeline1.build_dag()

        # Test with gVCF
        pipeline2 = self.create_basic_dnascope_pipeline()
        pipeline2.gvcf = True
        pipeline2.validate()
        pipeline2.configure()
        dag2 = pipeline2.build_dag()

        # Extract job names
        jobs1 = list(dag1.waiting_jobs.keys()) + list(dag1.ready_jobs.keys())
        jobs2 = list(dag2.waiting_jobs.keys()) + list(dag2.ready_jobs.keys())
        job_names1 = [job.name for job in jobs1]
        job_names2 = [job.name for job in jobs2]

        # gVCF mode should have gvcftyper job
        assert not any("gvcftyper" in name for name in job_names1)
        assert any("gvcftyper" in name for name in job_names2)


class TestDAGJobProperties:
    """Test properties of individual jobs in the DAG"""

    def test_job_thread_allocation(self):
        """Test that jobs have appropriate thread allocation"""
        # Test high-CPU jobs get multiple threads
        job = Job("sentieon driver --algo DNAscope", "variant-calling", 8)
        assert job.threads == 8

        # Test lightweight jobs get fewer threads
        job = Job("rm temp_file.vcf", "cleanup", 0)
        assert job.threads == 0

    def test_job_resource_requirements(self):
        """Test that jobs specify resource requirements correctly"""
        # Test NUMA-aware job
        job = Job("align-command", "alignment", 4, resources={"node0": 1})
        assert job.resources == {"node0": 1}

        # Test job without special resources
        job = Job("simple-command", "simple", 2)
        assert job.resources == {}

    def test_job_failure_tolerance(self):
        """Test that cleanup jobs are marked as failure-tolerant"""
        # Cleanup jobs should be failure tolerant
        job = Job("rm temp_files", "cleanup", 0, fail_ok=True)
        assert job.fail_ok is True

        # Critical jobs should not be failure tolerant
        job = Job("sentieon driver", "variant-calling", 4, fail_ok=False)
        assert job.fail_ok is False


class TestLongReadDAGConstruction:
    """Test DAG construction for long-read pipeline"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()

        # Create mock files
        self.mock_vcf = pathlib.Path(self.temp_dir) / "output.vcf.gz"
        self.mock_ref = pathlib.Path(self.temp_dir) / "reference.fa"
        self.mock_bam = pathlib.Path(self.temp_dir) / "longread.bam"
        self.mock_bundle = pathlib.Path(self.temp_dir) / "model.bundle"

        # Create empty files
        for file_path in [self.mock_ref, self.mock_bam, self.mock_bundle]:
            file_path.touch()

    def create_basic_longread_pipeline(self):
        """Create a basic long-read pipeline for testing"""
        pipeline = DNAscopeLRPipeline()

        # Setup logging first
        from tests.utils.test_helpers import create_mock_args
        args = create_mock_args()
        pipeline.setup_logging(args)

        pipeline.output_vcf = self.mock_vcf
        pipeline.reference = self.mock_ref
        pipeline.model_bundle = self.mock_bundle
        pipeline.sample_input = [self.mock_bam]
        pipeline.fastq = []
        pipeline.readgroups = []
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = pathlib.Path(self.temp_dir)

        return pipeline

    @patch('sentieon_cli.util.library_preloaded')
    def test_longread_dag_complexity(self, mock_lib_preloaded):
        """Test that long-read DAG is more complex due to phasing"""
        mock_lib_preloaded.return_value = True

        pipeline = self.create_basic_longread_pipeline()
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Long-read pipeline should have many jobs due to phasing complexity
        total_jobs = (len(dag.waiting_jobs) +
                     len(dag.ready_jobs) +
                     len(dag.finished_jobs))

        # Should have multiple phases of variant calling
        assert total_jobs > 10  # Rough estimate for complex phased calling

    @patch('sentieon_cli.util.library_preloaded')
    def test_ont_vs_hifi_dag_differences(self, mock_lib_preloaded):
        """Test DAG differences between ONT and HiFi technologies"""
        mock_lib_preloaded.return_value = True

        # Test HiFi pipeline
        pipeline_hifi = self.create_basic_longread_pipeline()
        pipeline_hifi.tech = "HiFi"
        pipeline_hifi.validate()
        pipeline_hifi.configure()
        dag_hifi = pipeline_hifi.build_dag()

        # Test ONT pipeline
        pipeline_ont = self.create_basic_longread_pipeline()
        pipeline_ont.tech = "ONT"
        pipeline_ont.validate()
        pipeline_ont.configure()
        dag_ont = pipeline_ont.build_dag()

        # Extract all jobs
        jobs_hifi = list(dag_hifi.waiting_jobs.keys()) + list(dag_hifi.ready_jobs.keys())
        jobs_ont = list(dag_ont.waiting_jobs.keys()) + list(dag_ont.ready_jobs.keys())

        # ONT should skip CNV calling (fewer jobs)
        assert len(jobs_hifi) >= len(jobs_ont)


class TestDAGValidation:
    """Test DAG validation and correctness"""

    def test_dag_has_no_cycles(self):
        """Test that generated DAGs have no cycles"""
        dag = DAG()

        # Create test jobs
        job1 = Job("command1", "job1")
        job2 = Job("command2", "job2")
        job3 = Job("command3", "job3")

        # Add jobs with dependencies: job1 -> job2 -> job3
        dag.add_job(job1)
        dag.add_job(job2, {job1})
        dag.add_job(job3, {job2})

        # This should be a valid DAG
        assert job1 in dag.ready_jobs
        assert job2 in dag.waiting_jobs
        assert job3 in dag.waiting_jobs

        # Dependencies should be correct
        assert dag.waiting_jobs[job2] == {job1}
        assert dag.waiting_jobs[job3] == {job2}

    def test_dag_execution_simulation(self):
        """Test simulated DAG execution"""
        dag = DAG()

        # Create test jobs
        job1 = Job("command1", "job1")
        job2 = Job("command2", "job2")
        job3 = Job("command3", "job3")

        # Add jobs: job1 and job2 can run in parallel, then job3
        dag.add_job(job1)
        dag.add_job(job2)
        dag.add_job(job3, {job1, job2})

        # Simulate execution
        dag_gen = dag.update_dag()
        ready_jobs = dag_gen.send(None)

        # Initially job1 and job2 should be ready
        assert job1 in ready_jobs
        assert job2 in ready_jobs
        assert job3 not in ready_jobs

        # Complete job1
        new_ready = dag_gen.send(job1)
        assert len(new_ready) == 0  # job3 still waiting for job2

        # Complete job2
        new_ready = dag_gen.send(job2)
        assert job3 in new_ready  # Now job3 is ready


if __name__ == "__main__":
    pytest.main([__file__])
