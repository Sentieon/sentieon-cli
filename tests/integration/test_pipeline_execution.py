"""
Integration tests for pipeline execution
"""

import pytest
import sys
import os
from unittest.mock import patch

# Add the parent directory to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from tests.utils.test_helpers import (
    CommandValidator,
    DAGAnalyzer,
    setup_basic_test_environment,
    teardown_test_environment
)


class TestPipelineDryRunExecution:
    """Test pipeline execution in dry-run mode"""

    def setup_method(self):
        """Setup test environment"""
        self.helper = setup_basic_test_environment()

    def teardown_method(self):
        """Cleanup test environment"""
        teardown_test_environment(self.helper)

    def test_dnascope_dry_run_execution(self):
        """Test DNAscope pipeline dry-run execution"""
        pipeline = self.helper.create_dnascope_pipeline(use_bam=True)
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Extract all commands that would be executed
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        commands = [job.shell for job in all_jobs]

        # Validate commands
        sentieon_commands = [cmd for cmd in commands if cmd.startswith("sentieon")]
        assert len(sentieon_commands) > 0, "Should have sentieon commands"

        for cmd in sentieon_commands:
            assert CommandValidator.validate_sentieon_command(cmd), f"Invalid command: {cmd}"

    def test_longread_dry_run_execution(self):
        """Test long-read pipeline dry-run execution"""
        pipeline = self.helper.create_longread_pipeline()
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Extract all commands
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        commands = [job.shell for job in all_jobs]

        # Should have complex phased calling commands
        variant_commands = [cmd for cmd in commands if "DNAscope" in cmd or "VariantPhaser" in cmd]
        assert len(variant_commands) > 1, "Should have multiple variant calling steps"

    def test_hybrid_dry_run_execution(self):
        """Test hybrid pipeline dry-run execution"""
        pipeline = self.helper.create_hybrid_pipeline()

        # Mock the archive loading and readgroup functions
        with patch('sentieon_cli.dnascope_hybrid.ar_load') as mock_ar_load, \
             patch('sentieon_cli.command_strings.get_rg_lines') as mock_get_rg:

            # Mock bundle info
            import json
            bundle_info = {
                "longReadPlatform": "HiFi",
                "shortReadPlatform": "Illumina",
                "minScriptVersion": "1.0.0",
                "pipeline": "DNAscope Hybrid"
            }
            mock_ar_load.side_effect = [
                json.dumps(bundle_info).encode(),  # First call for bundle_info.json
                ["longreadsv.model", "cnv.model", "bwa.model"]  # Second call for bundle members
            ]

            # Mock readgroup lines
            mock_get_rg.return_value = ["@RG\tID:test\tSM:sample"]

            pipeline.validate()
            pipeline.configure()

            dag = pipeline.build_dag()

            # Should have hybrid-specific jobs
            job_names = [job.name for job in list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())]

            # Look for hybrid-specific stages
            hybrid_stages = [name for name in job_names if "stage" in name.lower()]
            assert len(hybrid_stages) > 0, "Should have hybrid staging jobs"


class TestPipelineCommandGeneration:
    """Test command generation for different pipeline configurations"""

    def setup_method(self):
        """Setup test environment"""
        self.helper = setup_basic_test_environment()

    def teardown_method(self):
        """Cleanup test environment"""
        teardown_test_environment(self.helper)

    def test_dnascope_command_parameters(self):
        """Test DNAscope command parameter generation"""
        pipeline = self.helper.create_dnascope_pipeline(
            use_bam=True,
            dbsnp=self.helper.fs.create_file("dbsnp.vcf.gz"),
            pcr_free=True,
            gvcf=True
        )
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Find variant calling job
        variant_job = DAGAnalyzer.get_job_by_name(dag, "variant-calling")
        assert variant_job is not None, "Should have variant calling job"

        # Check that command includes expected parameters
        command = variant_job.shell
        assert "--algo DNAscope" in command
        assert "dbsnp.vcf.gz" in command, "Should include dbSNP file"
        assert "emit_mode" in command, "Should have emit mode for gVCF"

    def test_resource_allocation_in_commands(self):
        """Test that resource allocation is reflected in commands"""
        pipeline = self.helper.create_dnascope_pipeline(
            use_bam=True,
            cores=8
        )
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Check thread allocation
        thread_counts = DAGAnalyzer.count_jobs_by_thread_usage(dag)

        # Should have jobs with different thread counts
        assert len(thread_counts) > 1, "Should have jobs with different thread requirements"

        # High-CPU jobs should use multiple threads
        high_cpu_jobs = [job for job in list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
                        if job.threads >= 4]
        assert len(high_cpu_jobs) > 0, "Should have high-CPU jobs"

    def test_conditional_algorithm_inclusion(self):
        """Test conditional inclusion of algorithms based on parameters"""
        # Test with metrics enabled
        pipeline1 = self.helper.create_dnascope_pipeline(use_bam=True, skip_metrics=False)
        pipeline1.validate()
        pipeline1.configure()
        dag1 = pipeline1.build_dag()

        # Test with metrics disabled
        pipeline2 = self.helper.create_dnascope_pipeline(use_bam=True, skip_metrics=True)
        pipeline2.validate()
        pipeline2.configure()
        dag2 = pipeline2.build_dag()

        # Extract commands
        commands1 = [job.shell for job in list(dag1.waiting_jobs.keys()) + list(dag1.ready_jobs.keys())]
        commands2 = [job.shell for job in list(dag2.waiting_jobs.keys()) + list(dag2.ready_jobs.keys())]

        # Metrics pipeline should have more algorithms
        metrics_algos1 = [cmd for cmd in commands1 if "--algo" in cmd and any(
            algo in cmd for algo in ["InsertSizeMetricAlgo", "MeanQualityByCycle", "GCBias"]
        )]
        metrics_algos2 = [cmd for cmd in commands2 if "--algo" in cmd and any(
            algo in cmd for algo in ["InsertSizeMetricAlgo", "MeanQualityByCycle", "GCBias"]
        )]

        assert len(metrics_algos1) > len(metrics_algos2), "Metrics enabled should have more metric algorithms"


class TestPipelineErrorHandling:
    """Test pipeline error handling and validation"""

    def setup_method(self):
        """Setup test environment"""
        self.helper = setup_basic_test_environment()

    def teardown_method(self):
        """Cleanup test environment"""
        teardown_test_environment(self.helper)

    def test_missing_required_files(self):
        """Test handling of missing required files"""
        pipeline = self.helper.create_dnascope_pipeline(use_bam=True)

        # Remove a required file
        pipeline.reference.unlink()

        # Should handle missing file gracefully during validation
        # (In practice, this might be caught by argparse or during execution)
        try:
            pipeline.validate()
            pipeline.configure()
            dag = pipeline.build_dag()

            # Commands should still be generated, but execution would fail
            all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
            assert len(all_jobs) > 0
        except Exception:
            # Some pipelines might catch this during validation
            pass

    def test_invalid_parameter_combinations(self):
        """Test handling of invalid parameter combinations"""
        pipeline = self.helper.create_longread_pipeline()

        # Set invalid combination: haploid bed without diploid bed
        pipeline.haploid_bed = self.helper.fs.create_bed_file("haploid.bed")
        pipeline.bed = None
        pipeline.skip_small_variants = False

        with pytest.raises(SystemExit):
            pipeline.validate()

    @patch('sentieon_cli.util.check_version')
    def test_version_check_failure(self, mock_check_version):
        """Test behavior when version checks fail"""
        mock_check_version.return_value = False

        pipeline = self.helper.create_dnascope_pipeline(use_bam=True)
        pipeline.skip_version_check = False  # Enable version checking

        with pytest.raises(SystemExit):
            pipeline.validate()
            pipeline.configure()
            pipeline.build_dag()


class TestDAGDependencyValidation:
    """Test DAG dependency validation"""

    def setup_method(self):
        """Setup test environment"""
        self.helper = setup_basic_test_environment()

    def teardown_method(self):
        """Cleanup test environment"""
        teardown_test_environment(self.helper)

    def test_alignment_to_dedup_dependency(self):
        """Test that deduplication depends on alignment"""
        pipeline = self.helper.create_dnascope_pipeline()  # Uses FASTQ (requires alignment)
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Dedup should depend on LC and LC on alignment
        alignment_jobs = DAGAnalyzer.get_jobs_by_type(dag, "align")
        lc_jobs = DAGAnalyzer.get_jobs_by_type(dag, "locuscollector")
        dedup_jobs = DAGAnalyzer.get_jobs_by_type(dag, "dedup")

        lc_alignment_dep = False
        dedup_lc_dep = False
        for lc_job in lc_jobs:
            if lc_job in dag.waiting_jobs:
                deps = dag.waiting_jobs[lc_job]
                if any(align_job in deps for align_job in alignment_jobs):
                    lc_alignment_dep = True
                    break
        for dedup_job in dedup_jobs:
            if dedup_job in dag.waiting_jobs:
                deps = dag.waiting_jobs[dedup_job]
                if any(lc_job in deps for lc_job in lc_jobs):
                    dedup_lc_dep = True

        assert lc_alignment_dep, "LocusCollector should depend on alignment"
        assert dedup_lc_dep, "Dedup should depend on LocusCollector"

    def test_variant_calling_dependencies(self):
        """Test that variant calling has proper dependencies"""
        pipeline = self.helper.create_dnascope_pipeline(use_bam=True)
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Variant calling should depend on either dedup or input processing
        variant_jobs = DAGAnalyzer.get_jobs_by_type(dag, "variant")

        assert len(variant_jobs) > 0, "Should have variant calling jobs"

        for variant_job in variant_jobs:
            if variant_job in dag.waiting_jobs:
                # Should have some dependencies
                assert len(dag.waiting_jobs[variant_job]) > 0, "Variant calling should have dependencies"

    def test_cleanup_jobs_at_end(self):
        """Test that cleanup jobs are scheduled after main processing"""
        pipeline = self.helper.create_dnascope_pipeline(use_bam=True)
        pipeline.validate()
        pipeline.configure()

        dag = pipeline.build_dag()

        # Find cleanup jobs (typically have fail_ok=True)
        cleanup_jobs = [job for job in list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
                       if job.fail_ok]

        if cleanup_jobs:
            # Cleanup jobs should depend on main processing jobs
            for cleanup_job in cleanup_jobs:
                if cleanup_job in dag.waiting_jobs:
                    deps = dag.waiting_jobs[cleanup_job]
                    assert len(deps) > 0, "Cleanup jobs should have dependencies"


if __name__ == "__main__":
    pytest.main([__file__])
