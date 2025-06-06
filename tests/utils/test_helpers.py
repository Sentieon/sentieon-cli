"""
Test helper utilities for sentieon-cli testing
"""

import pathlib
import tempfile
import shutil
import argparse
from typing import Dict, List, Optional
from unittest.mock import patch
import sys
import os

# Add the parent directory to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline


def create_mock_args(loglevel="WARNING"):
    """Create a mock argparse.Namespace for pipeline testing"""
    args = argparse.Namespace()
    args.loglevel = loglevel
    return args


class MockFileSystem:
    """Helper class for creating mock file systems for testing"""

    def __init__(self):
        self.temp_dir = tempfile.mkdtemp()
        self.files = {}

    def create_file(self, path: str, content: bytes = b"", size: Optional[int] = None) -> pathlib.Path:
        """Create a mock file with optional content or specific size"""
        file_path = pathlib.Path(self.temp_dir) / path
        file_path.parent.mkdir(parents=True, exist_ok=True)

        if size is not None:
            content = b"x" * size

        file_path.write_bytes(content)
        self.files[path] = file_path
        return file_path

    def create_fastq(self, path: str, num_reads: int = 100) -> pathlib.Path:
        """Create a mock FASTQ file with specified number of reads"""
        content = ""
        for i in range(num_reads):
            content += f"@read_{i}\n"
            content += "ATCGATCGATCGATCG\n"
            content += "+\n"
            content += "IIIIIIIIIIIIIIII\n"

        return self.create_file(path, content.encode())

    def create_bam_header(self, path: str, readgroups: Optional[List[Dict[str, str]]] = None) -> pathlib.Path:
        """Create a mock BAM file"""
        if readgroups is None:
            readgroups = [{"ID": "rg1", "SM": "sample1", "PL": "ILLUMINA"}]

        # Create minimal BAM content (not a real BAM, just for testing)
        content = b"BAM_HEADER"
        for rg in readgroups:
            rg_line = "\t".join([f"{k}:{v}" for k, v in rg.items()])
            content += f"@RG\t{rg_line}\n".encode()

        return self.create_file(path, content)

    def create_model_bundle(self, path: str, pipeline_type: str = "DNAscope",
                           longread_tech: str = "HiFi", shortread_tech: str = "Illumina") -> pathlib.Path:
        """Create a mock model bundle file (ar archive will be mocked)"""
        # Create an empty file - the ar_load function will be mocked
        bundle_path = pathlib.Path(self.temp_dir) / path
        bundle_path.parent.mkdir(parents=True, exist_ok=True)

        # Just create an empty file - the actual content will be provided by mocked ar_load
        bundle_path.write_bytes(b"mock_ar_archive")

        self.files[path] = bundle_path
        return bundle_path

    def create_bed_file(self, path: str, regions: Optional[List[tuple]] = None) -> pathlib.Path:
        """Create a mock BED file"""
        if regions is None:
            regions = [("chr1", 1000, 2000), ("chr2", 5000, 6000)]

        content = ""
        for chrom, start, end in regions:
            content += f"{chrom}\t{start}\t{end}\n"

        return self.create_file(path, content.encode())

    def cleanup(self):
        """Clean up temporary files"""
        shutil.rmtree(self.temp_dir, ignore_errors=True)


class PipelineTestHelper:
    """Helper class for creating and configuring test pipelines"""

    def __init__(self):
        self.fs = MockFileSystem()

    def create_dnascope_pipeline(self, **kwargs) -> DNAscopePipeline:
        """Create a DNAscope pipeline with sensible defaults"""
        pipeline = DNAscopePipeline()

        # Setup logging first
        args = create_mock_args()
        pipeline.setup_logging(args)

        # Set default files
        pipeline.output_vcf = self.fs.create_file("output.vcf.gz")
        pipeline.reference = self.fs.create_file("reference.fa")
        pipeline.model_bundle = self.fs.create_model_bundle("model.bundle")

        # Create reference index
        self.fs.create_file("reference.fa.fai")

        # Default configuration
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = pathlib.Path(self.fs.temp_dir) / "tmp"
        pipeline.tmp_dir.mkdir(exist_ok=True)

        # Input files
        if "use_bam" in kwargs and kwargs["use_bam"]:
            pipeline.sample_input = [self.fs.create_bam_header("sample.bam")]
            pipeline.r1_fastq = []
            pipeline.r2_fastq = []
            pipeline.readgroups = []
        else:
            pipeline.sample_input = []
            pipeline.r1_fastq = [self.fs.create_fastq("sample_R1.fastq.gz")]
            pipeline.r2_fastq = [self.fs.create_fastq("sample_R2.fastq.gz")]
            pipeline.readgroups = ["@RG\\tID:rg1\\tSM:sample\\tPL:ILLUMINA"]

        # Override with any provided kwargs
        for key, value in kwargs.items():
            if key != "use_bam" and hasattr(pipeline, key):
                setattr(pipeline, key, value)

        return pipeline

    def create_longread_pipeline(self, **kwargs) -> DNAscopeLRPipeline:
        """Create a DNAscope LongRead pipeline with sensible defaults"""
        pipeline = DNAscopeLRPipeline()

        # Setup logging first
        args = create_mock_args()
        pipeline.setup_logging(args)

        # Set default files
        pipeline.output_vcf = self.fs.create_file("output.vcf.gz")
        pipeline.reference = self.fs.create_file("reference.fa")
        pipeline.model_bundle = self.fs.create_model_bundle("model.bundle")

        # Default configuration
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = pathlib.Path(self.fs.temp_dir) / "tmp"
        pipeline.tmp_dir.mkdir(exist_ok=True)

        # Input files
        if "use_fastq" in kwargs and kwargs["use_fastq"]:
            pipeline.sample_input = []
            pipeline.fastq = [self.fs.create_fastq("longread.fastq.gz", 50)]
            pipeline.readgroups = ["@RG\\tID:lr1\\tSM:sample\\tPL:PACBIO"]
        else:
            pipeline.sample_input = [self.fs.create_bam_header("longread.bam")]
            pipeline.fastq = []
            pipeline.readgroups = []

        # Override with any provided kwargs
        for key, value in kwargs.items():
            if key != "use_fastq" and hasattr(pipeline, key):
                setattr(pipeline, key, value)

        return pipeline

    def create_hybrid_pipeline(self, **kwargs) -> DNAscopeHybridPipeline:
        """Create a DNAscope Hybrid pipeline with sensible defaults"""
        pipeline = DNAscopeHybridPipeline()

        # Setup logging first
        args = create_mock_args()
        pipeline.setup_logging(args)

        # Set default files
        pipeline.output_vcf = self.fs.create_file("output.vcf.gz")
        pipeline.reference = self.fs.create_file("reference.fa")
        pipeline.model_bundle = self.fs.create_model_bundle("model.bundle", "DNAscope Hybrid")

        # Default configuration
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = pathlib.Path(self.fs.temp_dir) / "tmp"
        pipeline.tmp_dir.mkdir(exist_ok=True)

        # Input files
        pipeline.lr_aln = [self.fs.create_bam_header("longread.bam",
                                                    [{"ID": "lr1", "SM": "sample", "PL": "PACBIO"}])]
        pipeline.sr_aln = [self.fs.create_bam_header("shortread.bam",
                                                    [{"ID": "sr1", "SM": "sample", "PL": "ILLUMINA"}])]
        pipeline.sr_r1_fastq = []
        pipeline.sr_r2_fastq = []
        pipeline.sr_readgroups = []

        # Override with any provided kwargs
        for key, value in kwargs.items():
            if hasattr(pipeline, key):
                setattr(pipeline, key, value)

        return pipeline

    def cleanup(self):
        """Clean up test resources"""
        self.fs.cleanup()


class MockExternalTools:
    """Mock external tool dependencies for testing"""

    @staticmethod
    def mock_version_check():
        """Mock version checking to always pass"""
        return patch('sentieon_cli.util.check_version', return_value=True)

    @staticmethod
    def mock_library_preloaded():
        """Mock library preloading check"""
        return patch('sentieon_cli.util.library_preloaded', return_value=True)

    @staticmethod
    def mock_ar_load(bundle_info=None, bundle_members=None):
        """Mock ar_load function for bundle file testing"""
        if bundle_info is None:
            bundle_info = {
                "longReadPlatform": "HiFi",
                "shortReadPlatform": "Illumina",
                "minScriptVersion": "1.0.0",
                "pipeline": "DNAscope Hybrid"
            }

        if bundle_members is None:
            bundle_members = [
                "longreadsv.model",
                "cnv.model",
                "bwa.model",
                "dnascope.model",
                "diploid_model",
                "haploid_model",
                "gvcf_model",
                "hybrid.model",
                "HybridStage1.model",
                "HybridStage2.model",
                "HybridStage3.model"
            ]

        def mock_ar_load_func(path):
            # Parse the path to determine what to return
            if path.endswith("/bundle_info.json"):
                import json
                return json.dumps(bundle_info).encode()
            elif "/" in path and not path.endswith("/bundle_info.json"):
                # Specific member request (path/member)
                base_path, member = path.rsplit("/", 1)
                if member in bundle_members:
                    return b"mock_model_data"
                else:
                    raise RuntimeError(f"Member {member} not found in archive")
            else:
                # List all members - called with just the bundle path
                return bundle_members

        return patch('sentieon_cli.archive.ar_load', side_effect=mock_ar_load_func)


class CommandValidator:
    """Helper for validating generated commands"""

    @staticmethod
    def validate_sentieon_command(command: str) -> bool:
        """Validate that a sentieon command is well-formed"""
        if not command.startswith("sentieon"):
            return False

        parts = command.split()
        if len(parts) < 2:
            return False

        # Should have driver or util as second argument
        if parts[1] not in ["driver", "util"]:
            return False

        return True


class DAGAnalyzer:
    """Helper for analyzing DAG structure in tests"""

    @staticmethod
    def get_job_by_name(dag, name_pattern: str):
        """Find a job in the DAG by name pattern"""
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())

        for job in all_jobs:
            if name_pattern in job.name:
                return job
        return None

    @staticmethod
    def get_jobs_by_type(dag, job_type: str):
        """Get all jobs of a specific type"""
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())

        matching_jobs = []
        for job in all_jobs:
            if job_type in job.name:
                matching_jobs.append(job)

        return matching_jobs

    @staticmethod
    def count_jobs_by_thread_usage(dag) -> Dict[int, int]:
        """Count jobs by their thread usage"""
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())

        thread_counts = {}
        for job in all_jobs:
            threads = job.threads
            thread_counts[threads] = thread_counts.get(threads, 0) + 1

        return thread_counts


# Convenience function for common test setup
def setup_basic_test_environment():
    """Setup a basic test environment with common mocks"""
    helper = PipelineTestHelper()

    patches = [
        MockExternalTools.mock_version_check(),
        MockExternalTools.mock_library_preloaded(),
        MockExternalTools.mock_ar_load(),
    ]

    # Start all patches
    started_patches = [p.start() for p in patches]

    return helper, started_patches


def teardown_test_environment(helper: PipelineTestHelper, patches: List):
    """Teardown test environment"""
    # Stop all patches
    for patch in patches:
        patch.stop()

    # Cleanup filesystem
    helper.cleanup()
