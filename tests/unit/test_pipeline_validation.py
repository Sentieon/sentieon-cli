"""
Unit tests for pipeline validation logic
"""

import pathlib
import pytest
import tempfile
import sys
import os
import json
from unittest.mock import patch

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline
from tests.utils.test_helpers import create_mock_args


class TestDNAscopePipelineValidation:
    """Test DNAscope pipeline validation logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.pipeline = DNAscopePipeline()

        # Setup logging first
        args = create_mock_args()
        self.pipeline.setup_logging(args)

        self.temp_dir = tempfile.mkdtemp()

        # Create mock files
        self.mock_vcf = pathlib.Path(self.temp_dir) / "output.vcf.gz"
        self.mock_ref = pathlib.Path(self.temp_dir) / "reference.fa"
        self.mock_bam = pathlib.Path(self.temp_dir) / "sample.bam"
        self.mock_fastq = pathlib.Path(self.temp_dir) / "sample_R1.fastq.gz"
        self.mock_bundle = pathlib.Path(self.temp_dir) / "model.bundle"

        # Create empty files
        for file_path in [self.mock_ref, self.mock_bam, self.mock_fastq, self.mock_bundle]:
            file_path.touch()

    def test_valid_configuration_with_bam_input(self):
        """Test valid configuration with BAM input"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.sr_r1_fastq = []
        self.pipeline.sr_readgroups = []

        # Should not raise any exceptions
        try:
            self.pipeline.validate()
        except SystemExit:
            pytest.fail("Valid configuration should not raise SystemExit")

    def test_valid_configuration_with_fastq_input(self):
        """Test valid configuration with FASTQ input"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = []
        self.pipeline.r1_fastq = [self.mock_fastq]
        self.pipeline.readgroups = ["@RG\\tID:test\\tSM:sample"]

        # Should not raise any exceptions
        try:
            self.pipeline.validate()
        except SystemExit:
            pytest.fail("Valid configuration should not raise SystemExit")

    def test_missing_inputs_raises_error(self):
        """Test that missing inputs raise an error"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = []
        self.pipeline.sr_r1_fastq = []
        self.pipeline.sr_readgroups = []

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    def test_invalid_output_extension_raises_error(self):
        """Test that invalid output file extension raises an error"""
        invalid_vcf = pathlib.Path(self.temp_dir) / "output.vcf"  # Missing .gz
        self.pipeline.output_vcf = invalid_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.sr_r1_fastq = []
        self.pipeline.sr_readgroups = []

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    def test_mismatched_fastq_readgroups_raises_error(self):
        """Test that mismatched FASTQ and readgroup counts raise an error"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = []
        self.pipeline.sr_r1_fastq = [self.mock_fastq]
        self.pipeline.sr_readgroups = []  # Empty readgroups with non-empty fastq

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    def test_skip_multiqc_when_skip_metrics(self):
        """Test that skip_multiqc is set when skip_metrics is True"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.sr_r1_fastq = []
        self.pipeline.sr_readgroups = []
        self.pipeline.skip_metrics = True
        self.pipeline.skip_multiqc = False

        self.pipeline.validate()
        assert self.pipeline.skip_multiqc is True


class TestDNAscopeLRPipelineValidation:
    """Test DNAscope LongRead pipeline validation logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.pipeline = DNAscopeLRPipeline()

        # Setup logging first
        args = create_mock_args()
        self.pipeline.setup_logging(args)

        self.temp_dir = tempfile.mkdtemp()

        # Create mock files
        self.mock_vcf = pathlib.Path(self.temp_dir) / "output.vcf.gz"
        self.mock_ref = pathlib.Path(self.temp_dir) / "reference.fa"
        self.mock_bam = pathlib.Path(self.temp_dir) / "sample.bam"
        self.mock_fastq = pathlib.Path(self.temp_dir) / "sample.fastq.gz"
        self.mock_bundle = pathlib.Path(self.temp_dir) / "model.bundle"
        self.mock_bed = pathlib.Path(self.temp_dir) / "regions.bed"

        # Create empty files
        for file_path in [self.mock_ref, self.mock_bam, self.mock_fastq, self.mock_bundle, self.mock_bed]:
            file_path.touch()

    def test_valid_configuration_with_bam_input(self):
        """Test valid configuration with BAM input"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.fastq = []
        self.pipeline.readgroups = []

        try:
            self.pipeline.validate()
        except SystemExit:
            pytest.fail("Valid configuration should not raise SystemExit")

    def test_valid_configuration_with_fastq_input(self):
        """Test valid configuration with FASTQ input"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.lr_aln = []
        self.pipeline.fastq = [self.mock_fastq]
        self.pipeline.readgroups = ["@RG\\tID:test\\tSM:sample"]

        try:
            self.pipeline.validate()
        except SystemExit:
            pytest.fail("Valid configuration should not raise SystemExit")

    def test_ont_technology_skips_cnv(self):
        """Test that ONT technology automatically skips CNV calling"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.fastq = []
        self.pipeline.readgroups = []
        self.pipeline.tech = "ONT"
        self.pipeline.skip_cnv = False

        self.pipeline.validate()
        assert self.pipeline.skip_cnv is True

    def test_haploid_bed_requires_diploid_bed(self):
        """Test that haploid bed requires diploid bed"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.lr_aln = [self.mock_bam]
        self.pipeline.fastq = []
        self.pipeline.readgroups = []
        self.pipeline.haploid_bed = self.mock_bed
        self.pipeline.bed = None  # No diploid bed
        self.pipeline.skip_small_variants = False

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    def test_mismatched_fastq_readgroups_raises_error(self):
        """Test that mismatched FASTQ and readgroup counts raise an error"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.lr_aln = []
        self.pipeline.fastq = [self.mock_fastq]
        self.pipeline.readgroups = []  # Empty readgroups with non-empty fastq

        with pytest.raises(SystemExit):
            self.pipeline.validate()


class TestDNAscopeHybridPipelineValidation:
    """Test DNAscope Hybrid pipeline validation logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.pipeline = DNAscopeHybridPipeline()

        # Setup logging first
        args = create_mock_args()
        self.pipeline.setup_logging(args)

        self.temp_dir = tempfile.mkdtemp()

        # Create mock files
        self.mock_vcf = pathlib.Path(self.temp_dir) / "output.vcf.gz"
        self.mock_ref = pathlib.Path(self.temp_dir) / "reference.fa"
        self.mock_fai = pathlib.Path(self.temp_dir) / "reference.fa.fai"
        self.mock_lr_bam = pathlib.Path(self.temp_dir) / "longread.bam"
        self.mock_sr_bam = pathlib.Path(self.temp_dir) / "shortread.bam"
        self.mock_fastq = pathlib.Path(self.temp_dir) / "sample_R1.fastq.gz"

        # Create empty files
        for file_path in [self.mock_ref, self.mock_fai, self.mock_lr_bam, self.mock_sr_bam, self.mock_fastq]:
            file_path.touch()

        # Create mock bundle file (content will be mocked by ar_load)
        self.mock_bundle = pathlib.Path(self.temp_dir) / "model.bundle"
        self.mock_bundle.write_bytes(b"mock_ar_archive")

    @patch('sentieon_cli.dnascope_hybrid.ar_load')
    def test_valid_hybrid_configuration(self, mock_ar_load):
        """Test valid hybrid pipeline configuration"""
        # Setup ar_load mock
        bundle_info = {
            "longReadPlatform": "HiFi",
            "shortReadPlatform": "Illumina",
            "minScriptVersion": "1.0.0",
            "pipeline": "DNAscope Hybrid"
        }
        mock_ar_load.side_effect = [
            json.dumps(bundle_info).encode(),  # bundle_info.json
            ["longreadsv.model", "cnv.model", "bwa.model"]  # bundle members
        ]

        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.lr_aln = [self.mock_lr_bam]
        self.pipeline.sr_aln = [self.mock_sr_bam]
        self.pipeline.sr_r1_fastq = []
        self.pipeline.sr_readgroups = []

        try:
            # Only test bundle validation since full validation requires more mocking
            self.pipeline.validate_bundle()
        except SystemExit:
            pytest.fail("Valid bundle should not raise SystemExit")

    @patch('sentieon_cli.dnascope_hybrid.ar_load')
    @patch('sentieon_cli.command_strings.get_rg_lines')
    def test_missing_short_read_input_raises_error(self, mock_get_rg, mock_ar_load):
        """Test that missing short-read input raises an error"""
        # Setup ar_load mock
        bundle_info = {
            "longReadPlatform": "HiFi",
            "shortReadPlatform": "Illumina",
            "minScriptVersion": "1.0.0",
            "pipeline": "DNAscope Hybrid"
        }
        mock_ar_load.side_effect = [
            json.dumps(bundle_info).encode(),  # bundle_info.json
            ["longreadsv.model", "cnv.model", "bwa.model"]  # bundle members
        ]

        # Mock readgroup lines
        mock_get_rg.return_value = ["@RG\tID:lr1\tSM:sample"]

        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.lr_aln = [self.mock_lr_bam]
        self.pipeline.sr_aln = []  # No short-read input
        self.pipeline.sr_r1_fastq = []
        self.pipeline.sr_readgroups = []

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    @patch('sentieon_cli.dnascope_hybrid.ar_load')
    def test_invalid_bundle_pipeline_raises_error(self, mock_ar_load):
        """Test that invalid bundle pipeline type raises an error"""
        # Create bundle with wrong pipeline type
        bundle_info = {
            "longReadPlatform": "HiFi",
            "shortReadPlatform": "Illumina",
            "minScriptVersion": "1.0.0",
            "pipeline": "DNAscope"  # Wrong pipeline
        }
        # Mock both calls: first for bundle_info.json, second for member list
        mock_ar_load.side_effect = [
            json.dumps(bundle_info).encode(),  # bundle_info.json
            ["longreadsv.model", "cnv.model", "bwa.model"]  # bundle members
        ]

        self.pipeline.model_bundle = self.mock_bundle

        with pytest.raises(SystemExit):
            self.pipeline.validate_bundle()


class TestPipelineConfigurationHelpers:
    """Test helper methods used in pipeline configuration"""

    def test_total_input_size_calculation(self):
        """Test calculation of total input file size"""
        pipeline = DNAscopePipeline()

        # Setup logging first
        args = create_mock_args()
        pipeline.setup_logging(args)

        temp_dir = tempfile.mkdtemp()

        # Create test files with known sizes
        bam_file = pathlib.Path(temp_dir) / "test.bam"
        fastq_file = pathlib.Path(temp_dir) / "test.fastq.gz"

        bam_file.write_bytes(b"x" * 1000)  # 1KB
        fastq_file.write_bytes(b"y" * 500)  # 0.5KB

        pipeline.sample_input = [bam_file]
        pipeline.sr_r1_fastq = [fastq_file]
        pipeline.sr_r2_fastq = []

        total_size = pipeline.total_input_size()
        assert total_size == 1500  # 1KB + 0.5KB

    def test_memory_calculation_logic(self):
        """Test bwt_max_mem calculation logic"""
        pipeline = DNAscopePipeline()

        # Setup logging first
        args = create_mock_args()
        pipeline.setup_logging(args)

        # Test with specific memory values
        pipeline.bwt_max_mem = "10G"
        pipeline.set_bwt_max_mem(0, 1)

        import os
        assert os.environ.get("bwt_max_mem") == "10G"


if __name__ == "__main__":
    pytest.main([__file__])
