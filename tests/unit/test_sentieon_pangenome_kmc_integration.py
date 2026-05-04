"""
Tests for KMC patch integration in SentieonPangenome
"""

import sys
import unittest
from unittest.mock import MagicMock, patch

from sentieon_cli.sentieon_pangenome import SentieonPangenome


class TestSentieonPangenomeKMCIntegration(unittest.TestCase):
    def test_validate_kmc_patch_success(self):
        """Test validation passes when KMC patch check succeeds"""
        pipeline = SentieonPangenome()
        pipeline.sample_input = ["sample.bam"]
        pipeline.skip_version_check = False
        pipeline.skip_pangenome_name_checks = True # skip irrelevant checks
        pipeline.skip_contig_checks = True
        pipeline.bed = "bed" # mock existence
        pipeline.reference = "ref"
        pipeline.pop_vcf = "pop"
        pipeline.gbz = "ref.grch38.gbz"
        pipeline.fai_data = {}
        pipeline.pop_vcf_contigs = {}
        
        # Mock dependencies
        pipeline.handle_arguments = MagicMock()
        pipeline.setup_logging = MagicMock()
        pipeline.validate_bundle = MagicMock()
        pipeline.validate_fastq_rg = MagicMock()
        pipeline.validate_output_vcf = MagicMock()
        pipeline.validate_ref = MagicMock()
        pipeline.collect_readgroups = MagicMock()
        
        with patch("sentieon_cli.sentieon_pangenome.check_version", return_value=True), \
             patch("sentieon_cli.sentieon_pangenome.check_kmc_patch", return_value=True) as mock_check, \
             patch("sys.exit") as mock_exit:
            
            pipeline.validate()
            
            mock_check.assert_called_with("kmc")
            mock_exit.assert_not_called()

    def test_validate_kmc_patch_failure(self):
        """Test validation fails when KMC patch check fails"""
        pipeline = SentieonPangenome()
        pipeline.sample_input = ["sample.bam"]
        pipeline.skip_version_check = False
        # other attributes mock
        pipeline.validate_bundle = MagicMock()
        pipeline.validate_fastq_rg = MagicMock()
        pipeline.validate_output_vcf = MagicMock()
        pipeline.validate_ref = MagicMock()
        pipeline.collect_readgroups = MagicMock()
        pipeline.bed = "bed"
        pipeline.fai_data = {}
        pipeline.pop_vcf_contigs = {}

        pipeline.logger = MagicMock()

        with patch("sentieon_cli.sentieon_pangenome.check_version", return_value=True), \
             patch("sentieon_cli.sentieon_pangenome.check_kmc_patch", return_value=False) as mock_check, \
             patch("sys.exit") as mock_exit:
            
            pipeline.validate()
            
            mock_check.assert_called_with("kmc")
            pipeline.logger.error.assert_called()
            mock_exit.assert_called_with(2)
