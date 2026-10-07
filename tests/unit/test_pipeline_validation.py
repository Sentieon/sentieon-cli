"""
Unit tests for pipeline validation logic
"""

import pathlib
import pytest
import tempfile
import sys
import os
import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

# Add the parent directory to the path so we can import sentieon_cli
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline
from sentieon_cli.hybrid_pangenome import HybridPangenome
from sentieon_cli.pangenome_meta import (
    PangenomeMetadataError,
    PangenomeReference,
)
from sentieon_cli.sentieon_pangenome import SentieonPangenome
from sentieon_cli.util import VcfHeaderError, set_bwt_max_mem
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
        for file_path in [self.mock_ref, self.mock_bam, self.mock_fastq]:
            file_path.touch()

        # An empty ar archive: a bundle without a bundle_info.json
        self.mock_bundle.write_bytes(b"!<arch>\n")

        # Create BWA index files next to the reference
        for suf in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            pathlib.Path(str(self.mock_ref) + suf).touch()

    def test_valid_configuration_with_bam_input(self):
        """Test valid configuration with BAM input"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.r1_fastq = []
        self.pipeline.readgroups = []

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
        self.pipeline.r1_fastq = []
        self.pipeline.readgroups = []

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    def test_invalid_output_extension_raises_error(self):
        """Test that invalid output file extension raises an error"""
        invalid_vcf = pathlib.Path(self.temp_dir) / "output.vcf"  # Missing .gz
        self.pipeline.output_vcf = invalid_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.r1_fastq = []
        self.pipeline.readgroups = []

        with pytest.raises(SystemExit):
            self.pipeline.validate()

    def test_mismatched_fastq_readgroups_raises_error(self):
        """Test that mismatched FASTQ and readgroup counts raise an error"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = []
        # Two fastq files, but only one readgroup. Both must be non-empty to
        # get past the "supply --sample_input or --r1_fastq" guard.
        self.pipeline.r1_fastq = [self.mock_fastq, self.mock_fastq]
        self.pipeline.readgroups = ["@RG\\tID:a\\tSM:s"]

        with patch.object(self.pipeline.logger, "error") as mock_error:
            with pytest.raises(SystemExit) as excinfo:
                self.pipeline.validate()

        assert excinfo.value.code == 2
        mock_error.assert_any_call(
            "The number of readgroups does not equal the number of fastq files"
        )

    def test_skip_multiqc_when_skip_metrics(self):
        """Test that skip_multiqc is set when skip_metrics is True"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.r1_fastq = []
        self.pipeline.readgroups = []
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
        self.mock_fai = pathlib.Path(self.temp_dir) / "reference.fa.fai"
        self.mock_bam = pathlib.Path(self.temp_dir) / "sample.bam"
        self.mock_fastq = pathlib.Path(self.temp_dir) / "sample.fastq.gz"
        self.mock_bundle = pathlib.Path(self.temp_dir) / "model.bundle"
        self.mock_bed = pathlib.Path(self.temp_dir) / "regions.bed"

        # Create empty files
        for file_path in [self.mock_ref, self.mock_fai, self.mock_bam, self.mock_fastq, self.mock_bundle, self.mock_bed]:
            file_path.touch()

    def test_valid_configuration_with_bam_input(self):
        """Test valid configuration with BAM input"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.fastq = []
        self.pipeline.readgroups = []

        # Mock archive loading
        with patch('sentieon_cli.dnascope_longread.ar_load') as mock_ar_load:
            # Mock bundle info
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
        self.pipeline.fastq = [self.mock_fastq]
        self.pipeline.readgroups = ["@RG\\tID:test\\tSM:sample"]

        # Mock archive loading
        with patch('sentieon_cli.dnascope_longread.ar_load') as mock_ar_load:
            # Mock bundle info
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

        # Mock archive loading
        with patch('sentieon_cli.dnascope_longread.ar_load') as mock_ar_load:
            # Mock bundle info
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

            self.pipeline.validate()
            assert self.pipeline.skip_cnv is True

    def test_haploid_bed_requires_diploid_bed(self):
        """Test that haploid bed requires diploid bed"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = [self.mock_bam]
        self.pipeline.fastq = []
        self.pipeline.readgroups = []
        self.pipeline.haploid_bed = self.mock_bed
        self.pipeline.bed = None  # No diploid bed
        self.pipeline.skip_small_variants = False

        with patch.object(self.pipeline.logger, "error") as mock_error:
            with pytest.raises(SystemExit) as excinfo:
                self.pipeline.validate()

        assert excinfo.value.code == 2
        mock_error.assert_any_call(
            "Please supply a BED file of diploid regions to distinguish "
            "haploid and diploid regions of the genome."
        )

    def test_mismatched_fastq_readgroups_raises_error(self):
        """Test that mismatched FASTQ and readgroup counts raise an error"""
        self.pipeline.output_vcf = self.mock_vcf
        self.pipeline.reference = self.mock_ref
        self.pipeline.model_bundle = self.mock_bundle
        self.pipeline.sample_input = []
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
        pipeline.r1_fastq = [fastq_file]
        pipeline.r2_fastq = []

        total_size = pipeline.total_input_size()
        assert total_size == 1500  # 1KB + 0.5KB

    def test_bwt_max_mem_override(self, monkeypatch):
        """`--bwt_max_mem` short-circuits the calculation"""
        monkeypatch.delenv("bwt_max_mem", raising=False)

        assert set_bwt_max_mem(0, 1, override="10G") == "10G"
        assert os.environ.get("bwt_max_mem") == "10G"

    @pytest.mark.parametrize(
        "total_input_size,n_alignment_jobs,expected",
        [
            (0, 1, "54G"),  # 64 - 4 - 0 = 60; 60 / 1 - 6 = 54
            (0, 2, "24G"),  # 60 / 2 - 6 = 24
            (30 * 1024**3, 1, "0G"),  # 60 - 30 * 2.3 < 0, floored at 0
        ],
    )
    def test_bwt_max_mem_calculation(
        self, monkeypatch, total_input_size, n_alignment_jobs, expected
    ):
        """bwt_max_mem is derived from the memory and the input size"""
        monkeypatch.setenv("bwt_max_mem", "unset")

        with patch(
            "sentieon_cli.util.total_memory", return_value=64 * 1024**3
        ):
            result = set_bwt_max_mem(total_input_size, n_alignment_jobs)

        assert result == expected
        assert os.environ.get("bwt_max_mem") == expected


class TestValidateBwaIndex:
    """Test BasePipeline.validate_bwa_index()"""

    BWA_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

    def setup_method(self):
        self.pipeline = DNAscopePipeline()
        args = create_mock_args()
        self.pipeline.setup_logging(args)

        self.temp_dir = tempfile.mkdtemp()
        self.mock_ref = pathlib.Path(self.temp_dir) / "reference.fa"
        self.mock_ref.touch()
        self.pipeline.reference = self.mock_ref

    def _create_indexes(self, suffixes):
        for suf in suffixes:
            pathlib.Path(str(self.mock_ref) + suf).touch()

    def test_all_index_files_present(self):
        self._create_indexes(self.BWA_SUFFIXES)
        # Should not raise
        self.pipeline.validate_bwa_index()

    def test_all_index_files_missing_raises(self):
        with pytest.raises(SystemExit) as excinfo:
            self.pipeline.validate_bwa_index()
        assert excinfo.value.code == 2

    def test_partial_index_files_missing_raises(self):
        # All but .sa present
        self._create_indexes(self.BWA_SUFFIXES[:-1])
        with pytest.raises(SystemExit) as excinfo:
            self.pipeline.validate_bwa_index()
        assert excinfo.value.code == 2


PANGENOME_PIPELINES = [SentieonPangenome, HybridPangenome]

GRCH38_REFERENCE = PangenomeReference(
    ref_name="GRCh38",
    contig_prefix="GRCh38#0#",
    contigs=["chr1", "chr2"],
    reference_samples=["GRCh38", "CHM13"],
    n_paths=2,
)
CHM13_REFERENCE = PangenomeReference(
    ref_name="CHM13",
    contig_prefix="CHM13#0#",
    contigs=["chr1", "chr2"],
    reference_samples=["CHM13", "GRCh38"],
    n_paths=2,
)


@pytest.mark.parametrize(
    "cls", PANGENOME_PIPELINES, ids=lambda c: c.__name__
)
class TestResolvePangenomeReference:
    """`BasePangenome.resolve_pangenome_reference` in both pipelines.

    The reference name and contig prefix come from the graph's own
    metadata, so a graph built against the wrong reference cannot reach
    `vg haplotypes --set-reference`, which would silently drop every
    reference path.
    """

    def _pipeline(self, cls, tmp_path, **attributes):
        pipeline = cls()
        pipeline.logger = MagicMock()
        pipeline.gbz = tmp_path / "graph.gbz"
        pipeline.hapl = tmp_path / "graph.hapl"
        pipeline.gbz.touch()
        pipeline.hapl.touch()
        pipeline.fai_data = {"chr1": {"length": 10}, "chr2": {"length": 10}}
        for key, value in attributes.items():
            setattr(pipeline, key, value)
        return pipeline

    def _detect(self, reference):
        return patch(
            "sentieon_cli.base_pangenome.detect_pangenome_reference",
            return_value=reference,
        )

    def _build(self, build):
        return patch(
            "sentieon_cli.base_pangenome.detect_reference_build",
            return_value=build,
        )

    def _hapl(self, chains=2):
        return patch(
            "sentieon_cli.base_pangenome.read_hapl_header",
            return_value=SimpleNamespace(top_level_chains=chains),
        )

    def test_the_detected_reference_populates_both_attributes(
        self, cls, tmp_path
    ):
        pipeline = self._pipeline(cls, tmp_path)
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "GRCh38"
        assert pipeline.pangenome_contig_prefix == "GRCh38#0#"
        assert pipeline.pangenome_reference is GRCH38_REFERENCE

    def test_a_chm13_graph_resolves_to_the_chm13_prefix(
        self, cls, tmp_path
    ):
        pipeline = self._pipeline(cls, tmp_path)
        with self._detect(CHM13_REFERENCE), self._build("chm13"), (
            self._hapl()
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "CHM13"
        assert pipeline.pangenome_contig_prefix == "CHM13#0#"

    def test_a_matching_override_is_accepted(self, cls, tmp_path):
        pipeline = self._pipeline(
            cls,
            tmp_path,
            pangenome_ref_name="GRCh38",
            pangenome_contig_prefix="GRCh38#0#",
        )
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "GRCh38"

    def test_a_mismatched_name_override_exits(self, cls, tmp_path):
        pipeline = self._pipeline(cls, tmp_path, pangenome_ref_name="hg38")
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            with pytest.raises(SystemExit) as excinfo:
                pipeline.resolve_pangenome_reference()
        assert excinfo.value.code == 2
        pipeline.logger.error.assert_called()

    def test_a_mismatched_prefix_override_exits(self, cls, tmp_path):
        pipeline = self._pipeline(
            cls, tmp_path, pangenome_contig_prefix="GRCh38#1#"
        )
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            with pytest.raises(SystemExit) as excinfo:
                pipeline.resolve_pangenome_reference()
        assert excinfo.value.code == 2

    def test_the_skip_flag_honors_a_mismatched_override(
        self, cls, tmp_path
    ):
        pipeline = self._pipeline(
            cls,
            tmp_path,
            pangenome_ref_name="hg38",
            skip_pangenome_name_checks=True,
        )
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "hg38"
        assert pipeline.pangenome_contig_prefix == "GRCh38#0#"
        pipeline.logger.warning.assert_called()
        pipeline.logger.error.assert_not_called()

    def test_a_contig_missing_from_the_fai_exits(self, cls, tmp_path):
        pipeline = self._pipeline(cls, tmp_path)
        pipeline.fai_data = {"chr1": {"length": 10}}
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            with pytest.raises(SystemExit) as excinfo:
                pipeline.resolve_pangenome_reference()
        assert excinfo.value.code == 2
        message = pipeline.logger.error.call_args[0]
        assert "missing from the reference FASTA index" in message[0]
        assert "chr2" in message

    def test_the_skip_flag_downgrades_a_missing_contig(self, cls, tmp_path):
        pipeline = self._pipeline(
            cls, tmp_path, skip_pangenome_name_checks=True
        )
        pipeline.fai_data = {"chr1": {"length": 10}}
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            pipeline.resolve_pangenome_reference()
        pipeline.logger.warning.assert_called()

    def test_a_reference_build_mismatch_exits(self, cls, tmp_path):
        """A CHM13 graph with an hg38 reference FASTA.

        Contig names alone cannot catch this: both builds name their
        chromosomes chr1..chrM.
        """
        pipeline = self._pipeline(cls, tmp_path)
        with self._detect(CHM13_REFERENCE), self._build("hg38"), (
            self._hapl()
        ):
            with pytest.raises(SystemExit) as excinfo:
                pipeline.resolve_pangenome_reference()
        assert excinfo.value.code == 2
        assert "reference build" in pipeline.logger.error.call_args[0][0]

    def test_an_unknown_reference_build_is_only_logged(self, cls, tmp_path):
        pipeline = self._pipeline(cls, tmp_path)
        with self._detect(GRCH38_REFERENCE), self._build(None), (
            self._hapl()
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "GRCh38"
        pipeline.logger.error.assert_not_called()

    def test_a_hapl_chain_count_mismatch_only_warns(self, cls, tmp_path):
        pipeline = self._pipeline(cls, tmp_path)
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), (
            self._hapl(chains=195)
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "GRCh38"
        pipeline.logger.warning.assert_called()

    def test_an_unreadable_hapl_exits(self, cls, tmp_path):
        pipeline = self._pipeline(cls, tmp_path)
        hapl_error = PangenomeMetadataError("not a vg haplotypes file")
        with self._detect(GRCH38_REFERENCE), self._build("hg38"), patch(
            "sentieon_cli.base_pangenome.read_hapl_header",
            side_effect=hapl_error,
        ):
            with pytest.raises(SystemExit) as excinfo:
                pipeline.resolve_pangenome_reference()
        assert excinfo.value.code == 2

    def test_a_dry_run_falls_back_to_the_defaults(self, cls, tmp_path):
        """Mock graph files cannot be parsed, but dry runs still work"""
        pipeline = self._pipeline(cls, tmp_path, dry_run=True)
        error = PangenomeMetadataError("not a GBZ pangenome file")
        with patch(
            "sentieon_cli.base_pangenome.detect_pangenome_reference",
            side_effect=error,
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "GRCh38"
        assert pipeline.pangenome_contig_prefix == "GRCh38#0#"
        assert pipeline.pangenome_reference is None

    def test_a_dry_run_keeps_the_overrides(self, cls, tmp_path):
        pipeline = self._pipeline(
            cls,
            tmp_path,
            dry_run=True,
            pangenome_ref_name="CHM13",
            pangenome_contig_prefix="CHM13#0#",
        )
        error = PangenomeMetadataError("not a GBZ pangenome file")
        with patch(
            "sentieon_cli.base_pangenome.detect_pangenome_reference",
            side_effect=error,
        ):
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "CHM13"
        assert pipeline.pangenome_contig_prefix == "CHM13#0#"

    def test_an_unreadable_graph_exits_outside_a_dry_run(
        self, cls, tmp_path
    ):
        pipeline = self._pipeline(cls, tmp_path)
        error = PangenomeMetadataError("not a GBZ pangenome file")
        with patch(
            "sentieon_cli.base_pangenome.detect_pangenome_reference",
            side_effect=error,
        ):
            with pytest.raises(SystemExit) as excinfo:
                pipeline.resolve_pangenome_reference()
        assert excinfo.value.code == 2

    def test_both_overrides_survive_an_unreadable_graph(self, cls, tmp_path):
        pipeline = self._pipeline(
            cls,
            tmp_path,
            pangenome_ref_name="CHM13",
            pangenome_contig_prefix="CHM13#0#",
        )
        error = PangenomeMetadataError("not a GBZ pangenome file")
        with patch(
            "sentieon_cli.base_pangenome.detect_pangenome_reference",
            side_effect=error,
        ), self._hapl():
            pipeline.resolve_pangenome_reference()
        assert pipeline.pangenome_ref_name == "CHM13"
        assert pipeline.pangenome_contig_prefix == "CHM13#0#"
        pipeline.logger.warning.assert_called()


@pytest.mark.parametrize(
    "cls", PANGENOME_PIPELINES, ids=lambda c: c.__name__
)
class TestPangenomeContigLengthChecks:
    """`BasePangenome.validate_pangenome_contig_lengths` in both pipelines.

    The pop VCF is compared with the reference FASTA rather than with a
    table of GRCh38 lengths, so the check covers every backbone.
    """

    def _pipeline(self, cls, reference=GRCH38_REFERENCE, **attributes):
        pipeline = cls()
        pipeline.logger = MagicMock()
        pipeline.pangenome_reference = reference
        pipeline.pangenome_ref_name = (
            None if reference is None else reference.ref_name
        )
        pipeline.fai_data = {
            "chr1": {"length": 1000},
            "chr2": {"length": 2000},
        }
        pipeline.pop_vcf_contigs = {"chr1": 1000, "chr2": 2000}
        for key, value in attributes.items():
            setattr(pipeline, key, value)
        return pipeline

    def test_matching_lengths_pass(self, cls):
        pipeline = self._pipeline(cls)
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.error.assert_not_called()
        # The INFO line reports how many backbone contigs were compared
        assert pipeline.logger.info.call_args[0][1] == 2

    def test_a_mismatching_length_exits(self, cls):
        pipeline = self._pipeline(cls)
        pipeline.pop_vcf_contigs["chr2"] = 2001
        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_pangenome_contig_lengths()
        assert excinfo.value.code == 2
        reported = pipeline.logger.error.call_args[0][-1]
        assert "chr2" in reported
        assert "2000" in reported and "2001" in reported

    def test_a_contig_missing_from_the_pop_vcf_exits(self, cls):
        pipeline = self._pipeline(cls)
        del pipeline.pop_vcf_contigs["chr2"]
        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_pangenome_contig_lengths()
        assert excinfo.value.code == 2
        assert "absent" in pipeline.logger.error.call_args[0][-1]

    def test_a_contig_without_a_length_exits(self, cls):
        pipeline = self._pipeline(cls)
        pipeline.pop_vcf_contigs["chr2"] = None
        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_pangenome_contig_lengths()
        assert excinfo.value.code == 2
        assert "no length" in pipeline.logger.error.call_args[0][-1]

    def test_skip_contig_checks_bypasses_the_check(self, cls):
        pipeline = self._pipeline(cls, skip_contig_checks=True)
        pipeline.pop_vcf_contigs["chr1"] = 1
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.error.assert_not_called()

    def test_the_name_check_flag_does_not_bypass_the_check(self, cls):
        """Only `--skip_contig_checks` turns this check off"""
        pipeline = self._pipeline(cls, skip_pangenome_name_checks=True)
        pipeline.pop_vcf_contigs["chr1"] = 1
        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_pangenome_contig_lengths()
        assert excinfo.value.code == 2

    def test_a_chm13_backbone_is_checked(self, cls):
        """The old GRCh38-only check skipped every other backbone"""
        pipeline = self._pipeline(cls, reference=CHM13_REFERENCE)
        pipeline.pop_vcf_contigs["chr1"] = 248956422
        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_pangenome_contig_lengths()
        assert excinfo.value.code == 2
        assert "CHM13" in pipeline.logger.error.call_args[0]

    def test_a_matching_chm13_backbone_passes(self, cls):
        pipeline = self._pipeline(cls, reference=CHM13_REFERENCE)
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.error.assert_not_called()

    def test_empty_pop_vcf_contigs_warn_and_skip(self, cls):
        pipeline = self._pipeline(cls, pop_vcf_contigs={})
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.warning.assert_called()
        pipeline.logger.error.assert_not_called()

    def test_empty_pop_vcf_contigs_only_debug_in_a_dry_run(self, cls):
        pipeline = self._pipeline(cls, pop_vcf_contigs={}, dry_run=True)
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.warning.assert_not_called()
        pipeline.logger.debug.assert_called()

    def test_an_undetected_backbone_checks_the_shared_contigs(self, cls):
        """Without graph metadata, the fai/pop VCF intersection is used"""
        pipeline = self._pipeline(cls, reference=None)
        pipeline.pangenome_ref_name = "GRCh38"
        pipeline.pop_vcf_contigs = {"chr1": 1001, "chr3": 30}
        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_pangenome_contig_lengths()
        assert excinfo.value.code == 2
        reported = pipeline.logger.error.call_args[0][-1]
        assert "chr1" in reported and "chr2" not in reported

    def test_a_contig_missing_from_the_fai_is_not_reported_twice(self, cls):
        """`check_pangenome_reference_fasta` already reported it"""
        pipeline = self._pipeline(cls)
        del pipeline.fai_data["chr2"]
        del pipeline.pop_vcf_contigs["chr2"]
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.error.assert_not_called()
        assert pipeline.logger.info.call_args[0][1] == 1


@pytest.mark.parametrize(
    "cls", PANGENOME_PIPELINES, ids=lambda c: c.__name__
)
class TestPopVcfHeader:
    """`load_pop_vcf_header` and `check_pop_vcf_id` in both pipelines"""

    HEADER = (
        "##fileformat=VCFv4.2\n"
        "##SentieonVcfID=population-test-20260101\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    def _pipeline(self, cls, **attributes):
        pipeline = cls()
        pipeline.logger = MagicMock()
        for key, value in attributes.items():
            setattr(pipeline, key, value)
        return pipeline

    def _pop_vcf(self, tmp_path, text=None):
        pop_vcf = tmp_path / "pop.vcf.gz"
        pop_vcf.write_text(self.HEADER if text is None else text)
        return pop_vcf

    def test_the_header_and_the_contigs_are_parsed(self, cls, tmp_path):
        pipeline = self._pipeline(cls)
        pipeline.load_pop_vcf_header(self._pop_vcf(tmp_path))
        assert pipeline.pop_vcf_contigs == {"chr1": 1000}
        assert pipeline.pop_vcf_header is not None

    def test_no_pop_vcf_leaves_both_empty(self, cls):
        pipeline = self._pipeline(cls)
        pipeline.load_pop_vcf_header(None)
        assert pipeline.pop_vcf_contigs == {}
        assert pipeline.pop_vcf_header is None

    def test_a_dry_run_falls_back_on_an_unparsable_pop_vcf(
        self, cls, tmp_path
    ):
        """Unit tests and dry runs supply placeholder files"""
        pipeline = self._pipeline(cls, dry_run=True)
        placeholder = tmp_path / "placeholder.vcf.gz"
        placeholder.touch()
        pipeline.load_pop_vcf_header(placeholder)
        assert pipeline.pop_vcf_header is None
        assert pipeline.pop_vcf_contigs == {}
        pipeline.logger.debug.assert_called()
        pipeline.logger.error.assert_not_called()

    def test_an_unparsable_pop_vcf_exits_outside_a_dry_run(
        self, cls, tmp_path
    ):
        pipeline = self._pipeline(cls)
        placeholder = tmp_path / "placeholder.vcf.gz"
        placeholder.touch()
        with pytest.raises(SystemExit) as excinfo:
            pipeline.load_pop_vcf_header(placeholder)
        assert excinfo.value.code == 2

    def test_a_missing_pop_vcf_exits_outside_a_dry_run(self, cls, tmp_path):
        pipeline = self._pipeline(cls)
        with pytest.raises(SystemExit) as excinfo:
            pipeline.load_pop_vcf_header(tmp_path / "absent.vcf.gz")
        assert excinfo.value.code == 2

    def test_a_matching_vcf_id_passes(self, cls, tmp_path):
        pipeline = self._pipeline(cls)
        pipeline.load_pop_vcf_header(self._pop_vcf(tmp_path))
        pipeline.check_pop_vcf_id("population-test-20260101")
        pipeline.logger.error.assert_not_called()

    def test_a_mismatching_vcf_id_exits_in_a_dry_run(self, cls, tmp_path):
        """The ID check no longer waits for a real run"""
        pipeline = self._pipeline(cls, dry_run=True)
        pipeline.load_pop_vcf_header(self._pop_vcf(tmp_path))
        with pytest.raises(SystemExit) as excinfo:
            pipeline.check_pop_vcf_id("population-other-20240101")
        assert excinfo.value.code == 2

    def test_an_unread_header_skips_the_vcf_id_check(self, cls):
        pipeline = self._pipeline(cls, dry_run=True)
        pipeline.check_pop_vcf_id("population-test-20260101")
        pipeline.logger.error.assert_not_called()
        pipeline.logger.debug.assert_called()

    def test_an_unparsable_header_does_not_raise_in_a_dry_run(
        self, cls, tmp_path
    ):
        """The whole fall-back path, as a `--dry_run` walks it"""
        pipeline = self._pipeline(cls, dry_run=True)
        with patch(
            "sentieon_cli.base_pangenome.read_vcf_header",
            side_effect=VcfHeaderError("no VCF header"),
        ):
            pipeline.load_pop_vcf_header(tmp_path / "pop.vcf.gz")
        pipeline.check_pop_vcf_id("population-test-20260101")
        pipeline.validate_pangenome_contig_lengths()
        pipeline.logger.error.assert_not_called()


if __name__ == "__main__":
    pytest.main([__file__])
