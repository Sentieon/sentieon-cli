"""
Unit tests for the BAM/CRAM re-alignment command builders
"""

import os
import pathlib
import sys

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.command_strings import (
    cmd_samtools_fastq_bwa,
    cmd_samtools_fastq_minimap2,
)


REFERENCE = pathlib.Path("/ref/hg38.fa")
INPUT_REF = pathlib.Path("/ref/hs37d5.fa")
MODEL_BUNDLE = pathlib.Path("/bundle/sample.bundle")
INPUT_ALN = pathlib.Path("/data/sample.cram")
OUT_ALN = pathlib.Path("/out/sample_mm2_sorted_0.cram")
RG_HEADER = pathlib.Path("/tmp/sample.hdr")


class TestSamtoolsFastqMinimap2:
    """Test the minimap2 re-alignment command"""

    def build(self, **kwargs) -> str:
        args = {
            "out_aln": OUT_ALN,
            "input_aln": INPUT_ALN,
            "reference": REFERENCE,
            "model_bundle": MODEL_BUNDLE,
            "cores": 4,
            "rg_lines": ["@RG\tID:rg1\tSM:sample1"],
            "sample_name": "sample1",
        }
        args.update(kwargs)
        return str(cmd_samtools_fastq_minimap2(**args))

    def test_input_ref_decodes_the_input(self):
        """`input_ref` decodes the input file, not the target reference"""
        cmd_str = self.build(input_ref=INPUT_REF)
        assert f"samtools fastq --reference {INPUT_REF}" in cmd_str
        assert f"samtools fastq --reference {REFERENCE}" not in cmd_str
        # The target reference is still used for alignment and sorting
        assert f"--reference {REFERENCE}" in cmd_str

    def test_no_input_ref(self):
        """No `--reference` is passed to samtools without `input_ref`"""
        cmd_str = self.build()
        assert "--reference" not in cmd_str.split("sentieon minimap2")[0]
        assert "samtools fastq -@ 4" in cmd_str

    def test_default_minimap2_model(self):
        """The minimap2 model defaults to the bundle's minimap2.model"""
        cmd_str = self.build()
        assert f"-x {MODEL_BUNDLE}/minimap2.model" in cmd_str

    def test_explicit_minimap2_model(self):
        """An explicit minimap2 model overrides the default"""
        model = MODEL_BUNDLE.joinpath("minimap2_lr.model")
        cmd_str = self.build(minimap2_model=model)
        assert f"-x {model}" in cmd_str
        assert f"-x {MODEL_BUNDLE}/minimap2.model" not in cmd_str

    def test_sm_backfill(self):
        """A missing SM tag is backfilled with the sample name"""
        cmd_str = self.build(rg_lines=["@RG\tID:rg1"], sample_name="sample1")
        assert "addreplacerg" in cmd_str
        assert "SM:sample1" in cmd_str

    def test_sm_not_overwritten(self):
        """An existing SM tag is left alone"""
        cmd_str = self.build(
            rg_lines=["@RG\tID:rg1\tSM:from_header"],
            sample_name="sample1",
        )
        assert "SM:from_header" in cmd_str
        assert "SM:sample1" not in cmd_str

    def test_one_addreplacerg_per_readgroup(self):
        """Each readgroup gets its own addreplacerg command"""
        cmd_str = self.build(
            rg_lines=["@RG\tID:rg1\tSM:s1", "@RG\tID:rg2\tSM:s1"],
        )
        assert cmd_str.count("samtools addreplacerg") == 2

    def test_output_and_sort(self):
        """The realigned reads are sorted into the output file"""
        cmd_str = self.build()
        assert f"-o {OUT_ALN}" in cmd_str
        assert "sentieon util sort" in cmd_str
        assert "--sam2bam" in cmd_str


class TestSamtoolsFastqBwa:
    """Test the bwa re-alignment command"""

    def build(self, **kwargs) -> str:
        args = {
            "out_aln": OUT_ALN,
            "input_aln": INPUT_ALN,
            "reference": REFERENCE,
            "model_bundle": MODEL_BUNDLE,
            "cores": 4,
            "rg_header": RG_HEADER,
        }
        args.update(kwargs)
        return str(cmd_samtools_fastq_bwa(**args))

    def test_input_ref_without_collate(self):
        """`input_ref` decodes the input file read by samtools fastq"""
        cmd_str = self.build(input_ref=INPUT_REF)
        assert f"samtools fastq --reference {INPUT_REF}" in cmd_str
        assert "samtools collate" not in cmd_str

    def test_input_ref_with_collate(self):
        """`input_ref` decodes the input file read by samtools collate"""
        cmd_str = self.build(input_ref=INPUT_REF, collate=True)
        assert f"samtools collate --reference {INPUT_REF}" in cmd_str
        # The fastq stage reads the collated stream, not the input file
        assert f"samtools fastq --reference {INPUT_REF}" not in cmd_str

    def test_no_input_ref(self):
        """No `--reference` is passed to samtools without `input_ref`"""
        cmd_str = self.build()
        assert "--reference" not in cmd_str.split("sentieon bwa")[0]
        assert "samtools fastq -@ 4" in cmd_str

    def test_bwa_model_and_header(self):
        """bwa uses the bundle model and the readgroup header file"""
        cmd_str = self.build()
        assert f"-x {MODEL_BUNDLE}/bwa.model" in cmd_str
        assert f"-H {RG_HEADER}" in cmd_str
