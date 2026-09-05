"""
Unit tests for sex-aware CNV calling in the DNAscope hybrid pipeline
"""

import os
import pathlib
import sys
import tempfile
from unittest.mock import MagicMock, patch

import pytest

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline
from sentieon_cli.util import SampleSex


class TestDNAscopeHybridCnv:
    """CNV calling runs in a second DAG, after the sample sex is known"""

    def setup_method(self):
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_sr_aln = [self.mock_dir / "short.bam"]
        self.mock_lr_aln = [self.mock_dir / "long.bam"]
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_bed = self.mock_dir / "interval.bed"
        self.mock_par_bed = self.mock_dir / "par.bed"

        for file_path in [
            self.mock_ref,
            self.mock_sr_aln[0],
            self.mock_lr_aln[0],
            self.mock_bundle,
            self.mock_bed,
            self.mock_par_bed,
        ]:
            file_path.touch()

    def create_pipeline(self):
        """Create a DNAscopeHybridPipeline for testing"""
        with patch("sys.exit"):
            pipeline = DNAscopeHybridPipeline()

        pipeline.logger = MagicMock()

        pipeline.output_vcf = self.mock_vcf
        pipeline.reference = self.mock_ref
        pipeline.model_bundle = self.mock_bundle
        pipeline.bed = self.mock_bed
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.tmp_dir = self.mock_dir
        pipeline.pop_vcf = None

        # State normally set by validate()
        pipeline.fai_data = {"chr1": {"length": 1000}}
        pipeline.shards = [MagicMock()]
        pipeline.shards[0].contig = "chr1"
        pipeline.shards[0].start = 1
        pipeline.shards[0].stop = 1000
        pipeline.pop_vcf_contigs = {"chr1": 1000}
        pipeline.lr_aln_readgroups = [[{"ID": "lr_rg1", "SM": "sample1"}]]
        pipeline.sr_aln_readgroups = [[{"ID": "sr_rg1", "SM": "sample1"}]]
        pipeline.hybrid_rg_sm = "sample1"
        pipeline.hybrid_set_rg = False
        pipeline.shortread_tech = "Illumina"
        pipeline.longread_tech = "ONT"
        pipeline.sr_aln = self.mock_sr_aln
        pipeline.lr_aln = self.mock_lr_aln
        pipeline.skip_cnv = False

        return pipeline

    def build_dag(self, pipeline):
        """Build the first DAG and index the jobs by name"""
        dag = pipeline.build_dag()
        jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        return dag, {job.name: job for job in jobs}

    def build_second_dag(self, pipeline):
        """Build the second DAG and index the jobs by name"""
        dag = pipeline.build_second_dag()
        if dag is None:
            return None, {}
        jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        return dag, {job.name: job for job in jobs}

    def test_cnv_jobs_run_in_the_second_dag(self):
        pipeline = self.create_pipeline()
        pipeline.sample_sex = SampleSex.FEMALE
        _dag, jobs = self.build_dag(pipeline)

        assert "cnvscope" not in jobs
        assert "cnv-model-apply" not in jobs

        second_dag, second_jobs = self.build_second_dag(pipeline)
        assert "cnvscope" in second_jobs
        assert "cnv-model-apply" in second_jobs
        assert second_dag.waiting_jobs[second_jobs["cnv-model-apply"]] == {
            second_jobs["cnvscope"]
        }

    def test_the_first_dag_estimates_the_sample_ploidy(self):
        pipeline = self.create_pipeline()
        _dag, jobs = self.build_dag(pipeline)

        assert "estimate-ploidy" in jobs
        cmd_str = str(jobs["estimate-ploidy"].shell)
        assert "estimate_ploidy.py" in cmd_str
        # Ploidy is estimated from the deduped short-read alignment
        assert "output_deduped.cram" in cmd_str
        assert str(pipeline.ploidy_json).endswith("output_ploidy.json")
        # The chr-prefixed script defaults apply to this reference
        assert "--contigs" not in cmd_str
        assert "--x_contig" not in cmd_str

    def test_the_ploidy_job_runs_with_sample_sex(self):
        # `--sample_sex` overrides the estimate, which still runs
        pipeline = self.create_pipeline()
        pipeline.sample_sex = SampleSex.FEMALE
        _dag, jobs = self.build_dag(pipeline)

        assert "estimate-ploidy" in jobs

    def test_ploidy_contigs_follow_the_reference_build(self):
        pipeline = self.create_pipeline()
        pipeline.reference_build = "b37"
        _dag, jobs = self.build_dag(pipeline)

        cmd_str = str(jobs["estimate-ploidy"].shell)
        assert "--contigs 1 2 " in cmd_str
        assert "--autosomes 1 2 " in cmd_str
        assert "--x_contig X" in cmd_str
        assert "--y_contig Y" in cmd_str

    def test_the_ploidy_job_runs_when_cnvs_are_skipped(self):
        # The ploidy JSON is written by every run, but there is no
        # second DAG without CNV calling
        pipeline = self.create_pipeline()
        pipeline.skip_cnv = True
        _dag, jobs = self.build_dag(pipeline)

        assert "estimate-ploidy" in jobs
        assert self.build_second_dag(pipeline)[0] is None

    def test_a_female_sample_is_called_without_a_par_bed(self):
        pipeline = self.create_pipeline()
        pipeline.sample_sex = SampleSex.FEMALE
        self.build_dag(pipeline)
        _dag, jobs = self.build_second_dag(pipeline)

        cmd_str = str(jobs["cnvscope"].shell)
        assert "--sex F" in cmd_str
        assert "--par" not in cmd_str

    def test_a_male_sample_is_called_with_the_par_bed(self):
        pipeline = self.create_pipeline()
        pipeline.sample_sex = SampleSex.MALE
        pipeline.cnv_par_bed = self.mock_par_bed
        self.build_dag(pipeline)
        _dag, jobs = self.build_second_dag(pipeline)

        cmd_str = str(jobs["cnvscope"].shell)
        assert "--sex M" in cmd_str
        assert f"--par {self.mock_par_bed}" in cmd_str

    def test_a_dry_run_assumes_a_male_sample(self):
        # Validation has already confirmed a PAR BED file is available
        pipeline = self.create_pipeline()
        pipeline.cnv_par_bed = self.mock_par_bed
        self.build_dag(pipeline)
        _dag, jobs = self.build_second_dag(pipeline)

        assert pipeline.sample_sex is SampleSex.MALE
        cmd_str = str(jobs["cnvscope"].shell)
        assert "--sex M" in cmd_str
        assert f"--par {self.mock_par_bed}" in cmd_str

    def _validate_cnv(self, pipeline):
        """Run the CNV validation with the reference index parsed"""
        pipeline.fai_data = {"chr1": {"length": 1000}}
        pipeline.validate_cnv()

    @pytest.mark.parametrize(
        "sample_sex", [None, SampleSex.MALE, SampleSex.FEMALE]
    )
    def test_validation_requires_a_par_bed_for_cnv_calling(self, sample_sex):
        # The reference build is not recognized, so no packaged PAR BED
        # file can be selected and the run stops during validation
        pipeline = self.create_pipeline()
        pipeline.sample_sex = sample_sex

        with pytest.raises(SystemExit) as excinfo:
            self._validate_cnv(pipeline)
        assert excinfo.value.code == 2

    def test_validation_accepts_the_par_bed_argument(self):
        pipeline = self.create_pipeline()
        pipeline.par_bed = self.mock_par_bed

        self._validate_cnv(pipeline)  # no SystemExit

        assert pipeline.cnv_par_bed == self.mock_par_bed

    def test_validation_needs_no_par_bed_when_cnvs_are_skipped(self):
        pipeline = self.create_pipeline()
        pipeline.skip_cnv = True

        self._validate_cnv(pipeline)  # no SystemExit

        assert pipeline.cnv_par_bed is None
