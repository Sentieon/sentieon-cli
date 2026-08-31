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
from sentieon_cli.base_pangenome import SampleSex
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
        pipeline.has_cnv_model = True
        pipeline.pcr_free = False

        return pipeline

    def create_fastq_pipeline(self):
        """Create a fastq-input SentieonPangenome pipeline for testing.

        The default fixture uses BAM/CRAM input, which does not exercise the
        dedup/metrics branch of ``build_first_dag``. This variant provides
        FASTQ input so the bwa/mm2 dedup jobs are created.
        """
        mock_r1 = self.mock_dir / "sample_R1.fastq.gz"
        mock_r2 = self.mock_dir / "sample_R2.fastq.gz"
        for fq in (mock_r1, mock_r2):
            fq.touch()

        pipeline = self.create_pipeline()
        pipeline.sample_input = []
        pipeline.r1_fastq = [mock_r1]
        pipeline.r2_fastq = [mock_r2]
        pipeline.fastq_readgroup = {"ID": "rg1", "SM": "sample1"}
        pipeline.skip_metrics = False
        pipeline.skip_multiqc = True
        return pipeline

    def test_dedup_metrics_output(self):
        """Dedup on the primary (bwa) alignment emits a --metrics file that
        lands in the metrics directory scanned by MultiQC."""
        pipeline = self.create_fastq_pipeline()
        dag = pipeline.build_first_dag()

        job_names, all_jobs = self._get_all_job_names(dag)
        assert "dedup-bwa" in job_names
        assert "dedup-mm2" in job_names

        bwa_dedup = next(j for j in all_jobs if j.name == "dedup-bwa")
        bwa_cmd = str(bwa_dedup.shell)
        assert "--algo Dedup" in bwa_cmd
        assert "--metrics" in bwa_cmd
        assert "output_metrics/output.txt.dedup_metrics.txt" in bwa_cmd

        # The mm2 dedup does not emit a duplicate-metrics file
        mm2_dedup = next(j for j in all_jobs if j.name == "dedup-mm2")
        assert "--metrics" not in str(mm2_dedup.shell)

    def test_metrics_input_bwa_only(self):
        """The metrics job reads only the bwa alignment; including the mm2
        alignment would count the extracted reads it re-aligns twice."""
        pipeline = self.create_fastq_pipeline()
        dag = pipeline.build_first_dag()

        _, all_jobs = self._get_all_job_names(dag)
        metrics_job = next(j for j in all_jobs if j.name == "metrics")
        metrics_cmd = str(metrics_job.shell)
        assert "output_bwa_deduped.cram" in metrics_cmd
        assert "output_mm2_deduped.cram" not in metrics_cmd

    def test_dedup_metrics_skipped(self):
        """No Dedup --metrics output when metrics collection is skipped."""
        pipeline = self.create_fastq_pipeline()
        pipeline.skip_metrics = True
        dag = pipeline.build_first_dag()

        _, all_jobs = self._get_all_job_names(dag)
        bwa_dedup = next(j for j in all_jobs if j.name == "dedup-bwa")
        assert "--metrics" not in str(bwa_dedup.shell)

    def test_model_apply_default(self):
        """Model apply runs by default and writes the
        intermediate, not the final output VCF."""
        pipeline = self.create_pipeline()
        dag = pipeline.build_first_dag()

        # Verify DAG is created
        assert isinstance(dag, DAG)

        # Check for model-apply job
        job_names = [job.name for job in dag.waiting_jobs.keys()]
        assert "model-apply" in job_names

        apply_job = next(
            j for j in dag.waiting_jobs if j.name == "model-apply"
        )
        cmd_str = str(apply_job.shell)
        assert str(self.snv_apply_vcf(pipeline)) in cmd_str
        # The final output is written by the AD-update job
        assert str(pipeline.output_vcf) not in cmd_str

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

        # Verify that the concat job writes the intermediate,
        # which the AD-update job then rewrites to the final output
        # Find the merge-trim-concat job
        concat_job = None
        for job in dag.waiting_jobs:
            if job.name == "merge-trim-concat":
                concat_job = job
                break

        assert concat_job is not None
        # Check that the first argument (output file) is set to the intermediate
        args_str = " ".join([str(arg) for arg in concat_job.shell.nodes[0].args])
        assert str(self.snv_apply_vcf(pipeline)) in args_str
        assert str(pipeline.output_vcf) not in args_str

        # The AD-update job writes the final output and depends
        # on the concat job
        _, all_jobs = self._get_all_job_names(dag)
        ad_update_job = next(
            j for j in all_jobs if j.name == "sad-lad-update"
        )
        update_cmd = str(ad_update_job.shell)
        assert f"--input_vcf {self.snv_apply_vcf(pipeline)}" in update_cmd
        assert f"--output_vcf {pipeline.output_vcf}" in update_cmd
        assert dag.waiting_jobs[ad_update_job] == {concat_job}

    def test_gvcf_adds_gvcftyper_and_routes_outputs(self):
        """--gvcf produces a .g.vcf.gz from model-apply,
        updates it into the final .g.vcf.gz, then runs GVCFtyper on the
        updated gVCF to produce the final VCF at output_vcf."""
        pipeline = self.create_pipeline()
        pipeline.gvcf = True
        dag = pipeline.build_first_dag()

        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        job_names = [job.name for job in all_jobs]
        assert "model-apply" in job_names
        assert "sad-lad-update" in job_names
        assert "gvcftyper" in job_names

        expected_gvcf = str(pipeline.output_vcf).replace(
            ".vcf.gz", ".g.vcf.gz"
        )
        apply_gvcf = str(self.snv_apply_vcf(pipeline))
        assert apply_gvcf.endswith("sample-snv_apply.g.vcf.gz")

        # DNAscope emits gVCF; model-apply writes gVCF
        for name in ("dnascope-raw", "model-apply"):
            job = next(j for j in all_jobs if j.name == name)
            cmd_str = str(job.shell)
            if name == "dnascope-raw":
                assert "--emit_mode gvcf" in cmd_str
            else:
                assert apply_gvcf in cmd_str
                assert expected_gvcf not in cmd_str

        # The update job rewrites it to the final .g.vcf.gz
        apply_job = next(j for j in all_jobs if j.name == "model-apply")
        ad_update_job = next(
            j for j in all_jobs if j.name == "sad-lad-update"
        )
        update_cmd = str(ad_update_job.shell)
        assert f"--input_vcf {apply_gvcf}" in update_cmd
        assert f"--output_vcf {expected_gvcf}" in update_cmd
        assert dag.waiting_jobs[ad_update_job] == {apply_job}

        # GVCFtyper reads the updated gVCF and writes the final VCF
        gvcftyper_job = next(j for j in all_jobs if j.name == "gvcftyper")
        cmd_str = str(gvcftyper_job.shell)
        assert "--algo GVCFtyper" in cmd_str
        assert expected_gvcf in cmd_str
        assert apply_gvcf not in cmd_str
        assert str(pipeline.output_vcf) in cmd_str
        assert dag.waiting_jobs[gvcftyper_job] == {ad_update_job}

    def snv_apply_vcf(self, pipeline):
        """The small-variant model-apply intermediate for a pipeline"""
        name = (
            "sample-snv_apply.g.vcf.gz"
            if getattr(pipeline, "gvcf", False)
            else "sample-snv_apply.vcf.gz"
        )
        return pipeline.tmp_dir.joinpath(name)

    def test_ad_update_default(self):
        """A single AD-update job reads the model-apply
        intermediate and writes the final output VCF."""
        pipeline = self.create_pipeline()
        dag = pipeline.build_first_dag()

        job_names, all_jobs = self._get_all_job_names(dag)
        assert job_names.count("sad-lad-update") == 1

        ad_update_job = next(
            j for j in all_jobs if j.name == "sad-lad-update"
        )
        cmd_str = str(ad_update_job.shell)
        assert "sad_lad_update.py" in cmd_str
        assert f"--input_vcf {self.snv_apply_vcf(pipeline)}" in cmd_str
        assert f"--output_vcf {pipeline.output_vcf}" in cmd_str
        assert f"--threads {pipeline.cores}" in cmd_str
        assert ad_update_job.task_name == "ad-update"

        # Model apply writes the intermediate, not the final output
        apply_job = next(j for j in all_jobs if j.name == "model-apply")
        assert str(self.snv_apply_vcf(pipeline)) in str(apply_job.shell)

        # The AD-update job depends on model-apply
        assert dag.waiting_jobs[ad_update_job] == {apply_job}

    def test_no_ad_update_with_skip_small_variants(self):
        """No AD-update job when small variants are skipped"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        dag = pipeline.build_first_dag()

        job_names, _ = self._get_all_job_names(dag)
        assert "sad-lad-update" not in job_names

    def test_no_gvcftyper_without_gvcf(self):
        """Without --gvcf, no GVCFtyper job is added"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_first_dag()

        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        job_names = [job.name for job in all_jobs]
        assert "gvcftyper" not in job_names

    def build_segdup_cmd(self, tech):
        """Build the segdup-caller command string for a given platform"""
        pipeline = self.create_pipeline()
        pipeline.tech = tech
        pipeline.sample_sex = SampleSex.FEMALE
        job = pipeline.build_segdup_job(
            self.mock_dir / "output_segdups",
            self.mock_bam,
            self.mock_vcf,
            None,
        )
        return str(job.shell)

    def test_segdup_ultima_lowers_min_map_qual(self):
        """Ultima input overrides segdup-caller's default min_map_qual"""
        assert "--set main.min_map_qual=30" in self.build_segdup_cmd("Ultima")

    def test_segdup_no_override_for_short_reads(self):
        """Non-Ultima input leaves the segdup-caller defaults alone"""
        assert "--set" not in self.build_segdup_cmd("Illumina")

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

    def test_call_svs_without_cnv_model(self):
        """SVs run but CNV jobs are skipped when bundle lacks cnv.model"""
        pipeline = self.create_pipeline()
        pipeline.call_svs = True
        pipeline.has_cnv_model = False
        dag = pipeline.build_first_dag()

        job_names, all_jobs = self._get_all_job_names(dag)

        # PangenomeSV still runs
        assert "dnascope-raw" in job_names
        dnascope_job = next(j for j in all_jobs if j.name == "dnascope-raw")
        cmd_str = str(dnascope_job.shell)
        assert "--algo PangenomeSV" in cmd_str

        # CNV jobs are skipped
        assert "cnvscope" not in job_names
        assert "cnv-model-apply" not in job_names
        assert "indel2cnv" not in job_names
        assert "combine-cnv" not in job_names

    def test_call_svs_without_cnv_model_skip_small_variants(self):
        """SV-only mode: SVs run but CNV jobs are skipped without cnv.model"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        pipeline.call_svs = True
        pipeline.has_cnv_model = False
        dag = pipeline.build_first_dag()

        job_names, _ = self._get_all_job_names(dag)
        assert "dnascope-raw" in job_names
        assert "cnvscope" not in job_names
        assert "cnv-model-apply" not in job_names
        assert "indel2cnv" not in job_names
        assert "combine-cnv" not in job_names

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
        assert "cnv.model" in cmd_str
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
        assert "cnv.model" in cmd_str

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

