"""
Unit tests for the HybridPangenome pipeline logic
"""

from importlib.resources import files
import logging
import os
import pathlib
import sys
import tempfile
from typing import List
from unittest.mock import MagicMock

import packaging.version
import pytest

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli import command_strings as cmds
from sentieon_cli import hybrid_pangenome
from sentieon_cli.hybrid_pangenome import HybridPangenome
from sentieon_cli.command_strings import LONGREAD_SV_BED_AWK
from sentieon_cli.dag import DAG
from sentieon_cli.util import SampleSex

PACKAGE_LOGGER = "sentieon_cli"


class _RecordingHandler(logging.Handler):
    """Collects formatted messages from the package logger."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.messages: List[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


@pytest.fixture
def messages():
    """`caplog` does not see the package logger, which does not propagate"""
    logger = logging.getLogger(PACKAGE_LOGGER)
    handler = _RecordingHandler()
    logger.addHandler(handler)
    try:
        yield handler.messages
    finally:
        logger.removeHandler(handler)


class TestHybridPangenome:
    """Test the hybrid-pangenome pipeline logic"""

    def setup_method(self):
        """Setup test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.mock_dir = pathlib.Path(self.temp_dir)

        self.mock_vcf = self.mock_dir / "output.vcf.gz"
        self.mock_ref = self.mock_dir / "reference.fa"
        self.mock_lr_bam = self.mock_dir / "longreads.bam"
        self.mock_bundle = self.mock_dir / "model.bundle"
        self.mock_gbz = self.mock_dir / "pangenome.grch38.gbz"
        self.mock_hapl = self.mock_dir / "pangenome.grch38.hapl"
        self.mock_pop_vcf = self.mock_dir / "population.vcf.gz"
        self.mock_dbsnp = self.mock_dir / "dbsnp.vcf.gz"
        self.mock_bed = self.mock_dir / "autosomes.bed"
        self.mock_r1 = self.mock_dir / "sample_R1.fastq.gz"
        self.mock_r2 = self.mock_dir / "sample_R2.fastq.gz"
        self.mock_sr_bam = self.mock_dir / "shortreads.bam"
        self.mock_lr_ref = self.mock_dir / "lr_reference.fa"

        for file_path in [
            self.mock_ref,
            self.mock_lr_bam,
            self.mock_bundle,
            self.mock_gbz,
            self.mock_hapl,
            self.mock_pop_vcf,
            self.mock_dbsnp,
            self.mock_bed,
            self.mock_r1,
            self.mock_r2,
            self.mock_sr_bam,
            self.mock_lr_ref,
        ]:
            file_path.touch()

        with open(str(self.mock_ref) + ".fai", "w") as f:
            f.write("chr1\t1000\t0\t80\t81\n")

    def create_pipeline(self):
        """Create a HybridPangenome pipeline for testing"""
        pipeline = HybridPangenome()

        pipeline.logger = MagicMock()

        # Configure arguments
        pipeline.output_vcf = self.mock_vcf
        pipeline.reference = self.mock_ref
        pipeline.model_bundle = self.mock_bundle
        pipeline.r1_fastq = [self.mock_r1]
        pipeline.r2_fastq = [self.mock_r2]
        pipeline.lr_aln = [self.mock_lr_bam]
        pipeline.gbz = self.mock_gbz
        pipeline.hapl = self.mock_hapl
        pipeline.pop_vcf = self.mock_pop_vcf
        pipeline.dbsnp = self.mock_dbsnp
        pipeline.bed = self.mock_bed
        pipeline.cores = 2
        pipeline.dry_run = True
        pipeline.skip_version_check = True
        pipeline.skip_multiqc = True
        pipeline.tmp_dir = self.mock_dir

        # State normally set by validate()
        pipeline.fai_data = {"chr1": {"length": 1000}}
        pipeline.shards = [MagicMock()]
        pipeline.shards[0].contig = "chr1"
        pipeline.shards[0].start = 1
        pipeline.shards[0].stop = 1000
        pipeline.pop_vcf_contigs = {"chr1": 1000}
        pipeline.fastq_readgroup = {"ID": "rg1", "SM": "sample1"}
        pipeline.lr_readgroups = [[{"ID": "lr-rg1", "SM": "sample1"}]]
        pipeline.sample_sm = "sample1"

        return pipeline

    def create_aligned_pipeline(self):
        """Create a pipeline with aligned short-read input"""
        pipeline = self.create_pipeline()
        pipeline.r1_fastq = []
        pipeline.r2_fastq = []
        pipeline.readgroup = None
        pipeline.fastq_readgroup = None
        pipeline.sr_aln = [self.mock_sr_bam]
        pipeline.sr_readgroups = [[{"ID": "sr-rg1", "SM": "sample1"}]]
        return pipeline

    def create_lr_realign_pipeline(self, aligned_sr=False):
        """Create a pipeline that realigns the long-read input"""
        pipeline = (
            self.create_aligned_pipeline()
            if aligned_sr
            else self.create_pipeline()
        )
        pipeline.lr_align_input = True
        pipeline.lr_input_ref = self.mock_lr_ref
        return pipeline

    def _get_all_job_names(self, dag):
        """Helper to get all job names from a DAG"""
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        return [job.name for job in all_jobs], all_jobs

    def _get_job(self, all_jobs, name):
        return next(j for j in all_jobs if j.name == name)

    def _get_dep_names(self, dag, all_jobs, name):
        """Helper to get the dependency names of a job"""
        job = self._get_job(all_jobs, name)
        return {dep.name for dep in dag.waiting_jobs.get(job, set())}

    def test_dag_jobs_present(self):
        """All expected jobs are in the DAG"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        assert isinstance(dag, DAG)

        job_names, _ = self._get_all_job_names(dag)
        for name in (
            "kmc",
            "bwa-extract",
            "vg-haplotypes",
            "vg-convert-gfa",
            "graph-update-raw",
            "longreadsv",
            "longread-sv-bed",
            "graph-update",
            "gfa2fa",
            "faidx",
            "mm2-lift",
            "locuscollector-bwa",
            "dedup-bwa",
            "locuscollector-lift",
            "dedup-lift",
            "metrics",
            "estimate-ploidy",
            "pangenome-sv",
            "dnascope",
            "model-apply",
        ):
            assert name in job_names, f"missing job: {name}"

    def test_kmc_command(self):
        """KMC counts the short-read fastq and long-read fasta together"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "kmc").shell)
        assert "-k29" in cmd_str
        assert "-m30" in cmd_str
        assert "-okff" in cmd_str
        assert "-fa /dev/stdin" in cmd_str
        assert "samtools fasta" in cmd_str
        assert str(self.mock_lr_bam) in cmd_str
        assert str(self.mock_r1) in cmd_str
        assert str(self.mock_r2) in cmd_str

    def test_bwa_readgroup_lr0(self):
        """The bwa alignment carries the LR:0 readgroup attribute"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "bwa-extract").shell)
        # The user's RGID is used unchanged
        assert r"ID:rg1\t" in cmd_str
        assert "ID:rg1-bwa" not in cmd_str
        assert "LR:0" in cmd_str
        assert "pgutil extract" in cmd_str
        assert "bwa.model" in cmd_str
        # The alignment is written to the temporary directory for dedup
        bwa_bam = self.mock_dir / "sample-bwa.bam"
        assert f"-b {bwa_bam}" in cmd_str

    def test_mm2_lift_readgroup_lr2(self):
        """The lifted alignment carries the LR:2 readgroup attribute"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "mm2-lift").shell)
        assert "ID:rg1-pg" in cmd_str
        assert "LR:2" in cmd_str
        assert "--secondary=yes" in cmd_str
        assert "pgutil lift" in cmd_str
        assert "--prefix" in cmd_str
        assert "GRCh38#0#" in cmd_str
        assert "minimap2.model" in cmd_str
        # `util sort` writes the sorted alignment to the temporary directory
        assert "util sort" in cmd_str
        lift_bam = self.mock_dir / "sample-lift.bam"
        assert f"-o {lift_bam}" in cmd_str

    def test_dedup_commands(self):
        """The short-read alignments are deduplicated, with Dedup metrics
        from the bwa alignment only"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        bwa_bam = self.mock_dir / "sample-bwa.bam"
        lift_bam = self.mock_dir / "sample-lift.bam"
        out_bwa = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        out_lift = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
        dedup_metrics = str(
            self.mock_dir / "output_metrics" / "output.txt.dedup_metrics.txt"
        )

        lc_cmd = str(self._get_job(all_jobs, "locuscollector-bwa").shell)
        assert "--algo LocusCollector" in lc_cmd
        assert str(bwa_bam) in lc_cmd

        bwa_cmd = str(self._get_job(all_jobs, "dedup-bwa").shell)
        assert "--algo Dedup" in bwa_cmd
        assert str(bwa_bam) in bwa_cmd
        assert out_bwa in bwa_cmd
        assert f"--metrics {dedup_metrics}" in bwa_cmd
        assert "IndelLeftAlignReadTransform" not in bwa_cmd

        lift_cmd = str(self._get_job(all_jobs, "dedup-lift").shell)
        assert "--algo Dedup" in lift_cmd
        assert str(lift_bam) in lift_cmd
        assert out_lift in lift_cmd
        assert "IndelLeftAlignReadTransform" not in lift_cmd
        # Dedup metrics are only collected from the bwa alignment
        assert "--metrics" not in lift_cmd

        # The long-read input is assumed to be deduplicated already
        for job_name in ("locuscollector-bwa", "dedup-bwa", "dedup-lift"):
            assert str(self.mock_lr_bam) not in str(
                self._get_job(all_jobs, job_name).shell
            )

    def test_metrics_job(self):
        """Metrics are collected from the deduplicated alignments"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "metrics").shell)
        assert (
            str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
            in cmd_str
        )
        assert (
            str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
            in cmd_str
        )
        assert "--algo GCBias" in cmd_str
        assert "--algo WgsMetricsAlgo" in cmd_str
        assert "rehead-metrics" in job_names

    def test_skip_metrics(self):
        """skip_metrics removes the metrics jobs and the Dedup metrics"""
        pipeline = self.create_pipeline()
        pipeline.skip_metrics = True
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        assert "metrics" not in job_names
        assert "rehead-metrics" not in job_names
        assert "multiqc" not in job_names
        # Deduplication still runs
        assert "dedup-bwa" in job_names
        assert "dedup-lift" in job_names
        assert "--metrics" not in str(
            self._get_job(all_jobs, "dedup-bwa").shell
        )

    def test_graph_update_commands(self):
        """PGHapUpdateAlgo runs without and then with the SV BED"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        raw_cmd = str(self._get_job(all_jobs, "graph-update-raw").shell)
        assert "--algo PGHapUpdateAlgo" in raw_cmd
        assert "--gfa_file " in raw_cmd
        assert "sample-hap.raw.gfa" in raw_cmd
        assert "sample-pangenome-raw.gfa" in raw_cmd
        assert "--target_bed" not in raw_cmd
        assert str(self.mock_lr_bam) in raw_cmd
        assert "--prefix" in raw_cmd
        assert "GRCh38#0#" in raw_cmd

        update_cmd = str(self._get_job(all_jobs, "graph-update").shell)
        assert "--algo PGHapUpdateAlgo" in update_cmd
        assert "--gfa_file " in update_cmd
        assert "sample-pangenome-raw.gfa" in update_cmd
        assert "--target_bed " in update_cmd
        assert "sample-sv.bed" in update_cmd
        assert "--prefix" in update_cmd
        assert "GRCh38#0#" in update_cmd

    def test_longreadsv_and_bed(self):
        """LongReadSV runs on the long reads and the awk BED script is
        retained verbatim"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        sv_cmd = str(self._get_job(all_jobs, "longreadsv").shell)
        assert "--algo LongReadSV" in sv_cmd
        assert "longreadsv.model" in sv_cmd
        assert "--min_sv_size 20" in sv_cmd
        assert str(self.mock_lr_bam) in sv_cmd

        bed_cmd = str(self._get_job(all_jobs, "longread-sv-bed").shell)
        assert 'n=split($10,a,":")' in bed_cmd
        # Sort in reference contig order for bedtools merge
        assert f"bedtools sort -faidx {self.mock_ref}.fai -i -" in bed_cmd
        assert "bedtools merge" in bed_cmd
        assert "sample-sv.bed" in bed_cmd
        # The awk script matches the validated implementation
        assert "chr[^_]+_[0-9]+_[0-9]+" in LONGREAD_SV_BED_AWK

    def test_gfa2fa_and_faidx(self):
        """The updated graph is converted to an indexed FASTA"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        gfa2fa_cmd = str(self._get_job(all_jobs, "gfa2fa").shell)
        assert "pgutil gfa2fa" in gfa2fa_cmd
        assert str(self.mock_ref) + ".fai" in gfa2fa_cmd
        assert "sample-pangenome.gfa" in gfa2fa_cmd
        assert "sample-pangenome.fa" in gfa2fa_cmd

        faidx_cmd = str(self._get_job(all_jobs, "faidx").shell)
        assert "samtools faidx" in faidx_cmd
        assert "sample-pangenome.fa" in faidx_cmd

    def test_gfa2fa_with_vg_paths(self):
        """`vg paths` replaces `pgutil gfa2fa` when selected"""
        pipeline = self.create_pipeline()
        pipeline.gfa2fa_with_vg = True
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        gfa = self.mock_dir / "sample-pangenome.gfa"
        fasta = self.mock_dir / "sample-pangenome.fa"
        gfa2fa_cmd = str(self._get_job(all_jobs, "gfa2fa").shell)
        assert gfa2fa_cmd == (f"vg paths -x {gfa} -F >'{fasta}'")
        assert "pgutil gfa2fa" not in gfa2fa_cmd
        assert "-Q" not in gfa2fa_cmd

        # The DAG shape is unchanged
        assert self._get_dep_names(dag, all_jobs, "faidx") == {"gfa2fa"}
        assert "faidx" in self._get_dep_names(dag, all_jobs, "mm2-lift")

    def test_resolve_gfa2fa_tool_hg38(self, monkeypatch):
        """A GRCh38 reference keeps `pgutil gfa2fa`, unprobed"""
        pipeline = self.create_pipeline()
        pipeline.skip_version_check = False
        pipeline.fai_data = {
            "chr1": {"length": 248956422},
            "chr2": {"length": 242193529},
            "chrX": {"length": 156040895},
        }
        probe = MagicMock()
        monkeypatch.setattr(hybrid_pangenome, "executable_version", probe)

        pipeline.resolve_gfa2fa_tool()
        assert pipeline.reference_build == "hg38"
        assert pipeline.gfa2fa_with_vg is False
        probe.assert_not_called()

    def test_resolve_gfa2fa_tool_old_driver(self, monkeypatch):
        """A non-GRCh38 reference on the GRCh38-only driver uses vg"""
        pipeline = self.create_pipeline()
        pipeline.skip_version_check = False
        monkeypatch.setattr(
            hybrid_pangenome,
            "executable_version",
            lambda cmd: packaging.version.Version("202503.04"),
        )

        pipeline.resolve_gfa2fa_tool()
        assert pipeline.reference_build is None
        assert pipeline.gfa2fa_with_vg is True

    def test_resolve_gfa2fa_tool_new_driver(self, monkeypatch):
        """A newer driver keeps `pgutil gfa2fa` for any reference"""
        pipeline = self.create_pipeline()
        pipeline.skip_version_check = False
        monkeypatch.setattr(
            hybrid_pangenome,
            "executable_version",
            lambda cmd: packaging.version.Version("202503.05"),
        )

        pipeline.resolve_gfa2fa_tool()
        assert pipeline.gfa2fa_with_vg is False

    def test_resolve_gfa2fa_tool_skip_version_check(self, monkeypatch):
        """`--skip_version_check` does not probe and prefers vg"""
        pipeline = self.create_pipeline()
        pipeline.skip_version_check = True
        probe = MagicMock()
        monkeypatch.setattr(hybrid_pangenome, "executable_version", probe)

        pipeline.resolve_gfa2fa_tool()
        assert pipeline.gfa2fa_with_vg is True
        probe.assert_not_called()

    def test_resolve_gfa2fa_tool_probe_failure(self, monkeypatch):
        """An unknown driver version prefers vg"""
        pipeline = self.create_pipeline()
        pipeline.skip_version_check = False
        monkeypatch.setattr(
            hybrid_pangenome,
            "executable_version",
            lambda cmd: None,
        )

        pipeline.resolve_gfa2fa_tool()
        assert pipeline.gfa2fa_with_vg is True

    def test_calling_inputs_and_replace_rg(self):
        """PangenomeSV and DNAscope use the bwa, lifted, and long reads,
        rewriting only the long-read readgroups with LR:1"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        lift_aln = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
        replace_arg = r"lr-rg1=ID:lr-rg1\tSM:sample1\tLR:1"

        for job_name in ("pangenome-sv", "dnascope"):
            cmd_str = str(self._get_job(all_jobs, job_name).shell)
            assert bwa_aln in cmd_str
            assert lift_aln in cmd_str
            assert str(self.mock_lr_bam) in cmd_str
            assert replace_arg in cmd_str
            # The long-read input is preceded by its --replace_rg argument
            assert cmd_str.index(replace_arg) < cmd_str.index(
                str(self.mock_lr_bam)
            )
            # The pipeline-generated alignments are not rewritten
            assert cmd_str.count("--replace_rg") == 1

    def test_pangenomesv_command(self):
        """PangenomeSV uses the updated graph and min_af"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "pangenome-sv").shell)
        assert "--algo PangenomeSV" in cmd_str
        assert "--gfa_file" in cmd_str
        assert "sample-pangenome.gfa" in cmd_str
        assert "--min_af 0.1" in cmd_str
        assert "--prefix" in cmd_str
        assert "GRCh38#0#" in cmd_str
        sv_vcf = str(self.mock_vcf).replace(".vcf.gz", "_sv.vcf.gz")
        assert sv_vcf in cmd_str
        # SV calling is not restricted to the small-variant BED
        assert f"--interval {self.mock_bed}" not in cmd_str

    def test_pangenome_contig_prefix(self):
        """`--pangenome_contig_prefix` reaches every graph consumer"""
        pipeline = self.create_pipeline()
        pipeline.pangenome_contig_prefix = "CHM13#0#"
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        for job_name in (
            "mm2-lift",
            "graph-update-raw",
            "graph-update",
            "pangenome-sv",
        ):
            cmd_str = str(self._get_job(all_jobs, job_name).shell)
            assert "--prefix" in cmd_str, job_name
            assert "CHM13#0#" in cmd_str, job_name
            assert "GRCh38#0#" not in cmd_str, job_name

    def test_dnascope_command(self):
        """DNAscope runs with the model, interval, and pcr_indel_model"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert "--algo DNAscope" in cmd_str
        assert "dnascope.model" in cmd_str
        assert "--pcr_indel_model CONSERVATIVE" in cmd_str
        assert f"--interval {self.mock_bed}" in cmd_str
        assert f"--dbsnp {self.mock_dbsnp}" in cmd_str

    def test_pcr_free(self):
        """--pcr_free calls DNAscope with `--pcr_indel_model NONE`"""
        pipeline = self.create_pipeline()
        pipeline.pcr_free = True
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert "--pcr_indel_model NONE" in cmd_str

    def test_bam_format(self):
        """--bam_format switches the deduplicated outputs to BAM"""
        pipeline = self.create_pipeline()
        pipeline.bam_format = True
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        bwa_bam = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.bam")
        lift_bam = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.bam")

        assert bwa_bam in str(self._get_job(all_jobs, "dedup-bwa").shell)
        assert lift_bam in str(self._get_job(all_jobs, "dedup-lift").shell)

        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert bwa_bam in cmd_str
        assert lift_bam in cmd_str

    def test_model_apply_writes_output_vcf(self):
        """DNAModelApply produces the final output VCF"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "model-apply").shell)
        assert "--algo DNAModelApply" in cmd_str
        assert str(self.mock_vcf) in cmd_str

    def test_skip_model_apply(self):
        """With skip_model_apply the transfer writes the final VCF"""
        pipeline = self.create_pipeline()
        pipeline.skip_model_apply = True
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        assert "model-apply" not in job_names
        concat_job = self._get_job(all_jobs, "merge-trim-concat")
        assert str(pipeline.output_vcf) in str(concat_job.shell)

    def test_skip_svs(self):
        """skip_svs removes PangenomeSV but not the graph jobs"""
        pipeline = self.create_pipeline()
        pipeline.skip_svs = True
        dag = pipeline.build_dag()
        job_names, _ = self._get_all_job_names(dag)

        assert "pangenome-sv" not in job_names
        assert "dnascope" in job_names
        assert "graph-update" in job_names

    def test_skip_small_variants(self):
        """skip_small_variants removes DNAscope, transfer, and model-apply"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        dag = pipeline.build_dag()
        job_names, _ = self._get_all_job_names(dag)

        assert "dnascope" not in job_names
        assert "model-apply" not in job_names
        assert "merge-trim-concat" not in job_names
        assert "pangenome-sv" in job_names

    def test_multiple_lr_inputs(self):
        """Each long-read input gets its own --replace_rg arguments"""
        lr_bam2 = self.mock_dir / "longreads2.bam"
        lr_bam2.touch()

        pipeline = self.create_pipeline()
        pipeline.lr_aln = [self.mock_lr_bam, lr_bam2]
        pipeline.lr_readgroups = [
            [{"ID": "lr-rg1", "SM": "sample1"}],
            [{"ID": "lr-rg2", "SM": "sample1"}],
        ]
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert r"lr-rg1=ID:lr-rg1\tSM:sample1\tLR:1" in cmd_str
        assert r"lr-rg2=ID:lr-rg2\tSM:sample1\tLR:1" in cmd_str
        assert str(lr_bam2) in cmd_str

    def test_rgsm_overrides_sm(self):
        """--rgsm overrides the SM tag in the rewritten readgroups"""
        pipeline = self.create_pipeline()
        pipeline.rgsm = "override_sm"
        pipeline.sample_sm = "override_sm"
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert r"lr-rg1=ID:lr-rg1\tSM:override_sm\tLR:1" in cmd_str

        bwa_cmd = str(self._get_job(all_jobs, "bwa-extract").shell)
        assert "SM:override_sm" in bwa_cmd

    # Aligned short-read input

    def test_aligned_dag_jobs(self):
        """Aligned short-read input skips alignment, dedup, and metrics"""
        pipeline = self.create_aligned_pipeline()
        assert not pipeline.skip_metrics
        dag = pipeline.build_dag()
        job_names, _ = self._get_all_job_names(dag)

        for name in (
            "extract-kmc-symlink",
            "extract-kmc",
            "vg-haplotypes",
            "graph-update",
            "mm2-lift",
            "estimate-ploidy",
            "pangenome-sv",
            "dnascope",
            "model-apply",
        ):
            assert name in job_names, f"missing job: {name}"

        for name in (
            "kmc",
            "bwa-extract",
            "locuscollector-bwa",
            "dedup-bwa",
            "locuscollector-lift",
            "dedup-lift",
            "metrics",
            "rehead-metrics",
            "multiqc",
        ):
            assert name not in job_names, f"unexpected job: {name}"

    def test_aligned_extract_kmc_command(self):
        """Read extraction and k-mer counting run in a single pass"""
        pipeline = self.create_aligned_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        symlink_cmd = str(self._get_job(all_jobs, "extract-kmc-symlink").shell)
        rw_bam = self.mock_dir / "extract-kmc-rw.bam"
        assert f"ln -sf /dev/stdout {rw_bam}" in symlink_cmd

        cmd_str = str(self._get_job(all_jobs, "extract-kmc").shell)
        # One pass over the aligned short reads, without duplicate,
        # secondary, or supplementary reads
        assert "--algo ReadWriter" in cmd_str
        assert "--output_flag_filter 0xf00:0" in cmd_str
        assert str(self.mock_sr_bam) in cmd_str
        assert str(rw_bam) in cmd_str
        # writing the extracted reads to the fastq and the reads to kmc
        assert "pgutil extract" in cmd_str
        ext_fastq = self.mock_dir / "sample-extract.fq.gz"
        assert f"-o {ext_fastq}" in cmd_str
        assert "-a -" in cmd_str
        # concatenated with the long reads
        assert "cat " in cmd_str
        assert "samtools fasta" in cmd_str
        assert str(self.mock_lr_bam) in cmd_str
        assert "-fa /dev/stdin" in cmd_str
        assert "-k29" in cmd_str
        assert "-m30" in cmd_str

    def test_aligned_extract_kmc_memory(self):
        """`--kmer_memory` is passed to the single-pass KMC"""
        pipeline = self.create_aligned_pipeline()
        pipeline.kmer_memory = 64
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "extract-kmc").shell)
        assert "-m64" in cmd_str

    def test_aligned_dag_dependencies(self):
        """The pangenome is built from the extracted k-mers"""
        pipeline = self.create_aligned_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        assert self._get_dep_names(dag, all_jobs, "extract-kmc") == {
            "extract-kmc-symlink"
        }
        assert self._get_dep_names(dag, all_jobs, "vg-haplotypes") == {
            "extract-kmc"
        }
        assert "extract-kmc" in self._get_dep_names(dag, all_jobs, "mm2-lift")
        assert self._get_dep_names(dag, all_jobs, "dnascope") == {"mm2-lift"}

    def test_aligned_lift_output(self):
        """The lifted alignment is the final short-read output"""
        pipeline = self.create_aligned_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        lift_cram = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
        cmd_str = str(self._get_job(all_jobs, "mm2-lift").shell)
        assert f"-o {lift_cram}" in cmd_str
        # The readgroup is seeded from the first input readgroup
        assert "ID:sr-rg1-pg" in cmd_str
        assert "LR:2" in cmd_str

    def test_aligned_replace_rg(self):
        """The aligned short reads are rewritten with LR:0"""
        pipeline = self.create_aligned_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        lift_cram = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
        sr_arg = r"sr-rg1=ID:sr-rg1\tSM:sample1\tLR:0"
        lr_arg = r"lr-rg1=ID:lr-rg1\tSM:sample1\tLR:1"

        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert cmd_str.count("--replace_rg") == 2
        # Each row precedes the input file it applies to
        assert cmd_str.index(sr_arg) < cmd_str.index(str(self.mock_sr_bam))
        assert cmd_str.index(lr_arg) < cmd_str.index(str(self.mock_lr_bam))
        # The calling inputs are ordered short reads, lifted, long reads
        assert (
            cmd_str.index(str(self.mock_sr_bam))
            < cmd_str.index(lift_cram)
            < cmd_str.index(str(self.mock_lr_bam))
        )
        # The lifted alignment carries its LR tag already
        assert "LR:2" not in cmd_str

    def test_aligned_bam_format(self):
        """--bam_format switches the lifted output to BAM"""
        pipeline = self.create_aligned_pipeline()
        pipeline.bam_format = True
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        lift_bam = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.bam")
        assert lift_bam in str(self._get_job(all_jobs, "mm2-lift").shell)
        assert lift_bam in str(self._get_job(all_jobs, "dnascope").shell)

    # Short-read input validation

    def test_validate_sr_inputs_fastq(self):
        """fastq input with a readgroup is accepted"""
        pipeline = self.create_pipeline()
        pipeline.readgroup = r"@RG\tID:rg1\tSM:sample1"
        pipeline.validate_sr_inputs()

    def test_validate_sr_inputs_aligned(self):
        """Aligned input without a readgroup is accepted"""
        pipeline = self.create_aligned_pipeline()
        pipeline.validate_sr_inputs()

    def test_validate_sr_inputs_both(self):
        """fastq and aligned short reads cannot be combined"""
        pipeline = self.create_pipeline()
        pipeline.readgroup = r"@RG\tID:rg1\tSM:sample1"
        pipeline.sr_aln = [self.mock_sr_bam]
        with pytest.raises(SystemExit):
            pipeline.validate_sr_inputs()

    def test_validate_sr_inputs_neither(self):
        """Short reads are required"""
        pipeline = self.create_pipeline()
        pipeline.r1_fastq = []
        pipeline.r2_fastq = []
        with pytest.raises(SystemExit):
            pipeline.validate_sr_inputs()

    def test_validate_sr_inputs_fastq_without_readgroup(self):
        """fastq input requires a readgroup"""
        pipeline = self.create_pipeline()
        pipeline.readgroup = None
        with pytest.raises(SystemExit):
            pipeline.validate_sr_inputs()

    def test_validate_sr_inputs_fastq_length_mismatch(self):
        """The r1 and r2 fastq lists must have the same length"""
        pipeline = self.create_pipeline()
        pipeline.readgroup = r"@RG\tID:rg1\tSM:sample1"
        pipeline.r2_fastq = []
        with pytest.raises(SystemExit):
            pipeline.validate_sr_inputs()

    def test_validate_sr_inputs_aligned_with_readgroup(self):
        """`--readgroup` cannot be used with aligned input"""
        pipeline = self.create_aligned_pipeline()
        pipeline.readgroup = r"@RG\tID:rg1\tSM:sample1"
        with pytest.raises(SystemExit):
            pipeline.validate_sr_inputs()

    # Unaligned (uBAM/uCRAM) long-read input

    def test_lr_realign_job(self):
        """The long-read input is realigned with minimap2"""
        pipeline = self.create_lr_realign_pipeline()
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        assert "bam-realign-0" in job_names
        cmd_str = str(self._get_job(all_jobs, "bam-realign-0").shell)
        # The input reference decodes the input file
        assert f"samtools fastq --reference {self.mock_lr_ref}" in cmd_str
        # The long-read minimap2 model aligns to the linear reference
        assert "minimap2_lr.model" in cmd_str
        assert str(self.mock_ref) in cmd_str
        realigned = str(self.mock_vcf).replace(".vcf.gz", "_mm2_sorted_0.cram")
        assert f"-o {realigned}" in cmd_str
        assert self._get_dep_names(dag, all_jobs, "bam-realign-0") == set()

    def test_lr_realign_downstream(self):
        """Downstream jobs consume the realigned long reads"""
        pipeline = self.create_lr_realign_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        realigned = str(self.mock_vcf).replace(".vcf.gz", "_mm2_sorted_0.cram")
        for job_name in (
            "longreadsv",
            "graph-update-raw",
            "graph-update",
            "pangenome-sv",
            "dnascope",
        ):
            cmd_str = str(self._get_job(all_jobs, job_name).shell)
            assert realigned in cmd_str, job_name
            assert str(self.mock_lr_bam) not in cmd_str, job_name

        for job_name in ("longreadsv", "graph-update-raw", "dnascope"):
            assert "bam-realign-0" in self._get_dep_names(
                dag, all_jobs, job_name
            ), job_name

        # The readgroups of the realigned input are unchanged
        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert r"lr-rg1=ID:lr-rg1\tSM:sample1\tLR:1" in cmd_str

    def test_lr_realign_kmc_reads_original_input(self):
        """K-mer counting reads the original long-read input"""
        pipeline = self.create_lr_realign_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "kmc").shell)
        assert f"samtools fasta --reference {self.mock_lr_ref}" in cmd_str
        assert str(self.mock_lr_bam) in cmd_str
        assert "_mm2_sorted_0" not in cmd_str
        assert self._get_dep_names(dag, all_jobs, "kmc") == set()

    def test_lr_realign_without_input_ref(self):
        """The target reference decodes the input without `lr_input_ref`"""
        pipeline = self.create_lr_realign_pipeline()
        pipeline.lr_input_ref = None
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "kmc").shell)
        assert f"samtools fasta --reference {self.mock_ref}" in cmd_str

    def test_aligned_sr_with_lr_realign(self):
        """Aligned short reads combine with realigned long reads"""
        pipeline = self.create_lr_realign_pipeline(aligned_sr=True)
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        assert "bam-realign-0" in job_names
        assert "extract-kmc" in job_names
        assert "kmc" not in job_names
        assert "dedup-bwa" not in job_names

        realigned = str(self.mock_vcf).replace(".vcf.gz", "_mm2_sorted_0.cram")
        lift_cram = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
        cmd_str = str(self._get_job(all_jobs, "dnascope").shell)
        assert (
            cmd_str.index(str(self.mock_sr_bam))
            < cmd_str.index(lift_cram)
            < cmd_str.index(realigned)
        )
        assert self._get_dep_names(dag, all_jobs, "dnascope") == {
            "mm2-lift",
            "bam-realign-0",
        }

        # k-mer counting still reads the original long-read input
        extract_cmd = str(self._get_job(all_jobs, "extract-kmc").shell)
        assert str(self.mock_lr_bam) in extract_cmd
        assert f"--reference {self.mock_lr_ref}" in extract_cmd

    # Readgroup validation
    #
    # The readgroups are read from real input headers, so these tests set
    # the parsed readgroups directly (or patch `get_rg_lines`) rather than
    # relying on the synthetic readgroup of a dry run.

    def create_rg_pipeline(self):
        """A pipeline with aligned inputs and readgroup checks enabled"""
        pipeline = self.create_aligned_pipeline()
        pipeline.dry_run = False
        pipeline.sr_readgroups = [[{"ID": "sr-rg1", "SM": "sample1"}]]
        pipeline.lr_readgroups = [[{"ID": "lr-rg1", "SM": "sample1"}]]
        return pipeline

    def test_validate_readgroups_ok(self):
        """Unique readgroups with a shared SM tag are accepted"""
        pipeline = self.create_rg_pipeline()
        pipeline.validate_readgroups()
        assert pipeline.sample_sm == "sample1"

    def test_validate_readgroups_sr_without_rg_lines(self):
        """An aligned short-read input without @RG lines is rejected"""
        pipeline = self.create_rg_pipeline()
        pipeline.sr_readgroups = [[]]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_lr_without_rg_lines(self):
        """A long-read input without @RG lines is rejected"""
        pipeline = self.create_rg_pipeline()
        pipeline.lr_readgroups = [[]]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_no_rg_lines_with_rgsm(self):
        """`--rgsm` does not bypass the missing readgroup check"""
        pipeline = self.create_rg_pipeline()
        pipeline.rgsm = "override_sm"
        pipeline.sr_readgroups = [[]]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_duplicate_ids_across_inputs(self):
        """A readgroup ID shared by the short and long reads is rejected"""
        pipeline = self.create_rg_pipeline()
        pipeline.sr_readgroups = [[{"ID": "rg1", "SM": "sample1"}]]
        pipeline.lr_readgroups = [[{"ID": "rg1", "SM": "sample1"}]]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_duplicate_ids_within_input(self):
        """A readgroup ID repeated inside one input is rejected"""
        pipeline = self.create_rg_pipeline()
        pipeline.sr_readgroups = [
            [{"ID": "rg1", "SM": "sample1"}, {"ID": "rg1", "SM": "sample1"}]
        ]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_duplicate_with_fastq_readgroup(self):
        """A long-read ID matching the `--readgroup` ID is rejected"""
        pipeline = self.create_pipeline()
        pipeline.dry_run = False
        pipeline.fastq_readgroup = {"ID": "rg1", "SM": "sample1"}
        pipeline.lr_readgroups = [[{"ID": "rg1", "SM": "sample1"}]]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_lift_id_collision(self):
        """An input ID reserved for the lifted alignment is rejected"""
        pipeline = self.create_pipeline()
        pipeline.dry_run = False
        pipeline.fastq_readgroup = {"ID": "rg1", "SM": "sample1"}
        pipeline.lr_readgroups = [[{"ID": "rg1-pg", "SM": "sample1"}]]
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_validate_readgroups_without_sample_name(self):
        """A sample name is required for the output readgroups"""
        pipeline = self.create_rg_pipeline()
        pipeline.fastq_readgroup = None
        pipeline.sr_aln = []
        pipeline.lr_aln = []
        pipeline.sr_readgroups = []
        pipeline.lr_readgroups = []
        with pytest.raises(SystemExit):
            pipeline.validate_readgroups()

    def test_collect_readgroups_malformed_rg_line(self, monkeypatch):
        """A malformed @RG line in an input header is rejected"""
        pipeline = self.create_aligned_pipeline()
        pipeline.dry_run = False
        monkeypatch.setattr(
            cmds,
            "get_rg_lines",
            lambda aln, dry_run: ["@RG\tID:sr-rg1\t"],
        )
        with pytest.raises(SystemExit):
            pipeline.collect_readgroups()

    def test_collect_readgroups_parses_input_headers(self, monkeypatch):
        """The @RG lines of every input are parsed"""
        pipeline = self.create_aligned_pipeline()
        pipeline.dry_run = False
        headers = {
            str(self.mock_sr_bam): ["@RG\tID:sr-rg1\tSM:sample1"],
            str(self.mock_lr_bam): ["@RG\tID:lr-rg1\tSM:sample1"],
        }
        monkeypatch.setattr(
            cmds,
            "get_rg_lines",
            lambda aln, dry_run: headers[str(aln)],
        )
        pipeline.collect_readgroups()
        assert pipeline.sr_readgroups == [[{"ID": "sr-rg1", "SM": "sample1"}]]
        assert pipeline.lr_readgroups == [[{"ID": "lr-rg1", "SM": "sample1"}]]
        pipeline.validate_readgroups()
        assert pipeline.sample_sm == "sample1"

    # Ploidy estimation

    def test_estimate_ploidy_job_with_fastq_input(self):
        """Ploidy is estimated from the deduplicated bwa alignment"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        job = self._get_job(all_jobs, "estimate-ploidy")
        assert "estimate_ploidy.py" in str(job.shell)
        assert f"-i {bwa_aln}" in str(job.shell)
        assert job.task_name == "ploidy"
        assert job.threads == 0
        assert self._get_dep_names(dag, all_jobs, "estimate-ploidy") == {
            "dedup-bwa"
        }

    def test_estimate_ploidy_job_with_aligned_input(self):
        """Aligned short reads are used as-is, without dependencies"""
        pipeline = self.create_aligned_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        job = self._get_job(all_jobs, "estimate-ploidy")
        assert f"-i {self.mock_sr_bam}" in str(job.shell)
        assert self._get_dep_names(dag, all_jobs, "estimate-ploidy") == set()

    def test_ploidy_json_stashed(self):
        """The ploidy JSON path is stashed for the second DAG"""
        pipeline = self.create_pipeline()
        pipeline.build_dag()

        assert pipeline.ploidy_json == pathlib.Path(
            str(self.mock_vcf).replace(".vcf.gz", "_ploidy.json")
        )

    def test_estimate_ploidy_with_skip_small_variants(self):
        """Ploidy estimation runs before the small-variant early return"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        job_names, _ = self._get_all_job_names(pipeline.build_dag())

        assert "estimate-ploidy" in job_names

    # T1K HLA/KIR genotyping

    def enable_t1k(self, pipeline):
        """Supply the T1K reference files"""
        for attr, name in (
            ("t1k_hla_seq", "hla_seq.fa"),
            ("t1k_hla_coord", "hla_coord.fa"),
            ("t1k_kir_seq", "kir_seq.fa"),
            ("t1k_kir_coord", "kir_coord.fa"),
        ):
            path = self.mock_dir / name
            path.touch()
            setattr(pipeline, attr, path)
        return pipeline

    def test_no_t1k_jobs_by_default(self):
        """T1K runs only when its reference files are supplied"""
        pipeline = self.create_pipeline()
        job_names, _ = self._get_all_job_names(pipeline.build_dag())

        assert not [name for name in job_names if name.startswith("t1k")]

    def test_t1k_jobs_with_fastq_input(self):
        """T1K genotypes the deduplicated short reads"""
        pipeline = self.enable_t1k(self.create_pipeline())
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        for name in (
            "t1k-hla-extract",
            "t1k-hla",
            "t1k-kir-extract",
            "t1k-kir",
        ):
            assert name in job_names, f"missing job: {name}"

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        extract = self._get_job(all_jobs, "t1k-hla-extract")
        assert str(extract.shell) == (
            f"sentieon driver --input {bwa_aln} "
            f"--reference {self.mock_ref} --thread_count {pipeline.cores} "
            "--interval chr6:28510020-33480577 "
            f"--algo ReadWriter {self.mock_dir}/sample_hla.bam"
        )
        assert str(self._get_job(all_jobs, "t1k-hla").shell) == (
            f"run-t1k --abnormalUnmapFlag -t {pipeline.cores} "
            "--preset hla-wgs "
            f"-f {self.mock_dir}/hla_seq.fa "
            f"-c {self.mock_dir}/hla_coord.fa "
            f"--od {self.mock_dir}/output_hla "
            f"-b {self.mock_dir}/sample_hla.bam"
        )
        assert str(self._get_job(all_jobs, "t1k-kir").shell) == (
            f"run-t1k --abnormalUnmapFlag -t {pipeline.cores} "
            "--preset kir-wgs "
            f"-f {self.mock_dir}/kir_seq.fa "
            f"-c {self.mock_dir}/kir_coord.fa "
            f"--od {self.mock_dir}/output_kir "
            f"-b {self.mock_dir}/sample_kir.bam"
        )

        assert self._get_dep_names(dag, all_jobs, "t1k-hla") == {
            "t1k-hla-extract"
        }
        assert self._get_dep_names(dag, all_jobs, "t1k-hla-extract") == {
            "dedup-bwa"
        }

    def test_t1k_jobs_with_aligned_input(self):
        """T1K extracts from the `--sr_aln` input without dependencies"""
        pipeline = self.enable_t1k(self.create_aligned_pipeline())
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        extract = self._get_job(all_jobs, "t1k-kir-extract")
        assert str(extract.shell) == (
            f"sentieon driver --input {self.mock_sr_bam} "
            f"--reference {self.mock_ref} --thread_count {pipeline.cores} "
            "--interval chr19:53100000-55800000 "
            f"--algo ReadWriter {self.mock_dir}/sample_kir.bam"
        )
        assert self._get_dep_names(dag, all_jobs, "t1k-kir-extract") == set()
        assert self._get_dep_names(dag, all_jobs, "t1k-kir") == {
            "t1k-kir-extract"
        }

    # The second-DAG stashes

    def test_stashes_with_fastq_input(self):
        """The bwa alignment is stashed without a `--replace_rg` row"""
        pipeline = self.create_pipeline()
        pipeline.build_dag()

        bwa_aln = pathlib.Path(
            str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        )
        assert pipeline.sr_alignments == [bwa_aln]
        assert pipeline.sr_replace_rg is None
        assert pipeline.lr_alignment == self.mock_lr_bam

    def test_stashes_with_aligned_input(self):
        """The `--sr_aln` input is stashed with its `--replace_rg` row"""
        pipeline = self.create_aligned_pipeline()
        pipeline.build_dag()

        assert pipeline.sr_alignments == [self.mock_sr_bam]
        assert pipeline.sr_replace_rg == [
            [pipeline._replace_rg_arg({"ID": "sr-rg1", "SM": "sample1"}, "0")]
        ]
        assert pipeline.sr_replace_rg == [
            [r"sr-rg1=ID:sr-rg1\tSM:sample1\tLR:0"]
        ]

    def test_lr_alignment_stash_with_realignment(self):
        """The realigned long-read alignment is stashed"""
        pipeline = self.create_lr_realign_pipeline()
        pipeline.build_dag()

        assert pipeline.lr_alignment == pathlib.Path(
            str(self.mock_vcf).replace(".vcf.gz", "_mm2_sorted_0.cram")
        )

    def test_lr_alignment_stash_with_two_inputs(self):
        """segdup-caller takes a single long-read alignment"""
        lr_bam2 = self.mock_dir / "longreads2.bam"
        lr_bam2.touch()

        pipeline = self.create_pipeline()
        pipeline.lr_aln = [self.mock_lr_bam, lr_bam2]
        pipeline.lr_readgroups = [
            [{"ID": "lr-rg1", "SM": "sample1"}],
            [{"ID": "lr-rg2", "SM": "sample1"}],
        ]
        pipeline.build_dag()

        assert pipeline.lr_alignment is None

    # The second, sex-aware DAG

    def _build_second_dag(self, pipeline, sample_sex=SampleSex.FEMALE):
        """Build both DAGs and return the second.

        The first DAG has to be built first so that the short-read
        alignments and the ploidy JSON are stashed for the second.
        """
        pipeline.build_dag()
        pipeline.sample_sex = sample_sex
        return pipeline.build_second_dag()

    def _second_dag_job(self, pipeline, name, sample_sex=SampleSex.FEMALE):
        """Build both DAGs and return the named job of the second"""
        dag = self._build_second_dag(pipeline, sample_sex)
        _, all_jobs = self._get_all_job_names(dag)
        return self._get_job(all_jobs, name)

    def enable_cnv(self, pipeline):
        """Call CNVs with a bundle CNV model and a PAR BED file"""
        par_bed = self.mock_dir / "par.bed"
        par_bed.write_text("chrX\t10000\t2781479\n")
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True
        pipeline.cnv_par_bed = par_bed
        return pipeline

    def script_path(self, name):
        """A path to a packaged helper script, as the stages build it"""
        return pathlib.Path(str(files("sentieon_cli.scripts").joinpath(name)))

    CNV_JOBS = ("cnvscope", "cnv-model-apply", "indel2cnv", "combine-cnv")

    def test_no_second_dag_without_callers(self):
        """No second DAG when nothing consumes the sample sex"""
        pipeline = self.create_pipeline()

        assert self._build_second_dag(pipeline) is None

    def test_no_second_dag_without_a_cnv_model(self):
        """`--call_cnvs` without a bundle CNV model is rejected earlier"""
        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = False

        assert self._build_second_dag(pipeline) is None

    def test_second_dag_for_each_caller(self):
        """Every sex-aware caller on its own builds the second DAG"""
        pipeline = self.create_pipeline()
        pipeline.expansion_catalog = self.mock_bed
        job_names, _ = self._get_all_job_names(
            self._build_second_dag(pipeline)
        )
        assert job_names == ["expansion-hunter"]

        pipeline = self.create_pipeline()
        pipeline.segdup_caller = []
        job_names, _ = self._get_all_job_names(
            self._build_second_dag(pipeline)
        )
        assert job_names == ["segdup-caller"]

        pipeline = self.enable_cnv(self.create_pipeline())
        job_names, _ = self._get_all_job_names(
            self._build_second_dag(pipeline)
        )
        assert sorted(job_names) == sorted(self.CNV_JOBS)

    def test_cnv_jobs_in_the_second_dag(self):
        """The CNV jobs and their dependencies"""
        pipeline = self.enable_cnv(self.create_pipeline())
        first_dag = pipeline.build_dag()

        # CNV calling is sex-aware, so it is not in the first DAG
        first_job_names, _ = self._get_all_job_names(first_dag)
        for name in self.CNV_JOBS:
            assert name not in first_job_names

        pipeline.sample_sex = SampleSex.FEMALE
        dag = pipeline.build_second_dag()
        job_names, all_jobs = self._get_all_job_names(dag)
        assert sorted(job_names) == sorted(self.CNV_JOBS)

        assert self._get_dep_names(dag, all_jobs, "cnvscope") == set()
        assert self._get_dep_names(dag, all_jobs, "indel2cnv") == set()
        assert self._get_dep_names(dag, all_jobs, "cnv-model-apply") == {
            "cnvscope"
        }
        assert self._get_dep_names(dag, all_jobs, "combine-cnv") == {
            "cnv-model-apply",
            "indel2cnv",
        }

        for name in self.CNV_JOBS:
            job = self._get_job(all_jobs, name)
            assert job.task_name == "cnv", name
        for name in ("cnvscope", "cnv-model-apply"):
            assert self._get_job(all_jobs, name).threads == pipeline.cores
        for name in ("indel2cnv", "combine-cnv"):
            assert self._get_job(all_jobs, name).threads == 0

    def test_cnvscope_command_with_fastq_input(self):
        """CNVscope runs on the deduplicated short reads"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(pipeline, "cnvscope")

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        probes = str(self.mock_vcf).replace(".vcf.gz", "_cnv.probes")
        assert str(job.shell) == (
            f"sentieon driver --input {bwa_aln} "
            f"--reference {self.mock_ref} --thread_count {pipeline.cores} "
            f"--interval {self.mock_bed} "
            f"--algo CNVscope --model {self.mock_bundle}/cnv.model "
            "--sex F "
            f"--dump_probes {probes} "
            f"{self.mock_dir}/sample-cnvscope.vcf.gz"
        )
        # The pipeline-generated alignment carries the LR attribute already
        assert "--replace_rg" not in str(job.shell)

    def test_cnvscope_command_with_aligned_input(self):
        """The `--sr_aln` input is used with its `--replace_rg` row"""
        pipeline = self.enable_cnv(self.create_aligned_pipeline())
        job = self._second_dag_job(pipeline, "cnvscope")

        probes = str(self.mock_vcf).replace(".vcf.gz", "_cnv.probes")
        assert str(job.shell) == (
            "sentieon driver "
            r"--replace_rg 'sr-rg1=ID:sr-rg1\tSM:sample1\tLR:0' "
            f"--input {self.mock_sr_bam} "
            f"--reference {self.mock_ref} --thread_count {pipeline.cores} "
            f"--interval {self.mock_bed} "
            f"--algo CNVscope --model {self.mock_bundle}/cnv.model "
            "--sex F "
            f"--dump_probes {probes} "
            f"{self.mock_dir}/sample-cnvscope.vcf.gz"
        )

    def test_cnvscope_without_a_bed(self):
        """Without `--bed`, CNVscope runs across the whole reference"""
        pipeline = self.enable_cnv(self.create_pipeline())
        pipeline.bed = None
        job = self._second_dag_job(pipeline, "cnvscope")

        assert "--interval" not in str(job.shell)

    def test_cnv_model_apply_command(self):
        """CNVModelApply filters the CNVscope output"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(pipeline, "cnv-model-apply")

        assert str(job.shell) == (
            f"sentieon driver --reference {self.mock_ref} "
            f"--thread_count {pipeline.cores} "
            f"--algo CNVModelApply --model {self.mock_bundle}/cnv.model "
            f"--vcf {self.mock_dir}/sample-cnvscope.vcf.gz "
            f"{self.mock_dir}/sample-cnv_model_apply.vcf.gz"
        )
        # The BED restricts CNVscope, not the model apply
        assert "--interval" not in str(job.shell)

    def test_cnvscope_male_sample(self):
        """A male sample is called with the PAR BED file"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(
            pipeline, "cnvscope", sample_sex=SampleSex.MALE
        )

        assert "--sex M" in str(job.shell)
        assert f"--par {pipeline.cnv_par_bed}" in str(job.shell)

    def test_cnvscope_female_sample(self):
        """A female sample does not need the PAR regions"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(
            pipeline, "cnvscope", sample_sex=SampleSex.FEMALE
        )

        assert "--sex F" in str(job.shell)
        assert "--par" not in str(job.shell)

    def test_cnvscope_unknown_sex(self, messages):
        """An unknown sex calls a diploid genome, with a warning"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(
            pipeline, "cnvscope", sample_sex=SampleSex.UNKNOWN
        )

        assert "--sex" not in str(job.shell)
        assert "--par" not in str(job.shell)
        assert any("diploid" in msg for msg in messages)

    def test_indel2cnv_command(self):
        """The PangenomeSV output is converted to CNVs"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(pipeline, "indel2cnv")

        sv_vcf = str(self.mock_vcf).replace(".vcf.gz", "_sv.vcf.gz")
        assert str(job.shell) == (
            f"{sys.executable} {self.script_path('indel2cnv.py')} "
            f"{self.mock_ref} {sv_vcf} "
            f"{self.mock_dir}/sample-sv_cnv.vcf.gz -t {pipeline.cores}"
        )

    def test_combine_cnv_command(self):
        """The CNV and converted SV calls are combined"""
        pipeline = self.enable_cnv(self.create_pipeline())
        job = self._second_dag_job(pipeline, "combine-cnv")

        cnv_vcf = str(self.mock_vcf).replace(".vcf.gz", "_cnv.vcf.gz")
        assert str(job.shell) == (
            f"{sys.executable} {self.script_path('combine_cnv.py')} "
            f"--cnv {self.mock_dir}/sample-cnv_model_apply.vcf.gz "
            f"--converted {self.mock_dir}/sample-sv_cnv.vcf.gz "
            f"-o {cnv_vcf}"
        )

    def test_expansion_hunter_command(self):
        """ExpansionHunter genotypes the short reads"""
        catalog = self.mock_dir / "catalog.json"
        catalog.touch()

        pipeline = self.create_pipeline()
        pipeline.expansion_catalog = catalog
        dag = self._build_second_dag(pipeline)
        _, all_jobs = self._get_all_job_names(dag)
        job = self._get_job(all_jobs, "expansion-hunter")

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        assert str(job.shell) == (
            f"ExpansionHunter --reads {bwa_aln} "
            f"--reference {self.mock_ref} "
            f"--variant-catalog {catalog} "
            "--sex female "
            f"--threads {pipeline.cores} "
            f"--output-prefix {self.mock_dir}/output_expansion"
        )
        assert job.task_name == "expansion-hunter"
        assert job.threads == pipeline.cores
        assert self._get_dep_names(dag, all_jobs, "expansion-hunter") == set()

    def test_segdup_command_with_fastq_input(self):
        """segdup-caller uses the short and long reads, and one bundle"""
        pipeline = self.create_pipeline()
        pipeline.segdup_caller = ["CFH", "CYP2D6", "SMN1"]
        first_job_names, _ = self._get_all_job_names(pipeline.build_dag())
        assert "segdup-caller" not in first_job_names

        pipeline.sample_sex = SampleSex.MALE
        dag = pipeline.build_second_dag()
        _, all_jobs = self._get_all_job_names(dag)
        job = self._get_job(all_jobs, "segdup-caller")

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        assert str(job.shell) == (
            f"segdup-caller --short {bwa_aln} "
            f"--long {self.mock_lr_bam} "
            f"--reference {self.mock_ref} "
            f"--sr_model {self.mock_bundle} "
            f"--lr_model {self.mock_bundle} "
            f"--input_vcf {self.mock_vcf} "
            "--sex male "
            "--genes CFH,CYP2D6,SMN1 "
            f"--outdir {self.mock_dir}/output_segdups"
        )
        assert job.task_name == "segdup"
        assert job.threads == pipeline.cores
        assert self._get_dep_names(dag, all_jobs, "segdup-caller") == set()

    def test_segdup_without_genes(self):
        """An empty gene list runs segdup-caller's own default set"""
        pipeline = self.create_pipeline()
        pipeline.segdup_caller = []
        job = self._second_dag_job(pipeline, "segdup-caller")

        assert "--genes" not in str(job.shell)
        assert "--sex female" in str(job.shell)

    def test_segdup_command_with_aligned_input(self):
        """The `--sr_aln` input is passed to segdup-caller as-is"""
        pipeline = self.create_aligned_pipeline()
        pipeline.segdup_caller = []
        job = self._second_dag_job(pipeline, "segdup-caller")

        assert f"--short {self.mock_sr_bam} " in str(job.shell)

    def test_segdup_command_with_realigned_long_reads(self):
        """The realigned long-read alignment reaches segdup-caller"""
        pipeline = self.create_lr_realign_pipeline()
        pipeline.segdup_caller = []
        job = self._second_dag_job(pipeline, "segdup-caller")

        lr_aln = str(self.mock_vcf).replace(".vcf.gz", "_mm2_sorted_0.cram")
        assert f"--long {lr_aln} " in str(job.shell)

    # Second-DAG gates

    @pytest.mark.parametrize(
        "call_cnvs,has_cnv_model,expected",
        [
            (False, False, False),
            (False, True, False),
            (True, False, False),
            (True, True, True),
        ],
    )
    def test_cnv_in_second_dag(self, call_cnvs, has_cnv_model, expected):
        """CNVs are called with `--call_cnvs` and a bundle CNV model"""
        pipeline = self.create_pipeline()
        pipeline.call_cnvs = call_cnvs
        pipeline.has_cnv_model = has_cnv_model

        assert pipeline._cnv_in_second_dag() is expected

    def test_needs_second_dag(self):
        """Each sex-aware caller requires the second DAG"""
        pipeline = self.create_pipeline()
        assert pipeline._needs_second_dag() is False

        pipeline = self.create_pipeline()
        pipeline.expansion_catalog = self.mock_bed
        assert pipeline._needs_second_dag() is True

        pipeline = self.create_pipeline()
        pipeline.segdup_caller = []
        assert pipeline._needs_second_dag() is True

        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True
        assert pipeline._needs_second_dag() is True

        # `--call_cnvs` without a CNV model is rejected during validation
        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        assert pipeline._needs_second_dag() is False

    # Optional caller validation

    def test_validate_segdup_single_input(self):
        """One short-read and one long-read alignment are accepted"""
        pipeline = self.create_aligned_pipeline()
        pipeline.segdup_caller = []

        pipeline.validate_segdup()  # no SystemExit

    def test_validate_segdup_two_sr_inputs(self):
        """segdup-caller takes a single `--sr_aln` file"""
        pipeline = self.create_aligned_pipeline()
        pipeline.segdup_caller = []
        pipeline.sr_aln = [self.mock_sr_bam, self.mock_sr_bam]

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_segdup()
        assert excinfo.value.code == 2

    def test_validate_segdup_two_lr_inputs(self):
        """segdup-caller takes a single `--lr_aln` file"""
        pipeline = self.create_pipeline()
        pipeline.segdup_caller = ["CFH"]
        pipeline.lr_aln = [self.mock_lr_bam, self.mock_lr_bam]

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_segdup()
        assert excinfo.value.code == 2

    def test_validate_segdup_skip_small_variants(self):
        """segdup-caller reads the small-variant VCF"""
        pipeline = self.create_pipeline()
        pipeline.segdup_caller = []
        pipeline.skip_small_variants = True

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_segdup()
        assert excinfo.value.code == 2

    def test_validate_segdup_without_the_argument(self):
        """Nothing is validated without `--segdup_caller`"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True

        pipeline.validate_segdup()  # no SystemExit

    def test_validate_expansion_two_sr_inputs(self):
        """ExpansionHunter takes a single `--sr_aln` file"""
        pipeline = self.create_aligned_pipeline()
        pipeline.expansion_catalog = self.mock_bed
        pipeline.sr_aln = [self.mock_sr_bam, self.mock_sr_bam]

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_expansion()
        assert excinfo.value.code == 2

    def test_validate_expansion_single_input(self):
        """A single short-read alignment is accepted"""
        pipeline = self.create_aligned_pipeline()
        pipeline.expansion_catalog = self.mock_bed

        pipeline.validate_expansion()  # no SystemExit

    def test_validate_cnv_requires_sv_calling(self):
        """indel2cnv reads the PangenomeSV output"""
        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True
        pipeline.skip_svs = True

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_cnv()
        assert excinfo.value.code == 2

    def test_validate_cnv_requires_a_cnv_model(self):
        """`--call_cnvs` needs a bundle with a 'cnv.model' file"""
        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = False

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_cnv()
        assert excinfo.value.code == 2

    def test_validate_cnv_requires_a_par_bed(self):
        # The mock reference is not a recognized build, so no packaged
        # PAR BED file can be selected and validation stops the run
        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_cnv()
        assert excinfo.value.code == 2

    def test_validate_cnv_accepts_the_par_bed_argument(self):
        par_bed = self.mock_dir / "par.bed"
        par_bed.write_text("chrX\t10000\t2781479\n")

        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True
        pipeline.par_bed = par_bed

        pipeline.validate_cnv()  # no SystemExit

        assert pipeline.cnv_par_bed == par_bed

    def test_validate_cnv_needs_no_par_bed_without_call_cnvs(self):
        pipeline = self.create_pipeline()
        pipeline.has_cnv_model = True

        pipeline.validate_cnv()  # no SystemExit

        assert pipeline.cnv_par_bed is None

    def test_validate_cnv_warns_without_a_bed(self):
        """CNVscope runs across every contig without `--bed`"""
        par_bed = self.mock_dir / "par.bed"
        par_bed.write_text("chrX\t10000\t2781479\n")

        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True
        pipeline.par_bed = par_bed
        pipeline.bed = None

        pipeline.validate_cnv()

        assert any(
            "No `--bed` supplied" in str(call.args[0])
            for call in pipeline.logger.warning.call_args_list
        )

    def test_validate_cnv_does_not_warn_with_a_bed(self):
        par_bed = self.mock_dir / "par.bed"
        par_bed.write_text("chrX\t10000\t2781479\n")

        pipeline = self.create_pipeline()
        pipeline.call_cnvs = True
        pipeline.has_cnv_model = True
        pipeline.par_bed = par_bed

        pipeline.validate_cnv()

        pipeline.logger.warning.assert_not_called()

    # Model bundle validation

    def bundle_pipeline(self, monkeypatch, members):
        """A pipeline whose bundle holds `members`"""
        pipeline = self.create_pipeline()
        pipeline.skip_pop_vcf_id_check = True

        def fake_ar_load(path):
            if str(path).endswith("bundle_info.json"):
                return b'{"pipeline": "Hybrid pangenome"}'
            return list(members)

        monkeypatch.setattr(hybrid_pangenome, "ar_load", fake_ar_load)
        return pipeline

    BUNDLE_MEMBERS = [
        "dnascope.model",
        "longreadsv.model",
        "minimap2.model",
        "extract.model",
        "bwa.model",
    ]

    def test_validate_bundle_sets_has_cnv_model(self, monkeypatch):
        pipeline = self.bundle_pipeline(
            monkeypatch, self.BUNDLE_MEMBERS + ["cnv.model"]
        )
        pipeline.validate_bundle()
        assert pipeline.has_cnv_model is True

        pipeline = self.bundle_pipeline(monkeypatch, self.BUNDLE_MEMBERS)
        pipeline.validate_bundle()
        assert pipeline.has_cnv_model is False

    def test_validate_bundle_requires_the_diploid_model(self, monkeypatch):
        """segdup-caller reads the bundle's `diploid_model`"""
        pipeline = self.bundle_pipeline(monkeypatch, self.BUNDLE_MEMBERS)
        pipeline.segdup_caller = []

        with pytest.raises(SystemExit) as excinfo:
            pipeline.validate_bundle()
        assert excinfo.value.code == 2

        pipeline = self.bundle_pipeline(
            monkeypatch, self.BUNDLE_MEMBERS + ["diploid_model"]
        )
        pipeline.segdup_caller = []
        pipeline.validate_bundle()  # no SystemExit

    def test_validate_bundle_diploid_model_not_required_by_default(
        self, monkeypatch
    ):
        pipeline = self.bundle_pipeline(monkeypatch, self.BUNDLE_MEMBERS)
        pipeline.validate_bundle()  # no SystemExit
