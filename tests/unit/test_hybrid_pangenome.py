"""
Unit tests for the HybridPangenome pipeline logic
"""

import os
import pathlib
import sys
import tempfile
from unittest.mock import MagicMock

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli.hybrid_pangenome import HybridPangenome
from sentieon_cli.command_strings import LONGREAD_SV_BED_AWK
from sentieon_cli.dag import DAG


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

    def _get_all_job_names(self, dag):
        """Helper to get all job names from a DAG"""
        all_jobs = list(dag.waiting_jobs.keys()) + list(dag.ready_jobs.keys())
        return [job.name for job in all_jobs], all_jobs

    def _get_job(self, all_jobs, name):
        return next(j for j in all_jobs if j.name == name)

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
            "pangenome-sv",
            "dnascope-raw",
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
        assert "Rehead metrics" in job_names

    def test_skip_metrics(self):
        """skip_metrics removes the metrics jobs and the Dedup metrics"""
        pipeline = self.create_pipeline()
        pipeline.skip_metrics = True
        dag = pipeline.build_dag()
        job_names, all_jobs = self._get_all_job_names(dag)

        assert "metrics" not in job_names
        assert "Rehead metrics" not in job_names
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

        update_cmd = str(self._get_job(all_jobs, "graph-update").shell)
        assert "--algo PGHapUpdateAlgo" in update_cmd
        assert "--gfa_file " in update_cmd
        assert "sample-pangenome-raw.gfa" in update_cmd
        assert "--target_bed " in update_cmd
        assert "sample-sv.bed" in update_cmd

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
        assert "sort -k1,1 -k2,2n" in bed_cmd
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

    def test_calling_inputs_and_replace_rg(self):
        """PangenomeSV and DNAscope use the bwa, lifted, and long reads,
        rewriting only the long-read readgroups with LR:1"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        bwa_aln = str(self.mock_vcf).replace(".vcf.gz", "_bwa_deduped.cram")
        lift_aln = str(self.mock_vcf).replace(".vcf.gz", "_lift_deduped.cram")
        replace_arg = r"lr-rg1=ID:lr-rg1\tSM:sample1\tLR:1"

        for job_name in ("pangenome-sv", "dnascope-raw"):
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
        sv_vcf = str(self.mock_vcf).replace(".vcf.gz", "_sv.vcf.gz")
        assert sv_vcf in cmd_str
        # SV calling is not restricted to the small-variant BED
        assert f"--interval {self.mock_bed}" not in cmd_str

    def test_dnascope_command(self):
        """DNAscope runs with the model, interval, and pcr_indel_model"""
        pipeline = self.create_pipeline()
        dag = pipeline.build_dag()
        _, all_jobs = self._get_all_job_names(dag)

        cmd_str = str(self._get_job(all_jobs, "dnascope-raw").shell)
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

        cmd_str = str(self._get_job(all_jobs, "dnascope-raw").shell)
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

        cmd_str = str(self._get_job(all_jobs, "dnascope-raw").shell)
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
        assert "dnascope-raw" in job_names
        assert "graph-update" in job_names

    def test_skip_small_variants(self):
        """skip_small_variants removes DNAscope, transfer, and model-apply"""
        pipeline = self.create_pipeline()
        pipeline.skip_small_variants = True
        dag = pipeline.build_dag()
        job_names, _ = self._get_all_job_names(dag)

        assert "dnascope-raw" not in job_names
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

        cmd_str = str(self._get_job(all_jobs, "dnascope-raw").shell)
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

        cmd_str = str(self._get_job(all_jobs, "dnascope-raw").shell)
        assert r"lr-rg1=ID:lr-rg1\tSM:override_sm\tLR:1" in cmd_str

        bwa_cmd = str(self._get_job(all_jobs, "bwa-extract").shell)
        assert "SM:override_sm" in bwa_cmd
