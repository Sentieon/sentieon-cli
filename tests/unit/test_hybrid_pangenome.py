"""
Unit tests for the HybridPangenome pipeline logic
"""

import os
import pathlib
import sys
import tempfile
from unittest.mock import MagicMock

import pytest

# Add the parent directory to the path to import sentieon_cli
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from sentieon_cli import command_strings as cmds
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
        assert self._get_dep_names(dag, all_jobs, "dnascope") == {
            "mm2-lift"
        }

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
