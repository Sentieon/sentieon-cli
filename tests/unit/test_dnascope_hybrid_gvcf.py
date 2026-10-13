"""
Unit tests for the DNAscope hybrid pipeline's `--gvcf` output
"""

from unittest.mock import MagicMock, patch

import pytest

from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline


@pytest.fixture
def pipeline(tmp_path):
    """A DNAscopeHybridPipeline ready for `build_dag`"""
    for name in [
        "reference.fa",
        "short.bam",
        "long.bam",
        "population.vcf.gz",
        "interval.bed",
    ]:
        (tmp_path / name).touch()
    (tmp_path / "model.bundle").touch()

    with patch("sys.exit"):
        pipeline = DNAscopeHybridPipeline()
    pipeline.logger = MagicMock()

    pipeline.output_vcf = tmp_path / "output.vcf.gz"
    pipeline.reference = tmp_path / "reference.fa"
    pipeline.model_bundle = tmp_path / "model.bundle"
    pipeline.bed = tmp_path / "interval.bed"
    pipeline.cores = 2
    pipeline.dry_run = True
    pipeline.skip_version_check = True
    pipeline.tmp_dir = tmp_path
    pipeline.pop_vcf = None
    pipeline.gvcf = True

    pipeline.fai_data = {"chr1": {"length": 1000}}
    shard = MagicMock()
    shard.contig = "chr1"
    shard.start = 1
    shard.stop = 1000
    pipeline.shards = [shard]
    pipeline.pop_vcf_contigs = {"chr1": 1000}

    pipeline.lr_aln_readgroups = [[{"ID": "lr_rg1", "SM": "sample1"}]]
    pipeline.sr_aln_readgroups = [[{"ID": "sr_rg1", "SM": "sample1"}]]
    pipeline.hybrid_rg_sm = "sample1"
    pipeline.hybrid_set_rg = False
    pipeline.shortread_tech = "Illumina"
    pipeline.longread_tech = "ONT"
    pipeline.sr_aln = [tmp_path / "short.bam"]
    pipeline.lr_aln = [tmp_path / "long.bam"]
    return pipeline


def _jobs(pipeline):
    """Build the DAG, returning {job name: (shell, dependency names)}"""
    dag = pipeline.build_dag()
    jobs = {}
    for job, deps in dag.waiting_jobs.items():
        assert job.name not in jobs, f"duplicate job name {job.name}"
        jobs[job.name] = (str(job.shell), {d.name for d in deps})
    return jobs


def _out(pipeline, suffix):
    """A path next to the output VCF"""
    return str(pipeline.output_vcf).replace(".vcf.gz", suffix)


class TestHybridGvcf:
    def test_calling_passes_emit_gvcfs(self, pipeline):
        jobs = _jobs(pipeline)
        tmp = pipeline.tmp_dir
        call1, _ = jobs["dnascope-1"]
        assert "--emit_mode gvcf" in call1
        assert call1.endswith(str(tmp / "initial.g.vcf.gz"))
        call2, _ = jobs["dnascope-2"]
        assert "--emit_mode gvcf" in call2
        assert call2.endswith(str(tmp / "hybrid_pass2.g.vcf.gz"))

    def test_both_passes_are_clipped_before_the_concat(self, pipeline):
        jobs = _jobs(pipeline)
        tmp = pipeline.tmp_dir
        stage2_bed = str(tmp / "hybrid_stage2.bed")

        pass1, pass1_deps = jobs["subset-calls"]
        assert "gvcf_clip.py" in pass1
        assert f"-v {tmp / 'initial.g.vcf.gz'}" in pass1
        assert f"-b {stage2_bed}" in pass1
        assert f"-r {pipeline.reference}" in pass1
        assert "--exclude" in pass1
        assert pass1.endswith(str(tmp / "mix_subset.g.vcf.gz"))
        assert pass1_deps == {"second-stage"}

        pass2, pass2_deps = jobs["subset-pass2-calls"]
        assert "gvcf_clip.py" in pass2
        assert f"-v {tmp / 'hybrid_pass2.g.vcf.gz'}" in pass2
        assert f"-b {stage2_bed}" in pass2
        assert "--exclude" not in pass2
        assert pass2.endswith(str(tmp / "hybrid_pass2_subset.g.vcf.gz"))
        assert pass2_deps == {"dnascope-2"}

        concat, concat_deps = jobs["concat-calls"]
        assert concat.endswith(
            f"{tmp / 'mix_subset.g.vcf.gz'} "
            f"{tmp / 'hybrid_pass2_subset.g.vcf.gz'}"
        )
        assert concat_deps == {"subset-calls", "subset-pass2-calls"}

        cleanup, _ = jobs["rm-tmp5"]
        assert str(tmp / "hybrid_pass2_subset.g.vcf.gz") in cleanup

    def test_model_apply_writes_the_output_gvcf(self, pipeline):
        jobs = _jobs(pipeline)
        apply, apply_deps = jobs["model-apply"]
        assert f"--vcf {pipeline.tmp_dir / 'combined_tmp_anno.g.vcf.gz'}" in (
            apply
        )
        assert apply.endswith(_out(pipeline, ".g.vcf.gz"))
        assert apply_deps == {"anno-calls"}

    def test_gvcftyper_then_field_removal_then_final_norm(self, pipeline):
        jobs = _jobs(pipeline)
        tmp = pipeline.tmp_dir
        gvcftyper_vcf = str(tmp / "combined_gvcftyper.vcf.gz")
        typed_vcf = str(tmp / "combined_typed.vcf.gz")

        typer, typer_deps = jobs["gvcftyper"]
        assert f"--interval {pipeline.bed}" in typer
        assert f"--vcf {_out(pipeline, '.g.vcf.gz')}" in typer
        assert "--dbsnp" not in typer
        assert typer.endswith(gvcftyper_vcf)
        assert typer_deps == {"model-apply"}

        drop, drop_deps = jobs["drop-allele-fields"]
        assert drop == (
            "bcftools annotate -x FORMAT/LAD,FORMAT/LPL,FORMAT/SAD,FORMAT/SPL "
            f"-Oz -W=tbi -o {typed_vcf} {gvcftyper_vcf}"
        )
        assert drop_deps == {"gvcftyper"}

        norm, norm_deps = jobs["final-norm"]
        assert f"""-e 'GT="0/0"' {typed_vcf}""" in norm
        assert norm.endswith(str(pipeline.output_vcf))
        assert norm_deps == {"drop-allele-fields"}

    def test_gvcftyper_annotates_dbsnp(self, pipeline):
        pipeline.dbsnp = pipeline.tmp_dir / "dbsnp.vcf.gz"
        jobs = _jobs(pipeline)
        typer, _ = jobs["gvcftyper"]
        assert f"--dbsnp {pipeline.dbsnp}" in typer

    def test_pop_vcf_transfer_runs_on_the_gvcf(self, pipeline):
        pipeline.pop_vcf = pipeline.tmp_dir / "population.vcf.gz"
        jobs = _jobs(pipeline)
        transfer_vcf = str(pipeline.tmp_dir / "combined_tmp_transfer.g.vcf.gz")

        concat, _ = jobs["merge-trim-concat"]
        assert f"--output {transfer_vcf}" in concat
        apply, apply_deps = jobs["model-apply"]
        assert f"--vcf {transfer_vcf}" in apply
        assert apply.endswith(_out(pipeline, ".g.vcf.gz"))
        assert apply_deps == {"merge-trim-concat"}
        assert jobs["gvcftyper"][1] == {"model-apply"}

    def test_skip_model_apply(self, pipeline):
        pipeline.skip_model_apply = True
        jobs = _jobs(pipeline)
        gvcftyper_vcf = str(pipeline.tmp_dir / "combined_gvcftyper.vcf.gz")

        anno, _ = jobs["anno-calls"]
        assert anno.endswith(_out(pipeline, ".g.vcf.gz"))
        typer, typer_deps = jobs["gvcftyper"]
        assert f"--vcf {_out(pipeline, '.g.vcf.gz')}" in typer
        assert typer.endswith(gvcftyper_vcf)
        assert typer_deps == {"anno-calls"}
        drop, drop_deps = jobs["drop-allele-fields"]
        assert drop.endswith(f"-o {pipeline.output_vcf} {gvcftyper_vcf}")
        assert drop_deps == {"gvcftyper"}
        assert "model-apply" not in jobs
        assert "final-norm" not in jobs

    def test_skip_model_apply_with_pop_vcf(self, pipeline):
        pipeline.skip_model_apply = True
        pipeline.pop_vcf = pipeline.tmp_dir / "population.vcf.gz"
        jobs = _jobs(pipeline)

        concat, _ = jobs["merge-trim-concat"]
        assert f"--output {_out(pipeline, '.g.vcf.gz')}" in concat
        assert jobs["gvcftyper"][1] == {"merge-trim-concat"}
        drop, _ = jobs["drop-allele-fields"]
        assert f"-o {pipeline.output_vcf} " in drop
        assert "final-norm" not in jobs


class TestHybridWithoutGvcf:
    def test_vcf_path(self, pipeline):
        pipeline.gvcf = False
        jobs = _jobs(pipeline)
        tmp = pipeline.tmp_dir

        assert "--emit_mode" not in jobs["dnascope-1"][0]
        assert "--emit_mode" not in jobs["dnascope-2"][0]
        subset, _ = jobs["subset-calls"]
        assert "bcftools view -T" in subset
        assert "gvcf_clip.py" not in subset
        assert "subset-pass2-calls" not in jobs
        assert jobs["concat-calls"][1] == {"subset-calls", "dnascope-2"}
        assert "gvcftyper" not in jobs
        assert "drop-allele-fields" not in jobs

        apply, _ = jobs["model-apply"]
        assert apply.endswith(str(tmp / "combined_apply.vcf.gz"))
        norm, norm_deps = jobs["final-norm"]
        assert f"""-e 'GT="0/0"' {tmp / 'combined_apply.vcf.gz'}""" in norm
        assert norm.endswith(str(pipeline.output_vcf))
        assert norm_deps == {"model-apply"}
        assert not any(".g.vcf.gz" in shell for shell, _ in jobs.values())
