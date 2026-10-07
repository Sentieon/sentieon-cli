"""
Unit tests for the DNAscope pipeline's model-bundle handling: the
`--pop_vcf` annotation transfer and the platform's read filter
"""

import json
import pathlib
from typing import Dict, List, Optional
from unittest.mock import MagicMock

import pytest

from sentieon_cli.dnascope import DNAscopePipeline

BUNDLE_VCF_ID = "pop-vcf-v1"
CHR1_LENGTH = 25 * 1000 * 1000


def write_bundle(path: pathlib.Path, bundle_info: Optional[Dict]) -> None:
    """Write an ar archive holding at most a bundle_info.json member"""
    data = b"!<arch>\n"
    if bundle_info is not None:
        body = json.dumps(bundle_info).encode()
        name = b"bundle_info.json"  # exactly 16 bytes, so no "/" suffix
        header = (
            name.ljust(16)
            + b"0".ljust(12)
            + b"0".ljust(6)
            + b"0".ljust(6)
            + b"644".ljust(8)
            + str(len(body)).encode().ljust(10)
            + b"`\n"
        )
        data += header + body + (b"\n" if len(body) % 2 else b"")
    path.write_bytes(data)


def write_pop_vcf(
    path: pathlib.Path,
    vcf_id: Optional[str],
    contigs: Dict[str, int],
) -> None:
    """Write a header-only population VCF"""
    lines = ["##fileformat=VCFv4.2"]
    if vcf_id is not None:
        lines.append(f"##SentieonVcfID={vcf_id}")
    for ctg, length in contigs.items():
        lines.append(f"##contig=<ID={ctg},length={length}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def files(tmp_path):
    """The reference, alignment, bundle and pop VCF of a run"""
    ref = tmp_path / "reference.fa"
    ref.touch()
    pathlib.Path(str(ref) + ".fai").write_text(
        f"chr1\t{CHR1_LENGTH}\t6\t60\t61\nchrM\t16569\t25416684\t60\t61\n"
    )
    bam = tmp_path / "sample.bam"
    bam.touch()
    bundle = tmp_path / "model.bundle"
    write_bundle(bundle, {"SentieonVcfID": BUNDLE_VCF_ID})
    pop_vcf = tmp_path / "pop.vcf.gz"
    write_pop_vcf(pop_vcf, BUNDLE_VCF_ID, {"chr1": CHR1_LENGTH})
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    return {
        "ref": ref,
        "bam": bam,
        "bundle": bundle,
        "pop_vcf": pop_vcf,
        "tmp_dir": tmp_dir,
        "output_vcf": tmp_path / "output.vcf.gz",
    }


def make_pipeline(files, **overrides) -> DNAscopePipeline:
    """A dry-run DNAscope pipeline over a single BAM input"""
    pipeline = DNAscopePipeline()
    pipeline.logger = MagicMock()
    pipeline.output_vcf = files["output_vcf"]
    pipeline.reference = files["ref"]
    pipeline.model_bundle = files["bundle"]
    pipeline.sample_input = [files["bam"]]
    pipeline.pop_vcf = files["pop_vcf"]
    pipeline.cores = 2
    pipeline.dry_run = True
    pipeline.skip_version_check = True
    pipeline.tmp_dir = files["tmp_dir"]
    for key, value in overrides.items():
        setattr(pipeline, key, value)
    return pipeline


def build(pipeline: DNAscopePipeline):
    """Validate, configure and build the DAG"""
    pipeline.validate()
    pipeline.configure()
    dag = pipeline.build_dag()
    return dag


def all_jobs(dag) -> List:
    return list(dag.waiting_jobs) + list(dag.ready_jobs)


def job_named(dag, name: str):
    matches = [job for job in all_jobs(dag) if job.name == name]
    assert len(matches) == 1, (name, [job.name for job in all_jobs(dag)])
    return matches[0]


def dep_names(dag, job) -> List[str]:
    return sorted(dep.name for dep in dag.waiting_jobs.get(job, set()))


def exit_with_error(pipeline: DNAscopePipeline, message: str) -> None:
    """`validate()` exits 2 and logs an error containing `message`"""
    with pytest.raises(SystemExit) as excinfo:
        pipeline.validate()
    assert excinfo.value.code == 2
    logged = [str(call.args[0]) for call in pipeline.logger.error.call_args_list]
    assert any(message in line for line in logged), logged


class TestPopVcfValidation:
    """The `--pop_vcf` is checked against the model bundle"""

    def test_matching_pop_vcf_is_accepted(self, files):
        pipeline = make_pipeline(files)
        pipeline.validate()
        pipeline.logger.error.assert_not_called()

    def test_bundle_requiring_a_pop_vcf_rejects_a_run_without_one(
        self, files
    ):
        pipeline = make_pipeline(files, pop_vcf=None)
        exit_with_error(pipeline, "requires a population VCF")

    def test_mismatched_pop_vcf_id_is_rejected(self, files):
        write_pop_vcf(files["pop_vcf"], "other-id", {"chr1": CHR1_LENGTH})
        pipeline = make_pipeline(files)
        exit_with_error(pipeline, "does not match the population VCF")

    def test_pop_vcf_without_an_id_is_rejected(self, files):
        write_pop_vcf(files["pop_vcf"], None, {"chr1": CHR1_LENGTH})
        pipeline = make_pipeline(files)
        exit_with_error(pipeline, "does not match the population VCF")

    def test_skip_pop_vcf_id_check_accepts_a_mismatch(self, files):
        write_pop_vcf(files["pop_vcf"], "other-id", {"chr1": CHR1_LENGTH})
        pipeline = make_pipeline(files, skip_pop_vcf_id_check=True)
        pipeline.validate()
        pipeline.logger.error.assert_not_called()

    @pytest.mark.parametrize("bundle_info", [None, {"minScriptVersion": "1"}])
    def test_bundle_without_an_id_rejects_a_pop_vcf(self, files, bundle_info):
        write_bundle(files["bundle"], bundle_info)
        pipeline = make_pipeline(files)
        exit_with_error(pipeline, "does not require a population VCF")

    @pytest.mark.parametrize("bundle_info", [None, {"minScriptVersion": "1"}])
    def test_bundle_without_an_id_runs_without_a_pop_vcf(
        self, files, bundle_info
    ):
        write_bundle(files["bundle"], bundle_info)
        pipeline = make_pipeline(files, pop_vcf=None)
        pipeline.validate()
        pipeline.logger.error.assert_not_called()

    def test_transfer_inputs_are_collected(self, files):
        pipeline = make_pipeline(files)
        pipeline.validate()
        assert pipeline.pop_vcf_contigs == {"chr1": CHR1_LENGTH}
        assert set(pipeline.fai_data) == {"chr1", "chrM"}
        # 10 Mb shards: three for chr1, one for chrM
        assert [s.contig for s in pipeline.shards] == ["chr1"] * 3 + ["chrM"]

    def test_no_reference_index_is_read_without_a_pop_vcf(self, files):
        write_bundle(files["bundle"], None)
        pathlib.Path(str(files["ref"]) + ".fai").unlink()
        pipeline = make_pipeline(files, pop_vcf=None)
        pipeline.validate()
        assert pipeline.fai_data == {}
        assert pipeline.shards == []


class TestPopVcfDag:
    """Annotations are transferred before DNAModelApply"""

    def test_transfer_runs_between_dnascope_and_model_apply(self, files):
        dag = build(make_pipeline(files))

        names = [job.name for job in all_jobs(dag)]
        assert sorted(n for n in names if n.startswith("merge-trim")) == [
            "merge-trim-0",
            "merge-trim-1",
            "merge-trim-2",
            "merge-trim-concat",
            # chrM is not in the pop VCF, so it is copied, not merged
            "merge-trim-extra",
        ]
        for name in ("merge-trim-0", "merge-trim-extra"):
            assert dep_names(dag, job_named(dag, name)) == ["dnascope"]

        apply_job = job_named(dag, "model-apply")
        assert dep_names(dag, apply_job) == ["merge-trim-concat"]
        transfer_vcf = files["tmp_dir"] / "sample-dnascope_transfer.vcf.gz"
        assert str(transfer_vcf) in str(job_named(dag, "merge-trim-concat").shell)
        assert str(transfer_vcf) in str(apply_job.shell)
        assert str(files["output_vcf"]) in str(apply_job.shell)

        assert dep_names(dag, job_named(dag, "rm-tmp-vcf")) == ["model-apply"]

    def test_no_transfer_without_a_pop_vcf(self, files):
        write_bundle(files["bundle"], None)
        dag = build(make_pipeline(files, pop_vcf=None))

        names = [job.name for job in all_jobs(dag)]
        assert not any(n.startswith("merge-trim") for n in names)
        apply_job = job_named(dag, "model-apply")
        assert dep_names(dag, apply_job) == ["dnascope"]

    def test_gvcf_transfer(self, files):
        dag = build(make_pipeline(files, gvcf=True))

        transfer_gvcf = files["tmp_dir"] / "sample-dnascope_transfer.g.vcf.gz"
        apply_job = job_named(dag, "model-apply")
        assert str(transfer_gvcf) in str(apply_job.shell)
        out_gvcf = str(files["output_vcf"]).replace(".vcf.gz", ".g.vcf.gz")
        assert out_gvcf in str(apply_job.shell)
        assert dep_names(dag, job_named(dag, "gvcftyper")) == ["model-apply"]

    def test_skip_small_variants_has_no_transfer(self, files):
        dag = build(make_pipeline(files, skip_small_variants=True))
        names = [job.name for job in all_jobs(dag)]
        assert not any(n.startswith("merge-trim") for n in names)


class TestUltimaReadFilter:
    """An Ultima bundle adds `UltimaReadFilter` to the DNAscope call"""

    @staticmethod
    def build_with_bundle(files, bundle_info, **overrides):
        write_bundle(files["bundle"], bundle_info)
        return build(make_pipeline(files, pop_vcf=None, **overrides))

    @pytest.mark.parametrize("platform", ["Ultima", "ULTIMA", "ultima"])
    def test_ultima_bundle_filters_the_dnascope_call(self, files, platform):
        dag = self.build_with_bundle(files, {"platform": platform})

        cmd = str(job_named(dag, "dnascope").shell)
        assert cmd.count("--read_filter UltimaReadFilter") == 1
        # Duplicate marking and the rest of the pipeline are unfiltered
        others = [j for j in all_jobs(dag) if j.name != "dnascope"]
        assert others
        for job in others:
            assert "UltimaReadFilter" not in str(job.shell), job.name

    def test_ultima_filter_covers_sv_calling(self, files):
        dag = self.build_with_bundle(files, {"platform": "Ultima"})

        # Small variants and SVs come from the same driver call
        cmd = str(job_named(dag, "dnascope").shell)
        assert "--var_type BND" in cmd
        assert "--read_filter UltimaReadFilter" in cmd

    def test_ultima_bundle_with_a_pop_vcf(self, files):
        write_bundle(
            files["bundle"],
            {"platform": "Ultima", "SentieonVcfID": BUNDLE_VCF_ID},
        )
        dag = build(make_pipeline(files))

        cmd = str(job_named(dag, "dnascope").shell)
        assert "--read_filter UltimaReadFilter" in cmd
        assert dep_names(dag, job_named(dag, "model-apply")) == [
            "merge-trim-concat"
        ]

    @pytest.mark.parametrize(
        "bundle_info", [None, {}, {"platform": "Illumina"}]
    )
    def test_other_bundles_are_unfiltered(self, files, bundle_info):
        dag = self.build_with_bundle(files, bundle_info)

        for job in all_jobs(dag):
            assert "--read_filter" not in str(job.shell), job.name

    def test_platform_is_read_from_the_bundle(self, files):
        write_bundle(files["bundle"], {"platform": "Ultima"})
        pipeline = make_pipeline(files, pop_vcf=None)
        pipeline.validate()
        assert pipeline.tech == "ULTIMA"

        write_bundle(files["bundle"], None)
        pipeline = make_pipeline(files, pop_vcf=None)
        pipeline.validate()
        assert pipeline.tech == ""
