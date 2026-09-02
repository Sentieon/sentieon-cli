"""
The non-obvious edges of the DNAscope LongRead small-variant DAG
"""

import json
import pathlib
from typing import Dict, List, Set
from unittest.mock import patch

from sentieon_cli.dag import DAG
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.job import Job

from tests.utils.test_helpers import create_mock_args

BUNDLE_MEMBERS = [
    "diploid_hp_model",
    "diploid_model",
    "diploid_model_unphased",
    "gvcf_model",
    "haploid_hp_model",
    "haploid_model",
    "longreadsv.model",
    "minimap2.model",
]

BUNDLE_INFO = {
    "platform": "HiFi",
    "minScriptVersion": "1.5.2",
    "pipeline": "DNAscope LongRead",
}


def build_dag(tmp_path: pathlib.Path, **kwargs) -> DAG:
    """Build the first DAG of a HiFi long-read run"""
    reference = tmp_path / "reference.fa"
    reference.touch()
    (tmp_path / "reference.fa.fai").write_text(
        "chr1\t1000\t10\t60\t61\nchr2\t2000\t1030\t60\t61\n"
    )
    bam = tmp_path / "longread.bam"
    bam.touch()
    bundle = tmp_path / "model.bundle"
    bundle.touch()
    bed = tmp_path / "regions.bed"
    bed.touch()

    pipeline = DNAscopeLRPipeline()
    pipeline.setup_logging(create_mock_args())
    pipeline.output_vcf = tmp_path / "output.vcf.gz"
    pipeline.reference = reference
    pipeline.model_bundle = bundle
    pipeline.sample_input = [bam]
    pipeline.fastq = []
    pipeline.readgroups = []
    pipeline.bed = bed
    pipeline.cores = 2
    pipeline.dry_run = True
    pipeline.skip_version_check = True
    pipeline.tmp_dir = tmp_path
    # Only the small-variant jobs are under test
    pipeline.skip_mosdepth = True
    pipeline.skip_cnv = True
    pipeline.skip_svs = True
    for key, value in kwargs.items():
        setattr(pipeline, key, value)

    bundle_info = dict(BUNDLE_INFO)
    if pipeline.pop_vcf:
        bundle_info["SentieonVcfID"] = "test-pop-vcf"
        pipeline.skip_pop_vcf_id_check = True

    with patch("sentieon_cli.dnascope_longread.ar_load") as mock_ar_load:
        mock_ar_load.side_effect = [
            BUNDLE_MEMBERS,
            json.dumps(bundle_info).encode(),
        ]
        pipeline.validate()
    pipeline.configure()
    return pipeline.build_dag()


def deps_of(dag: DAG, name: str) -> List[str]:
    """The sorted dependency names of the single job called `name`"""
    return sorted(job.name for job in dag.waiting_jobs[job_named(dag, name)])


def job_named(dag: DAG, name: str) -> Job:
    """The one job called `name`"""
    matches = [job for job in all_jobs(dag) if job.name == name]
    assert len(matches) == 1, f"{name}: {len(matches)} jobs"
    return matches[0]


def all_jobs(dag: DAG) -> List[Job]:
    return list(dag.waiting_jobs) + list(dag.ready_jobs)


def roots(dag: DAG) -> Set[str]:
    """The names of the jobs with no dependencies"""
    names = {job.name for job in dag.ready_jobs}
    names.update(
        job.name for job, deps in dag.waiting_jobs.items() if not deps
    )
    return names


class TestSmallVariantEdges:
    """The edges `add_small_variant_calling` wires by hand"""

    def test_the_calling_passes_follow_the_phaser(self, tmp_path):
        dag = build_dag(tmp_path)

        assert deps_of(dag, "variantphaser") == ["model-apply-diploid"]
        for hap in ("dnascope-hap1", "dnascope-hap2"):
            assert deps_of(dag, hap) == ["repeatmodel", "variantphaser"]
        assert deps_of(dag, "dnascope-unphased") == [
            "bedtools-subtract",
            "repeatmodel",
        ]
        assert deps_of(dag, "patch") == ["dnascope-hap1", "dnascope-hap2"]

    def test_merge_waits_on_the_three_model_apply_jobs(self, tmp_path):
        dag = build_dag(tmp_path)

        assert deps_of(dag, "merge") == [
            "model-apply-hap1",
            "model-apply-hap2",
            "model-apply-unphased",
        ]

    def test_both_haplotype_applies_wait_on_both_transfers(self, tmp_path):
        pop_vcf = tmp_path / "population.vcf.gz"
        pop_vcf.touch()
        dag = build_dag(tmp_path, pop_vcf=pop_vcf)

        # Each model-apply reads only its own haplotype's transfer, but
        # both wait on both concat jobs, as the hand-wired DAG did
        for hap in ("model-apply-hap1", "model-apply-hap2"):
            assert deps_of(dag, hap) == [
                "merge-trim-hap1-concat",
                "merge-trim-hap2-concat",
                "patch",
            ]

    def test_each_transfer_is_tagged(self, tmp_path):
        pop_vcf = tmp_path / "population.vcf.gz"
        pop_vcf.touch()
        dag = build_dag(tmp_path, pop_vcf=pop_vcf)

        concats = {
            job.name for job in all_jobs(dag) if job.name.endswith("-concat")
        }
        assert concats == {
            "merge-trim-diploid-concat",
            "merge-trim-hap1-concat",
            "merge-trim-hap2-concat",
            "merge-trim-unphased-concat",
        }
        assert deps_of(dag, "model-apply-diploid") == [
            "dnascope-diploid",
            "merge-trim-diploid-concat",
        ]
        assert deps_of(dag, "model-apply-unphased") == [
            "diploid-patch",
            "merge-trim-unphased-concat",
        ]

    def test_fai_to_bed_is_a_root_without_a_bed(self, tmp_path):
        dag = build_dag(tmp_path, bed=None)

        assert "fai-to-bed" in roots(dag)
        assert deps_of(dag, "bedtools-subtract") == [
            "fai-to-bed",
            "variantphaser",
        ]

    def test_haploid_calling_only_follows_the_repeat_model(self, tmp_path):
        haploid_bed = tmp_path / "haploid.bed"
        haploid_bed.touch()
        dag = build_dag(tmp_path, haploid_bed=haploid_bed)

        # It does not wait on the alignment jobs, only on the repeat model
        assert deps_of(dag, "dnascope-haploid") == ["repeatmodel"]
        assert deps_of(dag, "haploid-patch2") == ["dnascope-haploid"]
        assert deps_of(dag, "haploid-diploid-concat") == [
            "haploid-patch2",
            "merge",
        ]

    def test_haploid_calling_is_a_root_with_a_repeat_model(self, tmp_path):
        haploid_bed = tmp_path / "haploid.bed"
        haploid_bed.touch()
        repeat_model = tmp_path / "given_repeat.model"
        repeat_model.touch()
        dag = build_dag(
            tmp_path, haploid_bed=haploid_bed, repeat_model=repeat_model
        )

        assert "repeatmodel" not in {job.name for job in all_jobs(dag)}
        assert "dnascope-haploid" in roots(dag)


class TestJobIdentity:
    """Job names and task names of the migrated stages"""

    def test_task_names(self, tmp_path):
        dag = build_dag(tmp_path)

        tasks: Dict[str, str] = {
            job.name: job.task_name for job in all_jobs(dag)
        }
        assert tasks["dnascope-diploid"] == "variant-calling"
        assert tasks["dnascope-hap1"] == "variant-calling"
        assert tasks["dnascope-unphased"] == "variant-calling"
        assert tasks["model-apply-diploid"] == "model-apply"
        assert tasks["model-apply-hap2"] == "model-apply"
        assert tasks["repeatmodel"] == "repeat-model"
