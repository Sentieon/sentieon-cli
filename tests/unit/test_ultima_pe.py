"""
Unit tests for the UltimaPE flow
"""

import json
import pathlib
from typing import Dict, List, Optional
from unittest.mock import MagicMock

import packaging.version
import pytest

from sentieon_cli import util
from sentieon_cli.dag import DAG
from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.driver import Driver, LuaScript
from sentieon_cli.sentieon_pangenome import SentieonPangenome
from sentieon_cli.shard import Shard, determine_shards_from_fai
from sentieon_cli.stages.base import StageContext, driver_job, rm_job
from sentieon_cli.stages.small_variants import (
    ApplySpec,
    PostSpec,
    TransferApplyStage,
    TransferSpec,
)
from sentieon_cli.stages.transfer import TransferConfig
from sentieon_cli.stages.ultima_pe import (
    ULTIMA_PE_MIN_VERSIONS,
    HapTableStage,
    PrimaryIndexStage,
    haplen_intervals,
    primary_contigs,
)
from tests.utils.test_helpers import write_bundle

VCF_ID = "pop-vcf-v1"
READ_FILTER = "UltimaReadFilter,remove_t0_tp=0,keep_oq=true"
REPEAT_MODEL_ARGS = (
    "--algo RepeatModel --min_map_qual 1 --min_group_count 5000 "
    "--read_flag_mask drop=supplementary --repeat_extension 5 "
    "--max_repeat_unit_size 2 --min_repeat_count 6"
)
MEMBERS = ["dnascope.model", "read_tables.tsv", "ultima.model"]
PANGENOME_MEMBERS = MEMBERS + ["extract.model", "minimap2.model"]

# A GRCh38-like index: primary, alt, HLA, decoy and EBV contigs
FAI = {
    "chr1": 25 * 1000 * 1000,
    "chr1_KI270762v1_alt": 354444,
    "chr2": 1000,
    "chrX": 2000,
    "chrY": 3000,
    "chrM": 16569,
    "chrUn_KI270302v1": 2274,
    "chrEBV": 171823,
    "chrUn_JTFH01000001v1_decoy": 25139,
    "HLA-A*01:01:01:01": 3503,
}


def write_ultima_pe_bundle(
    path: pathlib.Path,
    pipeline: str = "DNAscope",
    platform: str = "UltimaPE",
    members: Optional[List[str]] = None,
    vcf_id: Optional[str] = VCF_ID,
) -> None:
    """A bundle with a bundle_info.json and empty model members"""
    info: Dict[str, str] = {
        "minScriptVersion": "2.0.1",
        "pipeline": pipeline,
        "platform": platform,
    }
    if vcf_id:
        info["SentieonVcfID"] = vcf_id
    write_bundle(
        path,
        {
            "bundle_info.json": json.dumps(info).encode(),
            **{m: b"" for m in (MEMBERS if members is None else members)},
        },
    )


def write_pop_vcf(path: pathlib.Path) -> None:
    """A header-only population VCF with the bundle's ID"""
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"##SentieonVcfID={VCF_ID}\n"
        f"##contig=<ID=chr1,length={FAI['chr1']}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )


@pytest.fixture
def files(tmp_path):
    """The inputs of a run"""
    ref = tmp_path / "reference.fa"
    ref.touch()
    pathlib.Path(str(ref) + ".fai").write_text(
        "".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in FAI.items())
    )
    bam = tmp_path / "sample.bam"
    bam.touch()
    bundle = tmp_path / "model.bundle"
    write_ultima_pe_bundle(bundle)
    pop_vcf = tmp_path / "pop.vcf.gz"
    write_pop_vcf(pop_vcf)
    bed = tmp_path / "calls.bed"
    bed.write_text("chr1\t0\t1000\n")
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    return {
        "ref": ref,
        "bam": bam,
        "bundle": bundle,
        "pop_vcf": pop_vcf,
        "bed": bed,
        "tmp": tmp_dir,
        "out": tmp_path / "output.vcf.gz",
    }


def all_jobs(dag) -> List:
    return list(dag.waiting_jobs) + list(dag.ready_jobs)


def job_named(dag, name: str):
    matches = [job for job in all_jobs(dag) if job.name == name]
    assert len(matches) == 1, (name, [job.name for job in all_jobs(dag)])
    return matches[0]


def dep_names(dag, job) -> List[str]:
    return sorted(dep.name for dep in dag.waiting_jobs.get(job, set()))


def logged(mock_method) -> List[str]:
    """The messages logged through a mocked logger method, formatted"""
    return [
        str(call.args[0]) % call.args[1:]
        for call in mock_method.call_args_list
    ]


def exits_with(fn, logger, message: str) -> None:
    """`fn()` exits 2 and logs an error containing `message`"""
    with pytest.raises(SystemExit) as excinfo:
        fn()
    assert excinfo.value.code == 2
    assert any(message in line for line in logged(logger.error)), logged(
        logger.error
    )


# --- The linear pipeline -------------------------------------------------


def linear_pipeline(files, **overrides) -> DNAscopePipeline:
    """A dry-run DNAscope pipeline over a single BAM, with a BED"""
    pipeline = DNAscopePipeline()
    pipeline.logger = MagicMock()
    pipeline.output_vcf = files["out"]
    pipeline.reference = files["ref"]
    pipeline.model_bundle = files["bundle"]
    pipeline.sample_input = [files["bam"]]
    pipeline.pop_vcf = files["pop_vcf"]
    pipeline.bed = files["bed"]
    pipeline.pcr_free = True
    pipeline.duplicate_marking = "none"
    pipeline.skip_metrics = True
    pipeline.skip_svs = True
    pipeline.cores = 2
    pipeline.dry_run = True
    pipeline.skip_version_check = True
    pipeline.tmp_dir = files["tmp"]
    for key, value in overrides.items():
        setattr(pipeline, key, value)
    return pipeline


def build(pipeline):
    pipeline.validate()
    pipeline.configure()
    return pipeline.build_dag()


class TestLinearDag:
    """DNAscope with d1.lua's table and RepeatModel, then post.lua"""

    def test_dnascope_command(self, files):
        dag = build(linear_pipeline(files))

        b, t = files["bundle"], files["tmp"]
        assert str(job_named(dag, "dnascope").shell) == (
            f"env HMGL_TABLES={b} sentieon driver --input {files['bam']} "
            f"--reference {files['ref']} --thread_count 2 "
            f"--interval {files['bed']} --read_filter {READ_FILTER} "
            f"--algo DNAscope --emit_mode variant --model {b}/dnascope.model "
            f"--pcr_indel_model NONE {str(files['out'])[:-7]}_tmp.vcf.gz "
            f"{REPEAT_MODEL_ARGS} {t}/sample-repeat.model"
        )

    def test_env_is_a_prefix_not_a_replacement(self, files):
        dag = build(linear_pipeline(files))

        (command,) = job_named(dag, "dnascope").shell.nodes
        assert command.executable == "env"
        assert command.exec_kwargs == {}

    def test_post_lua_between_transfer_and_apply(self, files):
        dag = build(linear_pipeline(files))

        b, t = files["bundle"], files["tmp"]
        post = job_named(dag, "post-lua")
        assert str(post.shell) == (
            f"sentieon driver --reference {files['ref']} --thread_count 2 "
            f"--interval {files['bed']} --algo Lua {b}/ultima.model "
            f"--repeat-model {t}/sample-repeat.model --lambda 2 --nofilt "
            f"{t}/sample-dnascope_transfer.vcf.gz "
            f"{t}/sample-dnascope_nofilt.vcf.gz"
        )
        assert post.threads == 2
        assert dep_names(dag, post) == ["merge-trim-concat"]

        apply_job = job_named(dag, "model-apply")
        assert str(apply_job.shell) == (
            f"sentieon driver --reference {files['ref']} --thread_count 2 "
            f"--algo DNAModelApply --model {b}/dnascope.model "
            f"--vcf {t}/sample-dnascope_nofilt.vcf.gz {files['out']}"
        )
        assert dep_names(dag, apply_job) == ["post-lua"]
        assert dep_names(dag, job_named(dag, "rm-tmp-vcf")) == ["model-apply"]

    def test_interval_padding_reaches_post_lua(self, files):
        dag = build(linear_pipeline(files, interval_padding=100))

        for name in ("dnascope", "post-lua"):
            cmd = str(job_named(dag, name).shell)
            assert f"--interval {files['bed']} --interval_padding 100" in cmd

    def test_without_a_bed(self, files):
        pipeline = linear_pipeline(files, bed=None)
        dag = build(pipeline)

        assert any(
            "the repeat model also counts the alt and decoy" in line
            for line in logged(pipeline.logger.warning)
        )
        for name in ("dnascope", "post-lua"):
            assert "--interval" not in str(job_named(dag, name).shell)

    def test_sv_calling_follows_repeat_model(self, files):
        dag = build(linear_pipeline(files, skip_svs=False))

        cmd = str(job_named(dag, "dnascope").shell)
        repeat = cmd.index("--algo RepeatModel")
        assert cmd.index("--algo DNAscope") < repeat
        assert cmd.index("--var_type BND") > repeat
        assert "UltimaReadFilter" not in str(job_named(dag, "svsolver").shell)

    def test_without_a_pop_vcf_post_lua_reads_the_raw_vcf(self, files):
        write_ultima_pe_bundle(files["bundle"], vcf_id=None)
        dag = build(linear_pipeline(files, pop_vcf=None))

        post = job_named(dag, "post-lua")
        raw = str(files["out"])[:-7] + "_tmp.vcf.gz"
        assert f"--nofilt {raw} " in str(post.shell)
        assert dep_names(dag, post) == ["dnascope"]
        assert dep_names(dag, job_named(dag, "model-apply")) == ["post-lua"]

    @pytest.mark.parametrize("platform", ["UltimaPE", "ULTIMAPE", "ultimape"])
    def test_platform_is_case_insensitive(self, files, platform):
        write_ultima_pe_bundle(files["bundle"], platform=platform)
        dag = build(linear_pipeline(files))

        assert "post-lua" in [job.name for job in all_jobs(dag)]

    @pytest.mark.parametrize("platform", ["Ultima", "Illumina", "UltimaPE2"])
    def test_other_platforms_are_unchanged(self, files, platform):
        write_ultima_pe_bundle(files["bundle"], platform=platform)
        dag = build(linear_pipeline(files))

        names = [job.name for job in all_jobs(dag)]
        assert "post-lua" not in names
        cmd = str(job_named(dag, "dnascope").shell)
        assert cmd.startswith("sentieon driver ")
        assert "RepeatModel" not in cmd
        assert "remove_t0_tp" not in cmd
        assert "sample-dnascope_transfer.vcf.gz" in str(
            job_named(dag, "model-apply").shell
        )


class TestLinearValidation:
    """An UltimaPE bundle must fit the run"""

    def test_a_pangenome_bundle_is_refused(self, files):
        write_ultima_pe_bundle(files["bundle"], pipeline="Sentieon pangenome")
        pipeline = linear_pipeline(files)
        exits_with(pipeline.validate, pipeline.logger, "different pipeline")

    @pytest.mark.parametrize("member", MEMBERS)
    def test_a_missing_member_is_refused(self, files, member):
        write_ultima_pe_bundle(
            files["bundle"], members=[m for m in MEMBERS if m != member]
        )
        pipeline = linear_pipeline(files)
        exits_with(pipeline.validate, pipeline.logger, member)

    def test_gvcf_is_refused(self, files):
        pipeline = linear_pipeline(files, gvcf=True)
        exits_with(pipeline.validate, pipeline.logger, "`--gvcf`")

    def test_skip_small_variants_needs_no_flow_members(self, files):
        write_ultima_pe_bundle(files["bundle"], members=[])
        pipeline = linear_pipeline(files, skip_small_variants=True, gvcf=True)
        pipeline.validate()
        pipeline.logger.error.assert_not_called()

    def test_the_driver_version_is_checked(self, files, monkeypatch):
        checked = []

        def check_version(cmd, version):
            checked.append((cmd, version))
            return False

        monkeypatch.setattr(util, "check_version", check_version)
        pipeline = linear_pipeline(files, skip_version_check=False)
        with pytest.raises(SystemExit):
            pipeline.validate()
        assert checked == [
            ("sentieon driver", packaging.version.Version("202503.05rc4"))
        ]


class TestMinimumVersion:
    """rc4 is a pre-release of 202503.05: it and every later build pass"""

    @pytest.mark.parametrize(
        "build,ok",
        [
            ("202503.05.rc4", True),
            ("202503.05.rc5", True),
            ("202503.05", True),
            ("202503.06", True),
            ("202509", True),
            ("202503.05.rc3", False),
            ("202503.04", False),
        ],
    )
    def test_sentieon_driver(self, monkeypatch, build, ok):
        monkeypatch.setattr(
            util.sp,
            "check_output",
            MagicMock(return_value=f"sentieon-genomics-{build}\n".encode()),
        )
        monkeypatch.setattr(util.shutil, "which", lambda cmd: "/bin/" + cmd)
        (min_version,) = ULTIMA_PE_MIN_VERSIONS.values()
        assert util.check_version("sentieon driver", min_version) is ok


# --- The pangenome pipeline ----------------------------------------------


def pangenome_pipeline(files, **overrides) -> SentieonPangenome:
    """A dry-run pangenome pipeline over a BAM, set up as validate() would"""
    for name in ("pangenome.gbz", "pangenome.hapl"):
        (files["tmp"].parent / name).touch()
    write_ultima_pe_bundle(
        files["bundle"],
        pipeline="Sentieon pangenome",
        members=PANGENOME_MEMBERS,
    )
    pipeline = SentieonPangenome()
    pipeline.logger = MagicMock()
    pipeline.output_vcf = files["out"]
    pipeline.reference = files["ref"]
    pipeline.model_bundle = files["bundle"]
    pipeline.sample_input = [files["bam"]]
    pipeline.gbz = files["tmp"].parent / "pangenome.gbz"
    pipeline.hapl = files["tmp"].parent / "pangenome.hapl"
    pipeline.pop_vcf = files["pop_vcf"]
    pipeline.bed = files["bed"]
    pipeline.cores = 2
    pipeline.dry_run = True
    pipeline.skip_version_check = True
    pipeline.skip_pop_vcf_id_check = True
    pipeline.tmp_dir = files["tmp"]
    pipeline.pangenome_ref_name = "GRCh38"
    pipeline.pangenome_contig_prefix = "GRCh38#0#"
    pipeline.fai_data = {c: {"length": n} for c, n in FAI.items()}
    pipeline.shards = determine_shards_from_fai(pipeline.fai_data, 10**7)
    pipeline.pop_vcf_contigs = {"chr1": FAI["chr1"]}
    pipeline.bam_readgroups = [{"ID": "rg1", "SM": "sample1"}]
    pipeline.fastq_readgroup = {}
    pipeline.tech = "ULTIMAPE"
    for key, value in overrides.items():
        setattr(pipeline, key, value)
    return pipeline


class TestPangenomeDag:
    """The calling pass, the hap-length table and post.lua"""

    def test_dnascope_command(self, files):
        dag = pangenome_pipeline(files).build_dag()

        b, t = files["bundle"], files["tmp"]
        lift = str(files["out"])[:-7] + "_mm2_deduped.cram"
        assert str(job_named(dag, "dnascope").shell) == (
            f"env HMGL_TABLES={b} sentieon driver --input {files['bam']} "
            f"--input {lift} --reference {files['ref']} --thread_count 2 "
            f"--interval {files['bed']} --read_filter {READ_FILTER} "
            f"--algo DNAscope --emit_mode variant --model {b}/dnascope.model "
            f"--pcr_indel_model CONSERVATIVE {t}/sample-dnascope.vcf.gz "
            f"{REPEAT_MODEL_ARGS} {t}/sample-repeat.model"
        )

    def test_hap_table_commands(self, files):
        dag = pangenome_pipeline(files).build_dag()

        b, t, ref = files["bundle"], files["tmp"], files["ref"]
        assert str(job_named(dag, "primary-fasta").shell) == (
            f"samtools faidx {ref} -r {t}/primary.contigs -o {t}/primary.fa"
        )
        assert str(job_named(dag, "primary-index").shell) == (
            f"sentieon minimap2 -x asm5 -t 2 -d {t}/primary.asm5.mmi "
            f"{t}/primary.fa"
        )
        assert str(job_named(dag, "haps-bam").shell) == (
            f"sentieon driver --passthru --algo Lua -e hapslice "
            f"{b}/ultima.model {t}/sample-hap.fa - | sentieon minimap2 -a "
            f"-x asm5 -O17,81 -E3,1 --secondary=no -t 2 "
            f"-R '@RG\\tID:sample1\\tSM:sample1' {t}/primary.asm5.mmi - | "
            f"sentieon util sort -t 2 -i - --sam2bam -o {t}/haps.bam"
        )
        # The table covers chr1-22, X and Y, never the calling BED
        assert str(job_named(dag, "haplen").shell) == (
            f"sentieon driver --input {t}/haps.bam --reference {ref} "
            f"--thread_count 2 --interval {t}/haplen.bed --algo Lua "
            f"-e haplen {b}/ultima.model {t}/haplen.tsv.gz"
        )

    def test_post_lua_reads_the_hap_table(self, files):
        dag = pangenome_pipeline(files).build_dag()

        b, t = files["bundle"], files["tmp"]
        assert str(job_named(dag, "post-lua").shell) == (
            f"sentieon driver --reference {files['ref']} --thread_count 2 "
            f"--interval {files['bed']} --algo Lua {b}/ultima.model "
            f"--repeat-model {t}/sample-repeat.model "
            f"--hap-table {t}/haplen.tsv.gz --lambda 2 --nofilt "
            f"{t}/sample-dnascope_transfer.vcf.gz "
            f"{t}/sample-dnascope_nofilt.vcf.gz"
        )
        apply_cmd = str(job_named(dag, "model-apply").shell)
        assert (
            f"--vcf {t}/sample-dnascope_nofilt.vcf.gz "
            f"{t}/sample-snv_apply.vcf.gz"
        ) in apply_cmd

    def test_dependencies(self, files):
        dag = pangenome_pipeline(files).build_dag()

        def deps(name):
            return dep_names(dag, job_named(dag, name))

        # The index does not wait on the lift
        assert job_named(dag, "primary-fasta") in dag.ready_jobs
        assert deps("primary-index") == ["primary-fasta"]
        assert deps("haps-bam") == ["primary-index", "vg-paths-fasta"]
        assert deps("haplen") == ["haps-bam"]
        assert deps("post-lua") == ["haplen", "merge-trim-concat"]
        assert deps("model-apply") == ["post-lua"]
        assert deps("sad-lad-update") == ["model-apply"]

    def test_thread_reservations(self, files):
        dag = pangenome_pipeline(files, cores=8).build_dag()

        threads = {
            name: job_named(dag, name).threads
            for name in (
                "primary-fasta",
                "primary-index",
                "haps-bam",
                "haplen",
                "post-lua",
            )
        }
        assert threads == {
            "primary-fasta": 0,
            "primary-index": 0,
            "haps-bam": 8,
            "haplen": 8,
            "post-lua": 8,
        }

    def test_fastq_read_group_names_the_haps(self, files):
        r1, r2 = files["tmp"] / "r1.fq.gz", files["tmp"] / "r2.fq.gz"
        pipeline = pangenome_pipeline(
            files,
            sample_input=[],
            r1_fastq=[r1],
            r2_fastq=[r2],
            fastq_readgroup={"ID": "rg2", "SM": "sample2"},
            skip_metrics=True,
            skip_multiqc=True,
        )
        dag = pipeline.build_dag()

        cmd = str(job_named(dag, "haps-bam").shell)
        assert "'@RG\\tID:sample2\\tSM:sample2'" in cmd

    def test_skip_model_apply(self, files):
        dag = pangenome_pipeline(files, skip_model_apply=True).build_dag()

        t = files["tmp"]
        names = [job.name for job in all_jobs(dag)]
        assert "model-apply" not in names
        assert f"{t}/sample-dnascope_transfer.vcf.gz" in str(
            job_named(dag, "merge-trim-concat").shell
        )
        assert str(job_named(dag, "post-lua").shell).endswith(
            f"{t}/sample-snv_apply.vcf.gz"
        )
        update = job_named(dag, "sad-lad-update")
        assert f"--input_vcf {t}/sample-snv_apply.vcf.gz" in str(update.shell)
        assert dep_names(dag, update) == ["post-lua"]

    @pytest.mark.parametrize("tech", ["ULTIMA", "ILLUMINA"])
    def test_other_platforms_are_unchanged(self, files, tech):
        dag = pangenome_pipeline(files, tech=tech).build_dag()

        names = [job.name for job in all_jobs(dag)]
        for name in (
            "post-lua",
            "primary-fasta",
            "primary-index",
            "haps-bam",
            "haplen",
        ):
            assert name not in names
        cmd = str(job_named(dag, "dnascope").shell)
        assert cmd.startswith("sentieon driver ")
        assert "RepeatModel" not in cmd
        assert dep_names(dag, job_named(dag, "model-apply")) == [
            "merge-trim-concat"
        ]

    def test_skip_small_variants_has_no_flow(self, files):
        pipeline = pangenome_pipeline(
            files, skip_small_variants=True, call_svs=True
        )
        dag = pipeline.build_dag()

        names = [job.name for job in all_jobs(dag)]
        assert "haps-bam" not in names
        cmd = str(job_named(dag, "dnascope").shell)
        assert cmd.startswith("sentieon driver ")
        assert "--read_filter" not in cmd
        assert "RepeatModel" not in cmd


class TestPangenomeValidation:
    """The bundle check runs on the bundle on disk"""

    def test_the_flow_bundle_is_accepted(self, files):
        pipeline = pangenome_pipeline(files)
        pipeline.validate_bundle()
        pipeline.logger.error.assert_not_called()
        assert pipeline.tech == "ULTIMAPE"

    @pytest.mark.parametrize("member", ["read_tables.tsv", "ultima.model"])
    def test_a_missing_member_is_refused(self, files, member):
        pipeline = pangenome_pipeline(files)
        write_ultima_pe_bundle(
            files["bundle"],
            pipeline="Sentieon pangenome",
            members=[m for m in PANGENOME_MEMBERS if m != member],
        )
        exits_with(pipeline.validate_bundle, pipeline.logger, member)

    def test_a_linear_bundle_is_refused(self, files):
        pipeline = pangenome_pipeline(files)
        write_ultima_pe_bundle(
            files["bundle"], pipeline="DNAscope", members=PANGENOME_MEMBERS
        )
        exits_with(
            pipeline.validate_bundle, pipeline.logger, "different pipeline"
        )

    def test_gvcf_is_refused(self, files):
        pipeline = pangenome_pipeline(files, gvcf=True)
        exits_with(pipeline.validate_bundle, pipeline.logger, "`--gvcf`")

    def test_without_a_bed(self, files):
        pipeline = pangenome_pipeline(files, bed=None)
        pipeline.validate_bundle()
        assert any(
            "the repeat model also counts the alt and decoy" in line
            for line in logged(pipeline.logger.warning)
        )

    def test_skip_small_variants_needs_no_flow_members(self, files):
        pipeline = pangenome_pipeline(files, skip_small_variants=True)
        write_ultima_pe_bundle(
            files["bundle"],
            pipeline="Sentieon pangenome",
            members=["extract.model", "minimap2.model"],
        )
        pipeline.validate_bundle()
        pipeline.logger.error.assert_not_called()


# --- The building blocks -------------------------------------------------


def make_ctx(tmp_path: pathlib.Path, dry_run: bool = True) -> StageContext:
    return StageContext(
        reference=tmp_path / "ref.fa",
        output_vcf=tmp_path / "output.vcf.gz",
        tmp_dir=tmp_path,
        cores=4,
        dry_run=dry_run,
        skip_version_check=True,
    )


class TestLuaScript:
    """`--algo Lua [-e entry] <script> <args...>`"""

    def test_script_and_arguments(self, tmp_path):
        algo = LuaScript(tmp_path / "s.lua", ["--opt", tmp_path / "in", "-"])
        assert algo.build_cmd() == [
            "--algo",
            "Lua",
            f"{tmp_path}/s.lua",
            "--opt",
            f"{tmp_path}/in",
            "-",
        ]

    def test_entry(self, tmp_path):
        algo = LuaScript(tmp_path / "s.lua", ["out"], entry="haplen")
        assert algo.build_cmd()[:4] == ["--algo", "Lua", "-e", "haplen"]

    def test_passthru(self, tmp_path):
        driver = Driver(passthru=True)
        driver.add_algo(LuaScript(tmp_path / "s.lua", entry="hapslice"))
        assert driver.build_cmd() == [
            "sentieon",
            "driver",
            "--passthru",
            "--algo",
            "Lua",
            "-e",
            "hapslice",
            f"{tmp_path}/s.lua",
        ]

    def test_no_passthru_by_default(self, tmp_path):
        assert "--passthru" not in Driver(reference=tmp_path).build_cmd()


class TestDriverJobEnv:
    """A per-job environment is an `env` prefix"""

    def test_env_prefix(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [LuaScript(tmp_path / "s.lua")],
            name="lua",
            task_name="lua",
            env={"A": "1", "B": "x y"},
        )
        assert str(job.shell) == (
            f"env A=1 'B=x y' sentieon driver --reference {tmp_path}/ref.fa "
            f"--thread_count 4 --algo Lua {tmp_path}/s.lua"
        )

    def test_no_env(self, tmp_path):
        job = driver_job(
            make_ctx(tmp_path),
            [LuaScript(tmp_path / "s.lua")],
            name="lua",
            task_name="lua",
        )
        assert str(job.shell).startswith("sentieon driver ")


class TestTransferApplyPost:
    """post.lua between the transfer and DNAModelApply"""

    @staticmethod
    def post_spec(tmp_path, **kwargs) -> PostSpec:
        return PostSpec(
            script=tmp_path / "ultima.model",
            repeat_model=tmp_path / "repeat.model",
            out_vcf=tmp_path / "nofilt.vcf.gz",
            **kwargs,
        )

    @staticmethod
    def upstream(dag, tmp_path, name="upstream"):
        job = rm_job([tmp_path / name], name)
        dag.add_job(job)
        return job

    def test_post_only_reads_the_raw_vcf(self, tmp_path):
        dag = DAG()
        upstream = self.upstream(dag, tmp_path)
        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            post=self.post_spec(tmp_path),
        ).add_to(dag, [upstream])

        assert result.post_job is not None
        assert str(result.post_job.shell) == (
            f"sentieon driver --reference {tmp_path}/ref.fa --thread_count 4 "
            f"--algo Lua {tmp_path}/ultima.model --repeat-model "
            f"{tmp_path}/repeat.model --lambda 2 --nofilt "
            f"{tmp_path}/raw.vcf.gz {tmp_path}/nofilt.vcf.gz"
        )
        assert dag.waiting_jobs[result.post_job] == {upstream}
        assert result.terminal == {result.post_job}
        assert result.output == tmp_path / "nofilt.vcf.gz"

    def test_transfer_post_and_apply(self, tmp_path):
        dag = DAG()
        upstream = self.upstream(dag, tmp_path)
        table_job = self.upstream(dag, tmp_path, "table")
        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            transfer=TransferSpec(
                config=TransferConfig(
                    pop_vcf=tmp_path / "pop.vcf.gz",
                    shards=[Shard("chr1", 1, 1000)],
                    pop_vcf_contigs={"chr1": 1000},
                    fai_data={"chr1": {"length": 1000}},
                ),
                out_vcf=tmp_path / "transfer.vcf.gz",
            ),
            post=self.post_spec(
                tmp_path,
                hap_table=tmp_path / "haplen.tsv.gz",
                interval=tmp_path / "calls.bed",
                interval_padding=10,
                after=frozenset({table_job}),
            ),
            apply=ApplySpec(
                model=tmp_path / "dnascope.model",
                output=tmp_path / "applied.vcf.gz",
            ),
        ).add_to(dag, [upstream])

        assert result.transfer is not None
        post, apply_job = result.post_job, result.apply_job
        assert post is not None and apply_job is not None
        assert (
            f"--interval {tmp_path}/calls.bed --interval_padding 10 "
            f"--algo Lua {tmp_path}/ultima.model --repeat-model "
            f"{tmp_path}/repeat.model --hap-table {tmp_path}/haplen.tsv.gz "
            f"--lambda 2 --nofilt {tmp_path}/transfer.vcf.gz "
            f"{tmp_path}/nofilt.vcf.gz"
        ) in str(post.shell)
        assert dag.waiting_jobs[post] == {
            result.transfer.concat_job,
            table_job,
        }
        assert f"--vcf {tmp_path}/nofilt.vcf.gz" in str(apply_job.shell)
        assert dag.waiting_jobs[apply_job] == {post}
        assert result.jobs == [*result.transfer.jobs, post, apply_job]
        assert result.terminal == {apply_job}


class TestHapTableStages:
    """The primary-contig index and the hap-length table"""

    def test_primary_contigs(self):
        assert primary_contigs({c: {"length": n} for c, n in FAI.items()}) == [
            "chr1",
            "chr2",
            "chrX",
            "chrY",
            "chrM",
            "chrUn_KI270302v1",
        ]

    @pytest.mark.parametrize("prefix", ["chr", ""])
    def test_haplen_intervals(self, prefix):
        fai = {
            f"{prefix}{c}": {"length": n}
            for c, n in [("1", 10), ("22", 20), ("X", 30), ("Y", 40)]
        }
        fai.update(
            {
                f"{prefix}23": {"length": 1},
                f"{prefix}M": {"length": 1},
                f"{prefix}1_KI270706v1_random": {"length": 1},
            }
        )
        assert haplen_intervals(fai) == [
            f"{prefix}1\t0\t10\n",
            f"{prefix}22\t0\t20\n",
            f"{prefix}X\t0\t30\n",
            f"{prefix}Y\t0\t40\n",
        ]

    def test_files_are_written_outside_a_dry_run(self, tmp_path):
        fai_data = {c: {"length": n} for c, n in FAI.items()}
        dag = DAG()
        ctx = make_ctx(tmp_path, dry_run=False)
        index = PrimaryIndexStage(
            ctx=ctx, fai_data=fai_data, prefix=tmp_path / "primary"
        ).add_to(dag)
        HapTableStage(
            ctx=ctx,
            script=tmp_path / "ultima.model",
            haplotypes=tmp_path / "hap.fa",
            index=index.index,
            sample="s",
            fai_data=fai_data,
            haps_bam=tmp_path / "haps.bam",
            table=tmp_path / "haplen.tsv.gz",
        ).add_to(dag, index.terminal)

        assert (tmp_path / "primary.contigs").read_text() == (
            "chr1\nchr2\nchrX\nchrY\nchrM\nchrUn_KI270302v1\n"
        )
        assert (tmp_path / "haplen.bed").read_text() == (
            f"chr1\t0\t{FAI['chr1']}\nchr2\t0\t1000\nchrX\t0\t2000\n"
            "chrY\t0\t3000\n"
        )

    def test_nothing_is_written_in_a_dry_run(self, tmp_path):
        PrimaryIndexStage(
            ctx=make_ctx(tmp_path),
            fai_data={"chr1": {"length": 1}},
            prefix=tmp_path / "primary",
        ).add_to(DAG())
        assert not (tmp_path / "primary.contigs").exists()
