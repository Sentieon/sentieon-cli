"""
Unit tests for the small-variant stages
"""

import pathlib

import pytest

from sentieon_cli.dag import DAG
from sentieon_cli.driver import DNAscope, SVSolver
from sentieon_cli.shard import Shard
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.small_variants import (
    ApplySpec,
    DNAscopeStage,
    DriverStage,
    GVCFtyperStage,
    ModelApplyStage,
    TransferApplyStage,
    TransferSpec,
)
from sentieon_cli.stages.transfer import TransferConfig


def make_ctx(tmp_path: pathlib.Path, cores: int = 4) -> StageContext:
    """A StageContext over a temporary directory"""
    return StageContext(
        reference=tmp_path / "ref.fa",
        output_vcf=tmp_path / "output.vcf.gz",
        tmp_dir=tmp_path,
        cores=cores,
        # `build_transfer_jobs` shells out to bcftools unless this is set
        dry_run=True,
        skip_version_check=True,
    )


def make_call_stage(tmp_path: pathlib.Path, **kwargs) -> DNAscopeStage:
    """A DNAscopeStage calling SNVs and indels"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        algos=[DNAscope(tmp_path / "raw.vcf.gz")],
        inputs=[tmp_path / "sample.cram"],
    )
    defaults.update(kwargs)
    return DNAscopeStage(**defaults)  # type: ignore[arg-type]


def make_transfer_config(tmp_path: pathlib.Path) -> TransferConfig:
    """One shard the pop VCF has, one it lacks"""
    return TransferConfig(
        pop_vcf=tmp_path / "pop.vcf.gz",
        shards=[Shard("chr1", 0, 1000), Shard("chrUn_extra", 0, 500)],
        pop_vcf_contigs={"chr1": 1000},
        fai_data={
            "chr1": {"length": 1000},
            "chrUn_extra": {"length": 500},
        },
    )


class TestDriverStage:
    """The generic single-job driver stage"""

    def test_minimal_command(self, tmp_path):
        result = DriverStage(
            ctx=make_ctx(tmp_path),
            algos=[
                SVSolver(
                    output=tmp_path / "sv.vcf.gz",
                    vcf=tmp_path / "sv_tmp.vcf.gz",
                )
            ],
            name="svsolver",
            task_name="sv-calling",
        ).add_to(DAG())

        assert str(result.job.shell) == (
            f"sentieon driver --reference {tmp_path}/ref.fa "
            f"--thread_count 4 --algo SVSolver "
            f"--vcf {tmp_path}/sv_tmp.vcf.gz {tmp_path}/sv.vcf.gz"
        )

    def test_name_task_name_and_default_threads(self, tmp_path):
        result = DriverStage(
            ctx=make_ctx(tmp_path),
            algos=[DNAscope(tmp_path / "raw.vcf.gz")],
            name="a-name",
            task_name="a-task",
        ).add_to(DAG())

        assert result.job.name == "a-name"
        assert result.job.task_name == "a-task"
        assert result.job.threads == 4

    def test_threads_override_leaves_the_driver_alone(self, tmp_path):
        result = DriverStage(
            ctx=make_ctx(tmp_path),
            algos=[DNAscope(tmp_path / "raw.vcf.gz")],
            name="a-name",
            task_name="a-task",
            threads=1,
        ).add_to(DAG())

        assert result.job.threads == 1
        assert "--thread_count 4" in str(result.job.shell)

    def test_result_fields(self, tmp_path):
        result = make_call_stage(tmp_path).add_to(DAG())

        assert result.jobs == [result.job]
        assert result.terminal == {result.job}


class TestDNAscopeStage:
    """The DNAscope calling pass"""

    def test_default_name_and_task_name(self, tmp_path):
        result = make_call_stage(tmp_path).add_to(DAG())

        assert result.job.name == "dnascope"
        assert result.job.task_name == "variant-calling"

    def test_inputs_and_interval(self, tmp_path):
        bed = tmp_path / "regions.bed"
        result = make_call_stage(
            tmp_path,
            inputs=[tmp_path / "one.cram", tmp_path / "two.cram"],
            interval=bed,
        ).add_to(DAG())

        shell = str(result.job.shell)
        assert f"--input {tmp_path}/one.cram" in shell
        assert f"--input {tmp_path}/two.cram" in shell
        assert f"--interval {bed}" in shell

    def test_interval_padding_is_absent_unless_set(self, tmp_path):
        result = make_call_stage(tmp_path).add_to(DAG())

        assert "--interval_padding" not in str(result.job.shell)

    def test_interval_padding(self, tmp_path):
        result = make_call_stage(tmp_path, interval_padding=100).add_to(DAG())

        assert "--interval_padding 100" in str(result.job.shell)

    def test_read_filter(self, tmp_path):
        result = make_call_stage(
            tmp_path, read_filter=["UltimaReadFilter"]
        ).add_to(DAG())

        assert "--read_filter UltimaReadFilter" in str(result.job.shell)

    def test_empty_read_filter_emits_nothing(self, tmp_path):
        result = make_call_stage(tmp_path, read_filter=[]).add_to(DAG())

        assert "--read_filter" not in str(result.job.shell)

    def test_replace_rg_precedes_its_input(self, tmp_path):
        result = make_call_stage(
            tmp_path,
            inputs=[tmp_path / "one.cram", tmp_path / "two.cram"],
            replace_rg=[["rg1=rg1-lr"], ["rg2=rg2-lr"]],
        ).add_to(DAG())

        assert str(result.job.shell).startswith(
            f"sentieon driver --replace_rg rg1=rg1-lr "
            f"--input {tmp_path}/one.cram "
            f"--replace_rg rg2=rg2-lr --input {tmp_path}/two.cram"
        )

    def test_algos_keep_their_order(self, tmp_path):
        result = make_call_stage(
            tmp_path,
            algos=[
                DNAscope(tmp_path / "raw.vcf.gz"),
                DNAscope(tmp_path / "sv_tmp.vcf.gz", var_type="BND"),
            ],
        ).add_to(DAG())

        shell = str(result.job.shell)
        assert shell.index("raw.vcf.gz") < shell.index("--var_type BND")

    def test_upstream_becomes_the_dependency(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_call_stage(tmp_path).add_to(dag, [upstream])

        assert dag.waiting_jobs[result.job] == {upstream}


class TestModelApplyStage:
    """DNAModelApply"""

    def test_command(self, tmp_path):
        result = ModelApplyStage(
            ctx=make_ctx(tmp_path),
            model=tmp_path / "bundle" / "dnascope.model",
            vcf=tmp_path / "raw.vcf.gz",
            output=tmp_path / "applied.vcf.gz",
        ).add_to(DAG())

        assert str(result.job.shell) == (
            f"sentieon driver --reference {tmp_path}/ref.fa "
            f"--thread_count 4 --algo DNAModelApply "
            f"--model {tmp_path}/bundle/dnascope.model "
            f"--vcf {tmp_path}/raw.vcf.gz "
            f"{tmp_path}/applied.vcf.gz"
        )

    def test_defaults_and_result(self, tmp_path):
        output = tmp_path / "applied.vcf.gz"
        result = ModelApplyStage(
            ctx=make_ctx(tmp_path),
            model=tmp_path / "dnascope.model",
            vcf=tmp_path / "raw.vcf.gz",
            output=output,
        ).add_to(DAG())

        assert result.job.name == "model-apply"
        assert result.job.task_name == "model-apply"
        assert result.job.threads == 4
        assert result.output == output
        assert result.terminal == {result.job}

    def test_renamed_job(self, tmp_path):
        result = ModelApplyStage(
            ctx=make_ctx(tmp_path),
            model=tmp_path / "dnascope.model",
            vcf=tmp_path / "raw.vcf.gz",
            output=tmp_path / "applied.vcf.gz",
            name="snv-apply",
        ).add_to(DAG())

        assert result.job.name == "snv-apply"


class TestGVCFtyperStage:
    """Genotyping a gVCF"""

    def test_command_without_an_interval(self, tmp_path):
        result = GVCFtyperStage(
            ctx=make_ctx(tmp_path),
            gvcf=tmp_path / "sample.g.vcf.gz",
            output=tmp_path / "sample.vcf.gz",
        ).add_to(DAG())

        assert str(result.job.shell) == (
            f"sentieon driver --reference {tmp_path}/ref.fa "
            f"--thread_count 4 --algo GVCFtyper "
            f"--vcf {tmp_path}/sample.g.vcf.gz "
            f"{tmp_path}/sample.vcf.gz"
        )
        assert result.job.name == "gvcftyper"
        assert result.job.task_name == "gvcftyper"

    def test_interval(self, tmp_path):
        bed = tmp_path / "regions.bed"
        result = GVCFtyperStage(
            ctx=make_ctx(tmp_path),
            gvcf=tmp_path / "sample.g.vcf.gz",
            output=tmp_path / "sample.vcf.gz",
            interval=bed,
        ).add_to(DAG())

        assert f"--interval {bed}" in str(result.job.shell)
        assert result.output == tmp_path / "sample.vcf.gz"


class TestTransferApplyStage:
    """Transfer, apply, or both"""

    def test_neither_half_is_an_error(self, tmp_path):
        with pytest.raises(ValueError):
            TransferApplyStage(
                ctx=make_ctx(tmp_path),
                raw_vcf=tmp_path / "raw.vcf.gz",
            ).add_to(DAG())

    def test_apply_only_reads_the_raw_vcf(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            apply=ApplySpec(
                model=tmp_path / "dnascope.model",
                output=tmp_path / "applied.vcf.gz",
            ),
        ).add_to(dag, [upstream])

        assert result.transfer is None
        assert result.apply_job is not None
        assert f"--vcf {tmp_path}/raw.vcf.gz" in str(result.apply_job.shell)
        assert dag.waiting_jobs[result.apply_job] == {upstream}
        assert result.terminal == {result.apply_job}
        assert result.output == tmp_path / "applied.vcf.gz"

    def test_transfer_only(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            transfer=TransferSpec(
                config=make_transfer_config(tmp_path),
                out_vcf=tmp_path / "transfer.vcf.gz",
            ),
        ).add_to(dag, [upstream])

        assert result.apply_job is None
        assert result.transfer is not None
        assert result.terminal == {result.transfer.concat_job}
        assert result.output == tmp_path / "transfer.vcf.gz"
        for job in result.transfer.shard_jobs:
            assert dag.waiting_jobs[job] == {upstream}

    def test_both_halves(self, tmp_path):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            transfer=TransferSpec(
                config=make_transfer_config(tmp_path),
                out_vcf=tmp_path / "transfer.vcf.gz",
            ),
            apply=ApplySpec(
                model=tmp_path / "dnascope.model",
                output=tmp_path / "applied.vcf.gz",
            ),
        ).add_to(dag, [upstream])

        assert result.transfer is not None
        assert result.apply_job is not None
        # DNAModelApply reads the transfer output
        assert f"--vcf {tmp_path}/transfer.vcf.gz" in str(
            result.apply_job.shell
        )
        # ... and waits on the transfer alone; the concat already depends
        # on the raw-VCF producer
        assert dag.waiting_jobs[result.apply_job] == {
            result.transfer.concat_job,
        }
        assert result.terminal == {result.apply_job}
        assert result.jobs == [
            *result.transfer.jobs,
            result.apply_job,
        ]

    def test_transfer_tag_reaches_the_job_names(self, tmp_path):
        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            transfer=TransferSpec(
                config=make_transfer_config(tmp_path),
                out_vcf=tmp_path / "transfer.vcf.gz",
                tag="hp",
            ),
        ).add_to(DAG())

        assert result.transfer is not None
        assert result.transfer.concat_job.name == "merge-trim-hp-concat"

    def test_apply_name_override(self, tmp_path):
        result = TransferApplyStage(
            ctx=make_ctx(tmp_path),
            raw_vcf=tmp_path / "raw.vcf.gz",
            apply=ApplySpec(
                model=tmp_path / "dnascope.model",
                output=tmp_path / "applied.vcf.gz",
                name="snv-apply",
            ),
        ).add_to(DAG())

        assert result.apply_job is not None
        assert result.apply_job.name == "snv-apply"
