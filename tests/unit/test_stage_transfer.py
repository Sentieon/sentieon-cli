"""
Unit tests for the annotation-transfer stage
"""

import asyncio
import pathlib
import shutil
import subprocess
from importlib.resources import files
from types import SimpleNamespace

import pytest

from sentieon_cli.command_strings import (
    KEEP_CONTIG_AWK,
    bcftools_concat,
    cmd_bcftools_merge_trim,
    cmd_bcftools_subset_remove,
)
from sentieon_cli.dag import DAG
from sentieon_cli.shard import Shard
from sentieon_cli.shell_pipeline import Context, InputProcSub
from sentieon_cli.stages.base import StageContext, rm_job
from sentieon_cli.stages.transfer import TransferConfig, TransferStage


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


def make_config(tmp_path: pathlib.Path) -> TransferConfig:
    """Two shards on a contig the pop VCF has, one on a contig it lacks"""
    return TransferConfig(
        pop_vcf=tmp_path / "pop.vcf.gz",
        shards=[
            Shard("chr1", 0, 1000),
            Shard("chr1", 1000, 2000),
            Shard("chrUn_extra", 0, 500),
        ],
        pop_vcf_contigs={"chr1": 2000},
        fai_data={
            "chr1": {"length": 2000},
            "chrUn_extra": {"length": 500},
        },
    )


def make_stage(tmp_path: pathlib.Path, **kwargs) -> TransferStage:
    """A TransferStage over `make_config`"""
    defaults = dict(
        ctx=make_ctx(tmp_path),
        config=make_config(tmp_path),
        raw_vcf=tmp_path / "raw.vcf.gz",
        out_vcf=tmp_path / "transfer.vcf.gz",
    )
    defaults.update(kwargs)
    return TransferStage(**defaults)  # type: ignore[arg-type]


class TestJobs:
    """The jobs the fan-out produces"""

    def test_one_job_per_shard_plus_the_extra_contig(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert [job.name for job in result.shard_jobs] == [
            "merge-trim-0",
            "merge-trim-1",
            "merge-trim-extra",
        ]
        assert result.concat_job.name == "merge-trim-concat"

    def test_task_name_and_threads(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.jobs:
            assert job.task_name == "annotation-transfer"
        for job in result.shard_jobs:
            assert job.threads == 1
        # The concat runs with every core
        assert result.concat_job.threads == 4

    def test_shard_jobs_read_the_raw_vcf(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        for job in result.shard_jobs:
            assert str(tmp_path / "raw.vcf.gz") in str(job.shell)

    def test_concat_writes_the_output(self, tmp_path):
        out_vcf = tmp_path / "annotated.vcf.gz"
        result = make_stage(tmp_path, out_vcf=out_vcf).add_to(DAG())

        assert str(out_vcf) in str(result.concat_job.shell)
        assert result.out_vcf == out_vcf


class TestMergeCommand:
    """The shard merge keeps the input VCF's QUAL and drops FILTER"""

    @staticmethod
    def shard_job(tmp_path, name):
        result = make_stage(tmp_path).add_to(DAG())
        (job,) = [job for job in result.shard_jobs if job.name == name]
        return job

    def merge_cmd(self, tmp_path):
        return self.shard_job(tmp_path, "merge-trim-0").shell.nodes[0]

    def test_merge_reads_both_inputs_unindexed(self, tmp_path):
        merge = self.merge_cmd(tmp_path)

        assert merge.args[0] == "merge"
        assert "--no-index" in merge.args
        # `-R` cannot be combined with `--no-index`; the inputs are
        # subset as they are read instead
        assert "--regions-file" not in merge.args
        assert "--regions-overlap" not in merge.args
        raw_sub, pop_sub = merge.args[-2:]
        assert isinstance(raw_sub, InputProcSub)
        assert isinstance(pop_sub, InputProcSub)

    def test_both_inputs_are_subset_to_the_shard(self, tmp_path):
        merge = self.merge_cmd(tmp_path)

        beds = set()
        for sub, vcf in zip(merge.args[-2:], ("raw.vcf.gz", "pop.vcf.gz")):
            annotate, _ = sub.node.nodes
            assert annotate.args[0] == "annotate"
            assert annotate.args[-1] == str(tmp_path / vcf)
            i = annotate.args.index("--regions-file")
            beds.add(annotate.args[i + 1])
            i = annotate.args.index("--regions-overlap")
            assert annotate.args[i + 1] == "pos"
        # The same shard BED for both inputs
        (bed,) = beds
        assert bed.endswith("sample-dnascope_transfer-subset0.bed")

    def test_removed_annotations(self, tmp_path):
        merge = self.merge_cmd(tmp_path)
        removed = []
        for sub in merge.args[-2:]:
            annotate, _ = sub.node.nodes
            removed.append(annotate.args[annotate.args.index("-x") + 1])

        # QUAL goes from the pop VCF only; FILTER from both
        assert removed == ["FILTER", "QUAL,FILTER"]

    def test_both_inputs_keep_only_the_shard_contig(self, tmp_path):
        merge = self.merge_cmd(tmp_path)
        for sub in merge.args[-2:]:
            _, keep = sub.node.nodes
            assert keep.executable == "awk"
            assert keep.args == ["-v", "contig=chr1", KEEP_CONTIG_AWK]

    def test_extra_contig_drops_filter_too(self, tmp_path):
        (annotate,) = self.shard_job(tmp_path, "merge-trim-extra").shell.nodes

        assert annotate.args[0] == "annotate"
        assert annotate.args[annotate.args.index("-x") + 1] == "FILTER"
        assert annotate.args[-1] == str(tmp_path / "raw.vcf.gz")
        i = annotate.args.index("--regions-file")
        assert annotate.args[i + 1].endswith(
            "sample-dnascope_transfer-subset2.bed"
        )


HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    '##FILTER=<ID=AC0,Description="Allele count is zero">\n'
    '##FILTER=<ID=LowQual,Description="Low quality">\n'
)
RAW_VCF = (
    HEADER
    + '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
    + '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    + "chr1\t100\t.\tC\tG\t30\tPASS\tAF=0.5\tGT\t0/1\n"
    + "chr1\t200\t.\tC\tG\t150\tLowQual\tAF=0.5\tGT\t0/1\n"
    + "chr1\t300\t.\tA\tC\t.\tPASS\tAF=0.5\tGT\t0/1\n"
    + "chr1\t500\t.\tA\tC\t25\tPASS\tAF=0.5\tGT\t0/1\n"
    + "chr1\t600\t.\tA\tC\t35\tPASS\tAF=0.5\tGT\t0/1\n"
)
POP_VCF = (
    HEADER
    + '##INFO=<ID=AF_v20,Number=A,Type=Float,Description="pop AF">\n'
    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    + "chr1\t100\trs100\tC\tG\t100\tAC0\tAF_v20=0.1\n"
    + "chr1\t200\trs200\tC\tG\t50\tPASS\tAF_v20=0.2\n"
    + "chr1\t300\trs300\tA\tC\t60\tPASS\tAF_v20=0.3\n"
    + "chr1\t550\trs550\tG\tC\t999\tPASS\tAF_v20=0.05\n"
    + "chr1\t600\trs600\tA\tG\t500\tAC0\tAF_v20=0.6\n"
)


@pytest.mark.skipif(
    shutil.which("bcftools") is None, reason="bcftools is not installed"
)
class TestMergeTrimRun:
    """The transfer's shard commands run through the shell pipeline"""

    @staticmethod
    def bgzip(tmp_path, name, text):
        plain = tmp_path / f"{name}.vcf"
        plain.write_text(text)
        out = tmp_path / f"{name}.vcf.gz"
        subprocess.run(
            ["bcftools", "view", "-Oz", "-W=tbi", "-o", str(out), str(plain)],
            check=True,
        )
        return out

    @staticmethod
    async def run(pipeline):
        context = Context()
        try:
            proc = await pipeline.run(context)
            return await proc.async_wait()
        finally:
            await context.cleanup()

    def test_keeps_the_input_qual_and_drops_filter(self, tmp_path):
        raw_vcf = self.bgzip(tmp_path, "raw", RAW_VCF)
        pop_vcf = self.bgzip(tmp_path, "pop", POP_VCF)
        subset_bed = tmp_path / "shard.bed"
        subset_bed.write_text("chr1\t0\t1000\n")
        shard_vcf = tmp_path / "shard.vcf.gz"
        trim_script = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("trimalt.py"))
        )

        pipeline = cmd_bcftools_merge_trim(
            shard_vcf,
            raw_vcf,
            pop_vcf,
            trim_script,
            subset_bed,
            "chr1",
            merge_rules="AF_v20:sum",
            merge_xargs=["--no-version", "-m", "all"],
            view_xargs=["--no-version"],
        )
        assert asyncio.run(self.run(pipeline)) == 0

        query = subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%POS\t%ID\t%ALT\t%QUAL\t%FILTER\t%INFO/AF_v20\n",
                str(shard_vcf),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        rows = [line.split("\t") for line in query.stdout.splitlines()]
        assert rows == [
            # QUAL is the input's, whichever is higher, and FILTER is
            # dropped from both inputs
            ["100", "rs100", "G", "30", ".", "0.1"],
            ["200", "rs200", "G", "150", ".", "0.2"],
            # A missing input QUAL is not filled from the pop VCF
            ["300", "rs300", "C", ".", ".", "0.3"],
            ["500", ".", "C", "25", ".", "."],
            # The pop allele is trimmed and its QUAL does not survive
            ["600", "rs600", "C", "35", ".", "0"],
        ]

    def merge_shard(self, tmp_path, raw_vcf, pop_vcf, contig, i):
        """Run the merge of one whole-contig shard"""
        subset_bed = tmp_path / f"shard{i}.bed"
        subset_bed.write_text(f"{contig}\t0\t1000\n")
        shard_vcf = tmp_path / f"shard{i}.vcf.gz"
        trim_script = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("trimalt.py"))
        )
        pipeline = cmd_bcftools_merge_trim(
            shard_vcf,
            raw_vcf,
            pop_vcf,
            trim_script,
            subset_bed,
            contig,
            merge_rules="AF_v20:sum",
            merge_xargs=["--no-version", "-m", "all"],
            view_xargs=["--no-version"],
        )
        assert asyncio.run(self.run(pipeline)) == 0
        return shard_vcf

    @staticmethod
    def contig_lines(vcf):
        header = subprocess.run(
            ["bcftools", "view", "-h", str(vcf)],
            check=True,
            capture_output=True,
            text=True,
        )
        return [
            line
            for line in header.stdout.splitlines()
            if line.startswith("##contig")
        ]

    def test_contig_lists_that_differ(self, tmp_path):
        # The input VCF lists the reference's contigs; the pop VCF lists
        # other contigs, in another order
        raw_vcf = self.bgzip(
            tmp_path,
            "raw",
            RAW_VCF.replace(
                "##contig=<ID=chr1,length=1000>\n",
                "##contig=<ID=chrM,length=1000>\n"
                "##contig=<ID=chr1,length=1000>\n"
                "##contig=<ID=chr2,length=1000>\n",
            ),
        )
        pop_vcf = self.bgzip(
            tmp_path,
            "pop",
            POP_VCF.replace(
                "##contig=<ID=chr1,length=1000>\n",
                "##contig=<ID=chr2,length=1000,assembly=b38>\n"
                "##contig=<ID=chrUn_pop,length=1000>\n"
                "##contig=<length=1000,ID=chr1>\n",
            ),
        )
        shard_vcf = self.merge_shard(tmp_path, raw_vcf, pop_vcf, "chr1", 0)

        query = subprocess.run(
            ["bcftools", "query", "-f", "%POS\t%ID\t%QUAL\n", str(shard_vcf)],
            check=True,
            capture_output=True,
            text=True,
        )
        assert query.stdout.splitlines() == [
            "100\trs100\t30",
            "200\trs200\t150",
            "300\trs300\t.",
            "500\t.\t25",
            "600\trs600\t35",
        ]
        # The shard lists its own contig, with the input VCF's line
        assert self.contig_lines(shard_vcf) == [
            "##contig=<ID=chr1,length=1000>"
        ]

    def test_concat_lists_the_contigs_in_reference_order(self, tmp_path):
        # chrM, which the pop VCF lacks, sits between the merged contigs
        contigs = "".join(
            f"##contig=<ID={c},length=1000>\n" for c in ("chr1", "chrM", "chr2")
        )
        raw_vcf = self.bgzip(
            tmp_path,
            "raw",
            RAW_VCF.replace("##contig=<ID=chr1,length=1000>\n", contigs)
            + "chrM\t10\t.\tA\tC\t40\tPASS\tAF=0.5\tGT\t0/1\n"
            + "chr2\t10\t.\tA\tC\t40\tPASS\tAF=0.5\tGT\t0/1\n",
        )
        pop_vcf = self.bgzip(
            tmp_path,
            "pop",
            POP_VCF.replace(
                "##contig=<ID=chr1,length=1000>\n",
                "##contig=<ID=chr2,length=1000>\n"
                "##contig=<ID=chr1,length=1000>\n",
            ),
        )
        shard0 = self.merge_shard(tmp_path, raw_vcf, pop_vcf, "chr1", 0)
        extra_bed = tmp_path / "extra.bed"
        extra_bed.write_text("chrM\t0\t1000\n")
        shard1 = tmp_path / "extra.vcf.gz"
        assert (
            asyncio.run(
                self.run(
                    cmd_bcftools_subset_remove(
                        shard1, raw_vcf, regions_file=extra_bed, remove="FILTER"
                    )
                )
            )
            == 0
        )
        shard2 = self.merge_shard(tmp_path, raw_vcf, pop_vcf, "chr2", 2)

        out_vcf = tmp_path / "transfer.vcf.gz"
        concat = bcftools_concat(
            out_vcf, [shard0, shard1, shard2], xargs=["--no-version"]
        )
        assert asyncio.run(self.run(concat)) == 0
        assert self.contig_lines(out_vcf) == [
            "##contig=<ID=chr1,length=1000>",
            "##contig=<ID=chrM,length=1000>",
            "##contig=<ID=chr2,length=1000>",
        ]
        query = subprocess.run(
            ["bcftools", "query", "-f", "%CHROM:%POS\n", str(out_vcf)],
            check=True,
            capture_output=True,
            text=True,
        )
        assert query.stdout.splitlines()[-2:] == ["chrM:10", "chr2:10"]

    def test_extra_contig_subset_drops_filter(self, tmp_path):
        raw_vcf = self.bgzip(tmp_path, "raw", RAW_VCF)
        subset_bed = tmp_path / "shard.bed"
        subset_bed.write_text("chr1\t150\t1000\n")
        out_vcf = tmp_path / "subset.vcf.gz"

        pipeline = cmd_bcftools_subset_remove(
            out_vcf, raw_vcf, regions_file=subset_bed, remove="FILTER"
        )
        assert asyncio.run(self.run(pipeline)) == 0

        query = subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%POS\t%QUAL\t%FILTER\n",
                str(out_vcf),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert query.stdout.splitlines() == [
            "200\t150\t.",
            "300\t.\t.",
            "500\t25\t.",
            "600\t35\t.",
        ]
        assert (tmp_path / "subset.vcf.gz.tbi").exists()


class TestTag:
    """`tag` disambiguates two fan-outs in one run"""

    def test_tagged_names(self, tmp_path):
        result = make_stage(tmp_path, tag="hp").add_to(DAG())

        assert [job.name for job in result.shard_jobs] == [
            "merge-trim-hp-0",
            "merge-trim-hp-1",
            "merge-trim-hp-extra",
        ]
        assert result.concat_job.name == "merge-trim-hp-concat"


class TestDagWiring:
    """The edges the stage inserts"""

    def test_shards_depend_on_upstream_and_concat_on_the_shards(
        self, tmp_path
    ):
        dag = DAG()
        upstream = rm_job([tmp_path / "upstream"], "upstream")
        dag.add_job(upstream)

        result = make_stage(tmp_path).add_to(dag, [upstream])

        for job in result.shard_jobs:
            assert dag.waiting_jobs[job] == {upstream}
        assert dag.waiting_jobs[result.concat_job] == set(result.shard_jobs)

    def test_shards_are_roots_without_upstream(self, tmp_path):
        dag = DAG()
        result = make_stage(tmp_path).add_to(dag)

        for job in result.shard_jobs:
            assert job in dag.ready_jobs

    def test_result_fields(self, tmp_path):
        result = make_stage(tmp_path).add_to(DAG())

        assert result.jobs == [*result.shard_jobs, result.concat_job]
        assert result.terminal == {result.concat_job}


class TestFromPipeline:
    """`TransferConfig.from_pipeline` collecting a pipeline's attributes"""

    @staticmethod
    def make_pipeline(tmp_path: pathlib.Path, **kwargs) -> SimpleNamespace:
        """A stand-in carrying the four `TransferInputs` attributes"""
        config = make_config(tmp_path)
        defaults = dict(
            pop_vcf=config.pop_vcf,
            shards=config.shards,
            pop_vcf_contigs=config.pop_vcf_contigs,
            fai_data=config.fai_data,
        )
        defaults.update(kwargs)
        return SimpleNamespace(**defaults)

    def test_collects_the_four_values(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path)

        config = TransferConfig.from_pipeline(pipeline)

        assert config.pop_vcf == pipeline.pop_vcf
        assert config.shards == pipeline.shards
        assert config.pop_vcf_contigs == pipeline.pop_vcf_contigs
        assert config.fai_data == pipeline.fai_data

    def test_without_a_pop_vcf(self, tmp_path):
        pipeline = self.make_pipeline(tmp_path, pop_vcf=None)

        with pytest.raises(ValueError, match="population VCF"):
            TransferConfig.from_pipeline(pipeline)
