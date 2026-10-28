"""bwa-extract and mm2-lift run as one lane per NUMA node"""

import os
import sys
from unittest.mock import patch

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

import test_sentieon_pangenome  # noqa: E402

TWO_NODES = ["0-95", "96-191"]
LANES = ["0-47", "48-95", "96-143", "144-191"]


class TestNumaLanes:
    # Use the pangenome tests' fixtures without inheriting, and so
    # re-running, their tests
    def setup_method(self):
        self.base = test_sentieon_pangenome.TestSentieonPangenome()
        self.base.setup_method()

    def lane_pipeline(self, cores=192, nodes=TWO_NODES, no_split=False):
        pipeline = self.base.create_fastq_pipeline()
        pipeline.cores = cores
        pipeline.no_split_alignment = no_split
        with patch("sentieon_cli.util.find_numa_nodes", return_value=nodes):
            pipeline.configure()
        return pipeline

    def jobs(self, dag, name):
        """The step's jobs: `name`, or `name-<lane>` when split"""
        everything = list(dag.waiting_jobs) + list(dag.ready_jobs)
        return [
            j
            for j in everything
            if j.name == name or j.name.rpartition("-")[0] == name
        ]

    def test_configure_splits_like_dnascope(self):
        assert self.lane_pipeline().numa_nodes == LANES

    def test_bwa_extract_lanes(self):
        dag = self.lane_pipeline().build_dag()
        bwa = self.jobs(dag, "bwa-extract")
        assert len(bwa) == 4
        for j, cpus in enumerate(LANES):
            job = next(b for b in bwa if f"taskset -c {cpus} " in str(b.shell))
            cmd = str(job.shell)
            assert job.name == f"bwa-extract-{j}"
            assert job.threads == 48
            assert job.resources == {f"node{j}": 1}
            assert f"fqidx extract -F {j}/4 -K 10000000" in cmd
            assert "bwa mem" in cmd and cmd.count("-K 10000000") == 2
            assert " -p " in cmd and "-t 48" in cmd
            assert (
                f"sample-bwa-{j}.bam" in cmd and f"sample-bwa-{j}.fq.gz" in cmd
            )
            # pgutil extract is pinned too
            assert cmd.count(f"taskset -c {cpus}") == 2

    def test_mm2_lanes_follow_their_bwa_lane(self):
        pipeline = self.lane_pipeline()
        dag = pipeline.build_dag()
        mm2 = self.jobs(dag, "mm2-lift")
        assert len(mm2) == 4
        for j, cpus in enumerate(LANES):
            job = next(
                m for m in mm2 if f"sample-bwa-{j}.fq.gz" in str(m.shell)
            )
            cmd = str(job.shell)
            assert job.name == f"mm2-lift-{j}"
            assert job.threads == 48 and job.resources == {f"node{j}": 1}
            assert cmd.count(f"taskset -c {cpus}") == 3  # minimap2, lift, sort
            assert f"sample-mm2-{j}.bam" in cmd
            deps = self.base._get_deps(dag, job)
            lane_bwa = [d for d in deps if d.name.startswith("bwa-extract")]
            assert [d.name for d in lane_bwa] == [f"bwa-extract-{j}"]
            assert f"sample-bwa-{j}.fq.gz" in str(lane_bwa[0].shell)
            assert {"check-sample-gfa", "vg-paths-fasta"} <= {
                d.name for d in deps
            }

    def test_dedups_read_every_lane(self):
        dag = self.lane_pipeline().build_dag()
        everything = list(dag.waiting_jobs) + list(dag.ready_jobs)
        lc = [
            str(j.shell)
            for j in everything
            if "LocusCollector" in str(j.shell)
        ]
        bwa_lc = next(c for c in lc if "sample-bwa-0.bam" in c)
        mm2_lc = next(c for c in lc if "sample-mm2-0.bam" in c)
        for j in range(4):
            assert f"sample-bwa-{j}.bam" in bwa_lc
            assert f"sample-mm2-{j}.bam" in mm2_lc

    def test_no_split_alignment_keeps_one_job(self):
        dag = self.lane_pipeline(no_split=True).build_dag()
        bwa = self.jobs(dag, "bwa-extract")
        mm2 = self.jobs(dag, "mm2-lift")
        assert len(bwa) == 1 and len(mm2) == 1
        assert "taskset" not in str(bwa[0].shell) + str(mm2[0].shell)
        assert "fqidx" not in str(bwa[0].shell)
        assert "sample-bwa.bam" in str(bwa[0].shell)
        assert bwa[0].name == "bwa-extract" and mm2[0].name == "mm2-lift"
        # the same -K as a split run, so the alignments match
        assert "-K 10000000" in str(bwa[0].shell)
        assert "-K 2000000" not in str(bwa[0].shell)

    def test_small_machine_keeps_one_job(self):
        dag = self.lane_pipeline(cores=48, nodes=["0-47"]).build_dag()
        assert len(self.jobs(dag, "bwa-extract")) == 1
        assert len(self.jobs(dag, "mm2-lift")) == 1

    def test_job_names_are_unique(self):
        dag = self.lane_pipeline().build_dag()
        everything = list(dag.waiting_jobs) + list(dag.ready_jobs)
        names = [
            j.name
            for j in everything
            if j.name.startswith(("bwa-extract", "mm2-lift"))
        ]
        assert len(names) == 8 and len(set(names)) == 8

    def test_bam_input_is_unchanged(self):
        pipeline = self.base.create_pipeline()
        pipeline.cores = 192
        with patch(
            "sentieon_cli.util.find_numa_nodes", return_value=TWO_NODES
        ):
            pipeline.configure()
        assert pipeline.numa_nodes == []
        dag = pipeline.build_dag()
        mm2 = self.jobs(dag, "mm2-lift")
        assert len(mm2) == 1 and "taskset" not in str(mm2[0].shell)
