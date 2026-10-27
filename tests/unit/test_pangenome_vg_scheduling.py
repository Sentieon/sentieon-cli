"""vg-haplotypes runs alongside bwa-extract when both fit in memory"""

import os
import sys
from unittest.mock import patch

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
)

from test_sentieon_pangenome import TestSentieonPangenome  # noqa: E402

GB = 1024**3
TWO_NODES = ["0-95", "96-191"]


class TestVgHaplotypesScheduling(TestSentieonPangenome):
    def build(self, memory_gb, cores=64, nodes=None):
        pipeline = self.create_fastq_pipeline()
        pipeline.cores = cores
        with patch(
            "sentieon_cli.util.find_numa_nodes", return_value=nodes or []
        ):
            pipeline.configure()
        with patch(
            "sentieon_cli.sentieon_pangenome.total_memory",
            return_value=memory_gb * GB,
        ):
            dag = pipeline.build_dag()
        everything = list(dag.waiting_jobs) + list(dag.ready_jobs)
        vg = next(j for j in everything if j.name == "vg-haplotypes")
        return pipeline, dag, vg

    def test_runs_alongside_bwa_with_enough_memory(self):
        """No bwa dependency, and no threads reserved from the scheduler"""
        _, dag, vg = self.build(384, cores=192, nodes=TWO_NODES)
        deps = {d.name for d in self._get_deps(dag, vg)}
        assert "bwa-extract" not in deps
        assert "kmc" in deps
        assert vg.threads == 0
        # vg itself still uses every core
        assert "-t 192" in str(vg.shell)

    def test_waits_for_bwa_when_memory_is_short(self):
        """A 128 GB host has no room for vg next to a bwa lane"""
        pipeline, dag, vg = self.build(128)
        deps = {d.name for d in self._get_deps(dag, vg)}
        assert "bwa-extract" in deps
        assert vg.threads == pipeline.cores

    def test_waits_for_every_bwa_lane(self):
        """The memory check counts each lane, and vg waits for all of them"""
        _, dag, vg = self.build(200, cores=192, nodes=TWO_NODES)
        bwa = [d for d in self._get_deps(dag, vg) if d.name == "bwa-extract"]
        assert len(bwa) == 4
        assert vg.threads == 192
