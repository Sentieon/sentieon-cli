import os
import sys

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

import sentieon_cli.util  # NOQA


def test_split_numa_nodes():
    """Test the NUMA split func"""
    res = sentieon_cli.util.split_numa_nodes(["0-97"])
    assert res == ["0-48", "49-97"]
    res = sentieon_cli.util.split_numa_nodes(res)
    assert res == ["0-24", "25-48", "49-73", "74-97"]

    res = sentieon_cli.util.split_numa_nodes(["0-15,32-47", "16-31,48-63"])
    assert res == ["0-15", "32-47", "16-31", "48-63"]
    res = sentieon_cli.util.split_numa_nodes(res)
    assert res == [
        "0-7",
        "8-15",
        "32-39",
        "40-47",
        "16-23",
        "24-31",
        "48-55",
        "56-63",
    ]
