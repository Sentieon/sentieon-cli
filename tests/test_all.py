import sys
import os
import shutil

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

import sentieon_driver  # NOQA


def test_one():
    assert sentieon_driver


def test_tmp():
    """basic test for tmp()"""
    tmp_dir, tmp_base = sentieon_driver.tmp()
    assert os.path.exists(tmp_dir)
    assert os.path.isdir(tmp_dir)
    shutil.rmtree(tmp_dir)
