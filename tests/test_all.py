import sys
import os
import shutil

from importlib_resources import files

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

import sentieon_cli  # NOQA


def test_one():
    assert sentieon_cli


def test_tmp():
    """basic test for tmp()"""
    tmp_dir = sentieon_cli.tmp()
    assert os.path.exists(tmp_dir.name)
    assert os.path.isdir(tmp_dir.name)
    shutil.rmtree(tmp_dir.name)

def test_find_script():
    gvcf_combine = files('sentieon_cli.scripts').joinpath('gvcf_combine.py')
    vcf_mod = files('sentieon_cli.scripts').joinpath('vcf_mod.py')
    assert os.path.isfile(str(gvcf_combine))
    assert os.path.isfile(str(vcf_mod))
