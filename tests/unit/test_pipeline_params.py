"""
Every pipeline argument must be an attribute of a freshly built pipeline.

`BasePipeline.handle_arguments` raises when a `params` key is missing from
the instance, and argparse leaves an unset argument as `None`, which it
skips -- so the `__init__` value is the effective default. These tests keep
the two in step.
"""

from typing import Type

import pytest

from sentieon_cli.dnascope import DNAscopePipeline
from sentieon_cli.dnascope_hybrid import DNAscopeHybridPipeline
from sentieon_cli.dnascope_longread import DNAscopeLRPipeline
from sentieon_cli.hybrid_pangenome import HybridPangenome
from sentieon_cli.pipeline import BasePipeline
from sentieon_cli.sentieon_pangenome import SentieonPangenome

PIPELINES = [
    DNAscopePipeline,
    DNAscopeLRPipeline,
    DNAscopeHybridPipeline,
    SentieonPangenome,
    HybridPangenome,
]


@pytest.mark.parametrize("cls", PIPELINES, ids=lambda c: c.__name__)
def test_every_argument_is_an_attribute(cls: Type[BasePipeline]):
    """`handle_arguments` needs an attribute per param and positional"""
    attributes = set(vars(cls()))
    arguments = set(cls.params) | set(cls.positionals)
    assert arguments <= attributes, sorted(arguments - attributes)


@pytest.mark.parametrize("cls", PIPELINES, ids=lambda c: c.__name__)
def test_argparse_defaults_match_the_attributes(cls: Type[BasePipeline]):
    """An argparse default must match the value `__init__` assigns"""
    pipeline = cls()
    specs = list(cls.params.items()) + list(cls.positionals.items())
    mismatched = {
        key: (getattr(pipeline, key), spec["default"])
        for key, spec in specs
        if "default" in spec and getattr(pipeline, key) != spec["default"]
    }
    assert not mismatched
