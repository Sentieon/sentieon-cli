#!/bin/bash
# Regenerate the binary pangenome fixtures in this directory.
#
# The `.gfa` files are the checked-in sources; `vg` (tested with v1.72.0)
# turns them into the GBZ graphs and the `vg haplotypes` index that
# `tests/unit/test_pangenome_meta.py` parses. Everything here is a toy
# graph of 8 segments and two 30 bp contigs, so the outputs stay a few KB.
#
# The graphs differ only in their walks and in the `RS:Z:` header tag,
# which `vg` turns into the GBWT `reference_samples` tag:
#
#   tiny.gfa            GRCh38 is the unfragmented backbone, CHM13 is
#                       split into two fragments on chr1
#   tiny_chm13.gfa      the same graph with the two roles swapped
#   tiny_no_backbone.gfa    both reference samples are fragmented
#   tiny_no_ref_tag.gfa     no `RS:Z:` tag, so no reference sample at all
#
# `tiny_sampled.gbz` stands in for the graph `vg haplotypes
# --set-reference GRCh38` writes: one reference sample, as
# `check_sample_gbz` requires.
#
# The `rgfa_*.gfa` files are hand-written heads of a `vg convert -f -Q`
# rGFA and need no regeneration.
set -euo pipefail
cd "$(dirname "$0")"

for name in tiny tiny_chm13 tiny_no_backbone tiny_no_ref_tag; do
    vg gbwt -G "${name}.gfa" --gbz-format -g "${name}.gbz"
done

vg gbwt -Z tiny.gbz --set-reference GRCh38 -g tiny_sampled.gbz

vg index -j tiny.dist tiny.gbz
vg gbwt -Z tiny.gbz -r tiny.ri
vg haplotypes -d tiny.dist -r tiny.ri -H tiny.hapl tiny.gbz
rm -f tiny.dist tiny.ri

# What the fixtures should contain
vg gbwt -Z tiny.gbz --tags
vg paths -x tiny.gbz -L
