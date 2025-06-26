#!/usr/bin/env python

"""
Estimate sample ploidy from bam index statistics
"""

# Copyright (c) 2025 Sentieon Inc. All rights reserved

import argparse
import json
import os
import subprocess
import sys
from typing import Any, Dict, Tuple

DEFAULT_CONTIGS = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY",
]
DEFAULT_AUTOSOMES = DEFAULT_CONTIGS[:-2]

DEFAULT_SEXES = {
    "male": {"chrX": 1, "chrY": 1},
    "female": {"chrX": 2, "chrY": 0},
}

IDXSTATS_CMD = ["samtools", "idxstats"]


def process_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--input_bam", required=True, help="The input bam file"
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=DEFAULT_CONTIGS,
        help="Contigs to process",
    )
    parser.add_argument(
        "--autosomes",
        nargs="+",
        default=DEFAULT_AUTOSOMES,
        help="Autosome contigs",
    )
    parser.add_argument(
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("w"),
        help="The output json file",
    )
    return parser.parse_args()


def main(args: argparse.Namespace):
    contigs: Dict[str, bool] = {x: True for x in args.contigs}
    autosomes: Dict[str, bool] = {x: True for x in args.autosomes}

    idxstats_results: Dict[str, Tuple[int, int]] = {}
    cmd = IDXSTATS_CMD + [args.input_bam]
    res = subprocess.run(
        cmd, check=True, capture_output=True, universal_newlines=True
    )
    for line in res.stdout.split("\n"):
        line_split = line.split("\t")
        chrom = line_split[0]

        if chrom in contigs:
            idxstats_results[chrom] = (int(line_split[1]), int(line_split[2]))

    # Calculate the average reads per base over autosomes
    total_bases = sum(
        [x[0] for contig, x in idxstats_results.items() if contig in autosomes]
    )
    total_reads = sum(
        [x[1] for contig, x in idxstats_results.items() if contig in autosomes]
    )
    average_reads_per_base = total_reads / total_bases

    # Calculate read coverage by contig
    results: Dict[str, Any] = {"contigs": {}}
    for contig, (contig_l, contig_reads) in idxstats_results.items():
        contig_ploidy = "Unknown"
        contig_coverage = contig_reads / contig_l
        normalized_coverage = contig_coverage / average_reads_per_base * 2
        if normalized_coverage < 0.4:
            contig_ploidy = "0"
        elif normalized_coverage > 0.6 and normalized_coverage < 1.4:
            contig_ploidy = "1"
        elif normalized_coverage > 1.6 and normalized_coverage < 2.4:
            contig_ploidy = "2"
        elif normalized_coverage > 2.6 and normalized_coverage < 3.4:
            contig_ploidy = "3"
        elif normalized_coverage > 3.4:
            contig_ploidy = "3+"

        results["contigs"][contig] = {
            "normalized_coverage": str(normalized_coverage),
            "ploidy": contig_ploidy,
        }

    # Determine sex
    sex = "Unknown"
    for possible_sex, required_chroms in DEFAULT_SEXES.items():
        match = True
        for contig, req_ploidy in required_chroms.items():
            try:
                found_ploidy = int(results["contigs"][contig]["ploidy"])
                if found_ploidy != req_ploidy:
                    match = False
            except ValueError:
                match = False

        if match:
            sex = possible_sex

    results["sex"] = sex

    json.dump(results, args.outfile, indent=2)
    return sys.exit(os.EX_OK)


if __name__ == "__main__":
    args = process_args()
    main(args)
