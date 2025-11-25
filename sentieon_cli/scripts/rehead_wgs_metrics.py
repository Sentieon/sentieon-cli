#!/usr/bin/env python3

"""
Rehead the WGSMetricsAlgo output for MultiQC
"""

# Copyright (c) 2025 Sentieon Inc. All rights reserved

import argparse
import os.path
import shutil
import tempfile


def process_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metrics_file",
        required=True,
        help="The input metrics file",
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    metrics_file = args.metrics_file
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_metrics = os.path.join(tmp_dir, os.path.basename(metrics_file))
        shutil.move(metrics_file, temp_metrics)

        with (
            open(temp_metrics, 'r') as in_fh,
            open(metrics_file, 'w') as out_fh
        ):
            print("'## METRICS CLASS WgsMetrics'", file=out_fh)
            line = in_fh.readline()
            for line in in_fh:
                print(line.rstrip(), file=out_fh)


if __name__ == "__main__":
    args = process_args()
    main(args)
