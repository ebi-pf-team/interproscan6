#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) EMBL-EBI 2024


import logging
import json
import sys

from pathlib import Path
from typing import List, Optional

from src.plots import (
    plot_total_runtime,
    plot_process_runtime,
    plot_process_runtime_piechart,
    plot_overall_summary
)
from src.utilities import build_parser, load_data


logger = logging.getLogger(__name__)


def main(argv: Optional[List[str]] = None):
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if not args.data_files.exists():
        logger.error("Could not find file listing trace files at %s", args.data_files)
        sys.exit(1)

    if not args.outdir.is_dir():
        args.outdir.mkdir(parents=True, exist_ok=True)

    trace_file_paths = set()
    missing_files = set()
    group_order = []
    with open(args.data_files, "r") as fh:
        trace_file_paths = json.load(fh)
    for grp, paths in trace_file_paths.items():
        group_order.append(grp)
        for fp in paths:
            if not Path(fp).exists():
                missing_files.add(fp)

    if missing_files:
        missing_files = ''.join(missing_files)
        logger.error("Could not find the following trace file(s):\n%s", missing_files)
        sys.exit(1)

    all_data = load_data(trace_file_paths, args)
    all_data.to_csv("benchmarking/test.csv")

    plot_total_runtime(
        all_data,
        args.group_name,
        args.outdir,
        args.format,
        group_order
    )

    plot_process_runtime(
        all_data,
        args.group_name,
        args.outdir,
        args.format,
        group_order
    )

    plot_process_runtime_piechart(
        all_data,
        args.group_name,
        args.outdir,
        args.format
    )

    # plot overall memory usage
    plot_overall_summary(
        all_data,
        args.group_name,
        args.outdir,
        args.format,
        "overall_memory_usage",
        args.group_name if args.group_name else "Groups",
        args.group_name if args.group_name else "Groups",
        'raw_memory_MB',
        'Memory (MB)',
        'Memory Usuage',
        'groups',
        group_order,
    )

    # plot overall maximum memory usage
    plot_overall_summary(
        all_data,
        args.group_name,
        args.outdir,
        args.format,
        "overall_max_memory_usage",
        args.group_name if args.group_name else "Groups",
        args.group_name if args.group_name else "Groups",
        'raw_max_memory_MB',
        'Max. memory (MB)',
        'Maximum Memory Usage',
        'groups',
        group_order,
    )

    # plot mem per process
    plot_overall_summary(
        all_data,
        args.group_name,
        args.outdir,
        args.format,
        "memory_per_process",
        'raw_process',
        'Process',
        'raw_memory_MB',
        'Memory (MB)',
        'Memory Usage',
        'process',
        group_order,
    )

    # plot max mem per process
    plot_overall_summary(
        all_data,
        args.group_name,
        args.outdir,
        args.format,
        "max_memory_per_process",
        'raw_process',
        'Process',
        'raw_max_memory_MB',
        'Max. memory (MB)',
        'Memory Usage',
        'process',
        group_order,
    )


if __name__ == "__main__":
    main()
