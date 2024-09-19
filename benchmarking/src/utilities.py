#!/usr/bin/env python 


import argparse
import logging
import re
import sys

from pathlib import Path
from typing import List, Optional

import pandas as pd


logger = logging.getLogger(__name__)

    
class ValidateFormats(argparse.Action):
    """Check the user has provided file formats."""
    def __call__(self, parser, args, values, option_string=None):
        parsed_values = set()
        valid_formats = ("pdf", "png", "svg")
        invalid = False
        for value in values:
            if value.lower() not in valid_formats:
                invalid = True
                raise ValueError(
                    (f"Invalid file format '{value}' provided. Accepted formats: "
                      "pdf, png, svg")
                )
            parsed_values.add(value.lower())
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, parsed_values)


def build_parser(argv: Optional[List] = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="IPS6-benchmarking",
        description="Summarise performance of IPS6 across one or multiple runs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "data_files",
        type=Path,
        help="Path to plain text file listing paths to Nextflow trace files. One line per path"
    )

    parser.add_argument(
        "--group_name",
        type=Path,
        default=None,
        help="Overall name of the groups. Used as legend name. E.g. 'Batch sizes'"
    )

    parser.add_argument(
        "--fig_size",
        nargs=2,
        type=int,
        default=None,
        help="Size of the final plot, width then height (space separated), e.g. 10 5"
    )

    parser.add_argument(
        "--format",
        nargs="+",
        action=ValidateFormats,
        choices=["png", "svg", "pfd"],
        type=str,
        default={'pdf'},
        help=(
            "Figure file formats. Space separated list. "
            "Case insensitive. Default PDF. Accepted formats: png, svg, pdf"
        )
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("ips6-benchmark-report"),
        help="Output directory"
    )

    parser.add_argument(
        "--raw",
        dest="raw",
        action="store_true",
        default=False,
        help=(
            "Trace files contain raw values (e.g. mem in bytes)\n"
            "not human readable values (e.g. mem in GB and MB)"
        )
    )

    parser.add_argument(
        "--save_data",
        dest="save_data",
        action="store_true",
        default=False,
        help="Save the internal dataframe to a CSV file."
    )

    if argv is None:
        return parser
    else:
        return parser.parse_args(argv)


def load_data(data_files: dict[str, list[str]], args: argparse.ArgumentParser):
    """Load all trace files and combine into a single dataframe.
    Standardise the values also if raw values not provided.

    :param config_data: dict of config file content
    :param args: clu args parser
    """
    group_name = args.group_name if args.group_name else "Groups"
    first = True
    for group, file_paths in data_files.items():
        run = 1
        for file_path in file_paths:
            if first:
                all_data = load_dataframe(file_path, group, group_name, run)
                first = False
            else:
                all_data = pd.concat([all_data, load_dataframe(file_path, group, group_name, run)])
            run += 1

    if not args.raw:
        all_data["raw_realtime"] = convert_realtime(all_data["realtime"])
        all_data["raw_process"] = convert_process_names(all_data["process"])
        all_data["raw_memory_MB"] = convert_memory(all_data["rss"])
        all_data["raw_max_memory"] = convert_memory(all_data["peak_rss"])
    else:
        all_data["raw_realtime"] = all_data["realtime"]
        all_data["raw_process"] = all_data["process"]
        all_data["raw_memory_MB"] = all_data["rss"]
        all_data["raw_max_memory"] = all_data["peak_rss"]

    return all_data


def load_dataframe(data_file: str, group: str, group_name: str, run: int) -> pd.DataFrame:
    df = pd.read_table(data_file)
    df[group_name] = [group] * len(df)
    df['Run'] = [run] * len(df)
    df = df[~df['status'].isin(['FAILED', 'ABORTED'])]
    return df


def convert_memory(mem: list[str]) -> list[str]:
    """Convert human readable value into float"""
    raw_mem = []
    for item in mem:
        if item == '-' or item == '0':
            raw_mem.append(0)
        elif item.endswith('KB'):
            raw_mem.append(float(item.split()[0]) / 1000)
        elif item.endswith('MB'):
            raw_mem.append(float(item.split()[0]))
        elif item.endswith('GB'):
            raw_mem.append(float(item.split()[0]) * 1000)
        else:
            print('Not recognised', item)

    return raw_mem


def convert_process_names(processes: list[str]) -> list[str]:
    """Convert the process names to be shorter, and a standardised
    format with MEMBER_DB: RUNER/PARSER/POST_PROCESS/FILTER
    respectively."""
    shortened_processses = []
    unrecognised_processes = False
    for process in processes:
        if process.startswith(
            (
                'PARSE_SEQUENCE', 'SEQUENCE_PRECALC',
                'AGGREGATE_PARSED_SEQS', 'AGGREGATE_RESULTS',
                'REPRESENTATIVE_DOMAINS', 'XREFS', 'WRITE_RESULTS'
            )
        ):
            shortened_processses.append(process)
        elif process.startswith('SEQUENCE_ANALYSIS:'):
            shortened_processses.append(process.replace('SEQUENCE_ANALYSIS:', ''))
        else:
            print('Unrecognised process', process)
            unrecognised_processes = True

    if unrecognised_processes:
        sys.exit(1)

    return shortened_processses


def convert_realtime(realtimes) -> list[int]:
    """Convert human written values to standardised intergers
    (in seconds)."""
    sec_pattern = re.compile(r"^(\d+\.\d+|\d+)s$")
    min_pattern = re.compile(r"^(\d+)m\s(\d+)s$")
    min_only_pattern = re.compile(r"^(\d+)m$")
    hr_sec_pattern = re.compile(r"^(\d+)h\s(\d+)s$")
    hr_min_pattern = re.compile(r"^(\d+)h\s(\d+)m$")
    hr_min_sec_pattern = re.compile(r"^(\d+)h\s(\d+)m\s(\d+)s$")
    raw_realtime = []

    for value in realtimes:
        value = value.strip()

        if value == '-':
            raw_realtime.append(0)
            continue

        if value.endswith('ms'):
            raw_realtime.append(int(value.strip('ms'))/1000)
            continue

        min_data = min_only_pattern.match(value)
        if min_data:
            raw_realtime.append(int(min_data.group(1)) * 60)
            continue

        min_data = min_pattern.match(value)
        if min_data:
            mins = int(min_data.group(1)) * 60
            secs = int(min_data.group(2))
            raw_realtime.append(mins + secs)
            continue

        sec_data = sec_pattern.match(value)
        if sec_data:
            raw_realtime.append(float(sec_data.group(1)))
            continue

        hr_sec_data = hr_sec_pattern.match(value)
        if hr_sec_data:
            raw_realtime.append((int(hr_sec_data.group(1))) * 60 * 60 + int(hr_sec_data.group(2)))
            continue

        hr_min_data = hr_min_pattern.match(value)
        if hr_min_data:
            raw_realtime.append((int(hr_min_data.group(1))) * 60 * 60 + int(hr_min_data.group(2)) * 60)
            continue

        hr_min_sec_data = hr_min_sec_pattern.match(value)
        if hr_min_sec_data:
            raw_realtime.append(
                (
                    int(hr_min_sec_data.group(1))
                ) * 60 * 60 + int(
                    hr_min_sec_data.group(2)
                ) * 60 + int(
                    hr_min_sec_data.group(3)
                )
            )
            continue

        logger.error(f"Time value not recognised: {value}")

    return raw_realtime
