import json
import re
import sys

from filter_ips6_hits import (
    parse_cath,
    filter_matches
)

"""Matches up the corresponding IPS6 JSON file with the cath_resolve 
output file, and then runs filter_ips6_hits.py for each pair
of matches output files"""


CATH_PATTERN = re.compile(r"^(\d+\.\d+\.\d+).*\._\.(\d+\.\d+\.\d+\.\d+)\.cath\.resolved\.out$")
JSON_PATTERN = re.compile(r"^hmmer_parsed_.*\._\.(\d+\.\d+\.\d+\.\d+)\.json$")


def main():
    """
    Args include:
    All hmmer.out files from HMMER_PARSER
    All output files from cath-resolve hits
    post-processing params
    """
    release = None
    files = {}  # keyed by cath superfamily
    for input_arg in sys.argv[1:]:
        _file = CATH_PATTERN.match(input_arg)
        if _file:
            cath_superfam = _file.group(2)
            if cath_superfam not in files:
                files[cath_superfam] = {}
            files[cath_superfam]["ips6.json"] = input_arg
            release = _file.group(1)
            continue
        _file = JSON_PATTERN.match(input_arg)

        if _file:
            cath_superfam = _file.group(1)
            if cath_superfam not in files:
                files[cath_superfam] = {}
            files[_file.group(1)]["cath.resolve"] = input_arg
            continue

        print(f"Did not recognise this input file {input_arg}")
        sys.exit(1)

    if not release:
        print("Could not get release from cath resolve output files. Terminating")
        sys.exit(1)

    for cath_superfam, file_info in files.items():
        cath_out = parse_cath(file_info["cath.resolve"])
        processed_ips6 = filter_matches(
            file_info["ips6.json"],
            cath_out,
            release
        )

        with open(f"{file_info['ips6.json']}.processed.json", "w") as fh:
            json.dump(processed_ips6, fh, indent=2)


if __name__ == "__main__":
    main()
