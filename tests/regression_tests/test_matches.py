#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test for mismatches in hits (presence/absence) between two versions of IPS
# Default args assuming this is running in the root of the repo
import ast
import argparse
import difflib
import json

def main():
    parser = argparse.ArgumentParser(prog="IPS_match_regression_test", description="Check presence of matches", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--expected", type=str, default="tests/data/test_prot.fa.json", help="JSON with expected results")
    parser.add_argument("--observed", type=str, default="test_prot.fa.json", help="JSON output file from IPS6")
    parser.add_argument("--summary", action="store_true", help="Print only the summary message")
    parser.add_argument("--format", choices=["json", "tsv", "xml"], default="json", help="Format of input files. 'json' [default], 'tsv', or 'xml'")
    args = parser.parse_args()

    with open(args.expected, "r") as fh:
        if args.format == "json":
            expected = parse_json(json.load(fh))
        elif args.format == "tsv":
            expected = parse_tsv(fh)

    with open(args.observed, "r") as fh:
        if args.format == "json":
            observed = parse_json(json.load(fh))
        elif args.format == "tsv":
            observed = parse_tsv(fh)

    diff = difflib.ndiff(sorted(expected), sorted(observed))
    expected_only, observed_only, both = 0, 0, 0
    for line in diff:
        if line.startswith(("-", "+", " ")):  # skip summary lines such as "?    ^^^^  -    ^^  ^ ^"
            try:
                tab_separated = '\t'.join(map(str, ast.literal_eval(line[2:])))  # Safely eval str and convert to tab-sep str
                if line.startswith("-"):
                    if not args.summary:
                        print(f"< {tab_separated}")  # Only in expected
                    expected_only += 1
                elif line.startswith("+"):
                    if not args.summary:
                        print(f"> {tab_separated}")  # Only in observed
                    observed_only += 1
                elif line.startswith(" "):
                    if not args.summary:
                        print(f"- {tab_separated}")  # In both
                    both += 1
            except (SyntaxError, ValueError):
                print(f"Could not parse line {line}")

    print(
        f"============ Summary ============\n"
        f"Matches only in expected : {expected_only}\n"
        f"Matches only in observed : {observed_only}\n"
        f"Matches in both          : {both}"
    )


def parse_json(iprsn_dict: dict):
    """Convert all matches into lists, keeping only the data we care about"""
    matches = []  # [[md5, sig_acc, loc start, loc end]]
    for protein_dict in iprsn_dict["results"]:
        md5 = protein_dict["md5"].upper()  # account for diff between iprscn 5 and 6
        for match in protein_dict["matches"]:
            sig_acc = match["signature"]["accession"]
            for loc in match["locations"]:
                matches.append([md5, sig_acc, loc["start"], loc["end"]])
    return [repr(nested) for nested in matches]

def parse_tsv(file_handle):
    matches = []
    for line in file_handle:
        data = line.split("\t")
        matches.append([data[1], data[4], data[6], data[7]])
    return [repr(nested) for nested in matches]

if __name__ == "__main__":
    main()
