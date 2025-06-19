#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test for mismatches in hits (presence/absence) between two versions of IPS
# Default args assuming this is running in the root of the repo
import ast
import argparse
import difflib
import json
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(prog="IPS_match_regression_test", description="Check presence of matches", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--expected", type=str, default="tests/data/output/test.faa.json", help="JSON with expected results")
    parser.add_argument("--observed", type=str, default="tests/data/output/test.faa.json", help="JSON output file from IPS6")
    parser.add_argument("--summary", action="store_true", help="Print only the summary message")
    parser.add_argument("--format", choices=["gff3", "json", "tsv", "xml", "intermediate"], default="json", help=(
        "Format of input files.\n"
        "'gff3', 'json' [default], 'tsv', or 'xml' for final output files\n"
        "or 'intermediate' to compare the temporary working files of InterProScan6."
    ))
    parser.add_argument("--applications", type=str, default=None, help="Limit the comparison to a comma-separated list of applicaitons")
    args = parser.parse_args()

    if args.format == "gff3":
        expected = parse_gff3(args.expected)
        observed = parse_gff3(args.observed)
    elif args.format == "json":
        expected = parse_json(args.expected)
        observed = parse_json(args.observed)
    elif args.format == "tsv":
        expected = parse_tsv(args.expected)
        observed = parse_tsv(args.observed)
    elif args.format == "xml":
        expected = parse_xml(args.expected)
        observed = parse_xml(args.observed)
    else:
        expected = parse_intermediate(args.expected, args.applications)
        observed = parse_intermediate(args.observed, args.applications)

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


def parse_json(iprscan_path: str):
    """Convert all matches into lists, keeping only the data we care about"""
    matches = []  # [[md5, sig_acc, loc start, loc end]]
    with open(iprscan_path, "r") as fh:
        iprsn_dict = json.load(fh)
    for protein_dict in iprsn_dict["results"]:
        md5 = protein_dict["md5"].upper()  # account for diff between iprscn 5 and 6
        for match in protein_dict["matches"]:
            sig_acc = match["signature"]["accession"]
            for loc in match["locations"]:
                matches.append([md5, sig_acc, loc["start"], loc["end"]])
    return [repr(nested) for nested in matches]


def parse_gff3(iprscan_path: str):
    matches = []
    with open(iprscan_path, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                data = line.split("\t")
                matches.append([data[0], data[1], data[3], data[4], data[8]])
    return [repr(nested) for nested in matches]


def parse_tsv(iprscan_path: str):
    matches = []
    with open(iprscan_path, "r") as fh:
        for line in fh:
            data = line.split("\t")
            matches.append([data[1], data[4], data[6], data[7]])
    return [repr(nested) for nested in matches]


def parse_xml(iprscan_path: str):
    matches = []
    tree = ET.parse(iprscan_path)
    root = tree.getroot()

    for protein in root.findall(".//protein"):
        sequence = protein.find("sequence")
        md5 = sequence.attrib.get("md5") if sequence is not None else None
        for match in protein.findall(".//match"):
            signature = match.find("signature")
            sig_ac = signature.attrib.get("ac") if signature is not None else None
            locations = match.find("locations")
            if locations is not None:
                for location in locations.findall("location"):
                    start = location.attrib.get("start")
                    end = location.attrib.get("end")
                    matches.append([md5, sig_ac, start, end])
    return [repr(match) for match in matches]


def parse_intermediate(iprscan_path: str, applications: str):
    matches = []
    apps = {app.strip().lower() for app in applications.split(",")} if applications else None
    with open(iprscan_path, "r") as fh:
        matches_dict = json.load(fh)

    for upi, data in matches_dict.items():
        for model_acc, model_data in data.items():
            library = model_data.get("signature", {}).get("signatureLibraryRelease", {}).get("library", "").lower()
            if not apps or library in apps:
                for location in model_data.get("locations", []):
                    matches.append([upi, model_acc, location.get("start"), location.get("end")])

    return [repr(nested) for nested in matches]

if __name__ == "__main__":
    main()
