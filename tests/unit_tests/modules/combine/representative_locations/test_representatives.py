#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test for mismatches in identifying representative domains and locations
# Default args assuming this is running in the root of the repo

import ast
import argparse
import difflib
import json

def main():
    parser = argparse.ArgumentParser(prog="IPS6_test_representatives", description="Check classification of representative domains and families", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--expected", type=str, default="tests/data/channels/matches_with_representative.json", help="JSON with expected results")
    parser.add_argument("--observed", type=str, default="tests/data/channels/matches_with_representative.json", help="JSON with generated results")
    parser.add_argument("--summary", action="store_true", help="Print only the summary message")
    args = parser.parse_args()

    expected = parse_json(args.expected)
    observed = parse_json(args.observed)

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


def parse_json(json_path: str):
    matches = []
    with open(json_path, "r") as fh:
        matches_dict = json.load(fh)
    for upi, data in matches_dict.items():
        for model_acc in data:
            if (data[model_acc]["representativeInfo"]):
                if (data[model_acc]["representativeInfo"]["type"]):
                    for location in data[model_acc]["locations"]:
                        matches.append([upi, model_acc, data[model_acc]["representativeInfo"]["type"], str(location["representative"])])
    return [repr(nested) for nested in matches]

if __name__ == "__main__":
    main()
