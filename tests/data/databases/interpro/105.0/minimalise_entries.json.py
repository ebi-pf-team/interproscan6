#!/usr/bin/env python3
# Make a smaller copy of an entries.json file containing only the
# interpro entries we need for the unit tests to reduce the storage size

import json

ENTRIES = "tests/data/interpro/entries.json" # path to the original entries file
OUTPUT = "tests/data/interpro/entries.json" # path to the output json file
MATCHES = "tests/data/channels/matches2xrefs_minimal.json" # path to the JSON file containing matches with interpro entries

def main():
    with open(ENTRIES, "r") as fh:
        entries = json.load(fh)

    with open(MATCHES, "r") as fh:
        matches = json.load(fh)

    minimalised_entries = {}

    for upi, matches_dict in matches.items():
        for model_ac, match in matches_dict.items():
            if match["signature"]["entry"]:
                interpro_ac = match["signature"]["entry"]["accession"]
                minimalised_entries[interpro_ac] = entries[interpro_ac]

    with open(OUTPUT, "w") as fh:
        json.dump(minimalised_entries, fh)

if __name__ == "__main__":
    main()
