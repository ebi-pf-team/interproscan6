"""Use to build a entries.json for the unit tests that only
contains the necessary InterPro Entries.
Initially copy across the entire entries.json file to
tests/unit_tests/test_inputs/xrefs/entries/entries.json. 
Then run this python file from the root of the repo. It 
will trim down the entries.json.
"""

import json

ENTRIES = "tests/unit_tests/test_inputs/xrefs/entries/entries.json"
MATCH_ENTRIES = "tests/unit_tests/test_inputs/xrefs/match_entries.json"

with open(ENTRIES, "r") as fh:
    entries = json.load(fh)
with open(MATCH_ENTRIES, "r") as fh:
    results = json.load(fh)

entries_of_interest = set()
for seq_id, matches in results.items():
    for sig_acc, sig_data in matches.items():
        if sig_data["entry"]["accession"]:
            entries_of_interest.add(sig_data["entry"]["accession"])

parsed_entries = {"databases": entries["databases"], "entries": {}}
for ent in entries_of_interest:
    parsed_entries["entries"][ent] = entries["entries"][ent]

with open(ENTRIES, "w") as fh:
    json.dump(parsed_entries, fh)
