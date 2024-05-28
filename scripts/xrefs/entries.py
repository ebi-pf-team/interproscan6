import json
import sys


def add_entries(matches_path: str, entries_path: str) -> dict:
    matches2entries = {}
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(entries_path + ".ipr.json", "r") as fh:
        entries = json.load(fh)

    for seq_id, match_info in matches_info.items():
        for match_key, data in match_info.items():
            acc_id = match_key.split(".")[0]
            try:
                entry = entries[acc_id]
                match_info[match_key]["entry"] = {
                    "accession": entry[0] if entry[0] is not None else "-",
                    "short_name": entry[1] if entry[1] is not None else "-",
                    "name": entry[2] if entry[2] is not None else "-",
                    "description": entry[3] if entry[3] is not None else "-",
                    "type": entry[4],
                    "goXRefs": [],
                    "pathwayXRefs": []
                }
            except KeyError:
                acc_id = match_key  # some accs need the .
                try:
                    entry = entries[acc_id]
                    match_info[match_key]["entry"] = {
                        "accession": entry[0] if entry[0] is not None else "-",
                        "short_name": entry[1] if entry[1] is not None else "-",
                        "name": entry[2] if entry[2] is not None else "-",
                        "description": entry[3] if entry[3] is not None else "-",
                        "type": entry[4],
                        "goXRefs": [],
                        "pathwayXRefs": []
                    }
                except KeyError:
                    match_info[match_key]["entry"] = None
        matches2entries[seq_id] = match_info

    return matches2entries


def main():
    args = sys.argv[1:]

    matches = args[0]
    entries = args[1]

    matches2entries = add_entries(matches, entries)

    print(json.dumps(matches2entries, indent=2))


if __name__ == "__main__":
    main()
