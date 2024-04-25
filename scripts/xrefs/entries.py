import json
import sys


def add_entries(matches_path: str, entries_path: str) -> dict:
    matches2entries = {}
    with open(matches_path, "r") as matches:
        with open(entries_path + ".ipr.json", "r") as fh:
            entries = json.load(fh)
            matches_info = json.load(matches)
            for seq_id, match_info in matches_info.items():
                for match in match_info:
                    for domain in match["domains"]:
                        acc_id = match["accession"].split(".")[0]
                        try:
                            entry = entries[acc_id]
                            domain["interpro_annotations_desc"] = entry[0]
                            domain["signature_desc"] = entry[1]
                            domain["interpro_annotations_acc"] = entry[2]
                        except KeyError:
                            pass
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
