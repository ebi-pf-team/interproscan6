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
            if data["member_db"] == "mobidb":
                match_info[match_key]["entry"] = None
            else:
                try:
                    entry = entries[acc_id]
                    match_info[match_key]["entry"] = {
                        "accession": entry[0],
                        "short_name": entry[1],
                        "name": entry[2],
                        "description": entry[3],
                        "type": entry[4],
                        "goXRefs": [],
                        "pathwayXRefs": []
                    }
                except KeyError:
                    acc_id = match_key  # some accs need the '.'  , e.g. Gene3D
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

            if data["member_db"].upper() == "PANTHER":
                acc_id_family = data["accession"]
                try:
                    match_info[match_key]["entry"]["family_name"] = entries[acc_id_family][3]
                    match_info[match_key]["entry"]["family_type"] = entries[acc_id_family][4]
                except KeyError:
                    pass

        matches2entries[seq_id] = match_info

    return matches2entries


def main():
    """CL input:
    0. Str repr of the path to the internal IPS6 JSON file
    1. Str repr of the path to the XREFS entries JSON file
    2. Str repr of the path for the output file"""
    args = sys.argv[1:]

    matches = args[0]
    entries = args[1]

    matches2entries = add_entries(matches, entries)
    with open(args[2], "w") as fh:
        json.dump(matches2entries, fh)


if __name__ == "__main__":
    main()
