import json
import sys


def add_entries(matches_path: str, entries_path: str) -> dict:
    matches2entries = {}
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(entries_path, "r") as fh:
        entries = json.load(fh)

    for seq_id, match_info in matches_info.items():
        for match_key, data in match_info.items():
            acc_id = match_key.split(".")[0]
            match_info[match_key]['member_db'] = data['member_db']
            match_info[match_key]['version'] = entries[data["member_db"]]['version']
            try:
                entry = entries[acc_id]
                match_info[match_key]["entry"] = {
                    "accession": entry["integrated"],
                    "name": entry["name"],
                    "description": entry["description"],
                    "type": entry["type"],
                    "version": entry["database"]["version"],
                    "member_db": entry["database"]["name"],
                    "goXRefs": [],
                    "pathwayXRefs": []
                }
            except KeyError:
                acc_id = match_key  # some accs need the '.'  , e.g. Gene3D
                try:
                    entry = entries[acc_id]
                    match_info[match_key]["entry"] = {
                        "accession": entry["integrated"],
                        "name": entry["name"],
                        "description": entry["description"],
                        "type": entry["type"],
                        "version": entry["database"]["version"],
                        "member_db": entry["database"]["name"],
                        "goXRefs": [],
                        "pathwayXRefs": []
                    }
                except KeyError:  # members with no match in entries
                    member_db_lower = data["member_db"].lower()
                    for key in entries['databases']:
                        if key.lower() == member_db_lower:
                            version = entries['databases'][key]

                    match_info[match_key]["entry"] = {
                        "accession": None,
                        "name": None,
                        "description": None,
                        "type": None,
                        "database": data["member_db"],
                        "goXRefs": [],
                        "pathwayXRefs": [],
                        "version": version
                    }

            if data["member_db"].upper() == "PANTHER":
                acc_id_family = data["accession"]
                try:
                    match_info[match_key]["entry"]["family_name"] = entries[acc_id_family]["name"]
                    match_info[match_key]["entry"]["family_type"] = entries[acc_id_family]["type"]
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
