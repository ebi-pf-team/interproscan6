import json
import sys

map_databases = {
    "SFLD": "SFLD",
    "PRINTS": "PRINTS",
    "PFAM": "Pfam",
    "INTERPRO": "InterPro",
    "CDD": "CDD",
    "PROSITE_PROFILES": "PROSITE profiles",
    "NCBIFAM": "NCBIfam",
    "PROSITE_PATTERNS": "PROSITE patterns",
    "HAMAP": "HAMAP",
    "SMART": "SMART",
    "PIRSF": "PIRSF",
    "PANTHER": "PANTHER",
    "GENE3D": "CATH-Gene3D",
    "SUPERFAMILY": "SUPERFAMILY",
    "ANTIFAM": "AntiFam",
    "FUNFAM": "FunFam",
    "MOBIDB": "MobiDB Lite",
    "PHOBIUS": "Phobius",
    "SIGNALP": ["SignalP_Euk", "SignalP_Gram_positive", "SignalP_Gram_negative"],
    "TMHMM": "TMHMM",
    "COILS": "COILS",
    "PIRSR": "PIRSR"
}


def add_entries(matches_path: str, entries_path: str) -> dict:
    matches2entries = {}
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(entries_path, "r") as fh:
        entries = json.load(fh)

    for seq_id, match_info in matches_info.items():
        for match_key, data in match_info.items():
            acc_id = match_key.split(".")[0]
            databases_versions = entries["databases"]
            member_db = data["member_db"].upper()
            version = databases_versions[map_databases[member_db]]
            match_info[match_key]['member_db'] = member_db
            match_info[match_key]['library'] = map_databases[member_db]
            match_info[match_key]['version'] = version
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

            if member_db == "PANTHER":
                acc_subfamily = data["accession"]
                try:
                    match_info[match_key]["entry"]["subfamily_name"] = entries[acc_subfamily]["name"]
                    match_info[match_key]["entry"]["subfamily_type"] = entries[acc_subfamily]["type"]
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
