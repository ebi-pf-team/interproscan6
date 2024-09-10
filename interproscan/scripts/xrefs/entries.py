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
    "SIGNALP": "SignalP",
    "SIGNALP_EUK": "SignalP_Euk",
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
            databases_versions = entries["databases"]
            entries_info = entries['entries']
            acc_id = match_key.split(".")[0]
            member_db = data["member_db"].upper()

            match_info[match_key]['member_db'] = member_db
            match_info[match_key]['library'] = map_databases[member_db]
            if 'version' not in match_info[match_key]:  # mobidb, signalp and tmhmm get version from members.config
                match_info[match_key]['version'] = databases_versions[map_databases[member_db]]

            entry = entries_info.get(acc_id) or entries_info.get(match_key)
            if entry:
                interpro_key = entry['integrated']
                match_info[match_key]["entry"] = {
                    "accession": interpro_key,
                    "name": entry["name"],
                    "description": entry["description"],
                    "type": entry["type"],
                    "representative": entry["representative"],
                    "member_db": member_db,
                    "goXRefs": [],
                    "pathwayXRefs": []
                }
                ipr_info = entries_info.get(interpro_key)
                if ipr_info:
                    match_info[match_key]["entry"]["ipr_name"] = ipr_info["name"]
                    match_info[match_key]["entry"]["ipr_description"] = ipr_info["description"]
                    match_info[match_key]["entry"]["ipr_type"] = ipr_info["type"]

            else:  # members with no match in entries. e.g. PIRSR
                match_info[match_key]["entry"] = {
                    "accession": None,
                    "name": None,
                    "description": None,
                    "type": None,
                    "database": member_db,
                    "goXRefs": [],
                    "pathwayXRefs": []
                }

            if member_db == "PANTHER":
                acc_subfamily = data["accession"]
                try:
                    match_info[match_key]["entry"]["subfamily_name"] = entries_info[acc_subfamily]["name"]
                    match_info[match_key]["entry"]["subfamily_description"] = entries_info[acc_subfamily]["description"]
                    match_info[match_key]["entry"]["subfamily_type"] = entries_info[acc_subfamily]["type"]
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
