import json
import sys


def add_goterms(matches_path: str, goterm_path: str):
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(goterm_path + ".ipr.json", "r") as ipr:
        ipr2go = json.load(ipr)
    with open(goterm_path + ".json", "r") as go:
        go_info = json.load(go)

    db_pattern = {"P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"}
    for seq_id, match_info in matches_info.items():
        for match_key, data in match_info.items():
            if data["entry"]:
                ipr_id = data["entry"]["accession"]
                try:
                    go_ids = ipr2go[ipr_id]
                    for go_id in go_ids:
                        go_dict = {
                            "name": go_info[go_id][0],
                            "databaseName": "GO",
                            "category": db_pattern[go_info[go_id][1]],
                            "id": go_id
                        }
                        match_info[match_key]["entry"]["goXRefs"].append(go_dict)
                except KeyError:
                    pass
        matches_info[seq_id] = match_info

    return matches_info


def main():
    args = sys.argv[1:]

    matches = args[0]
    pathways = args[1]

    matches2goterms = add_goterms(matches, pathways)

    print(json.dumps(matches2goterms, indent=2))


if __name__ == "__main__":
    main()
