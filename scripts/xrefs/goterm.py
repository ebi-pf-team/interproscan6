import json
import sys


def add_goterms(matches_path: str, goterm_path: str):
    matches2go = {}
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(goterm_path + ".ipr.json", "r") as ipr:
        ipr2go = json.load(ipr)
    with open(goterm_path + ".json", "r") as go:
        go_info = json.load(go)

    for seq_id, matches in matches_info.items():
        for match in matches:
            for domain in match["domains"]:
                try:
                    acc = domain["interpro_annotations_acc"]
                    go_ids = ipr2go[acc]
                    for go_id in go_ids:
                        matches2go[go_id] = go_info[go_id]
                        domain["GOTERMS"] = matches2go
                except KeyError:
                    pass
    return matches_info


def main():
    args = sys.argv[1:]

    matches = args[0]
    pathways = args[1]

    matches2goterms = add_goterms(matches, pathways)

    print(json.dumps(matches2goterms, indent=2))


if __name__ == "__main__":
    main()
