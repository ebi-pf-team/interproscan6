import argparse
import json


def add_entries_info(matches_path: str, entries_path: str) -> list[dict]:
    matches2entries = []
    with open(matches_path, "r") as matches:
        with open(entries_path + ".ipr.json", "r") as fh:
            entries = json.load(fh)
            for m in matches:
                matches_info = json.loads(m)
                for info in matches_info:
                    try:
                        acc = info["signature_acc"]
                        entry = entries[acc]
                        info["interpro_annotations_desc"] = entry[0]
                        info["signature_desc"] = entry[1]
                        info["interpro_annotations_acc"] = entry[2]
                    except KeyError:
                        pass
                    matches2entries.append(info)
    return matches2entries


def add_goterms_info(matches_info: list[dict], goterm_path: str):
    matches2go = {}
    with open(goterm_path + ".ipr.json", "r") as ipr:
        with open(goterm_path + ".json", "r") as go:
            ipr2go = json.load(ipr)
            go_info = json.load(go)
            for match in matches_info:
                try:
                    acc = match["interpro_annotations_acc"]
                    go_ids = ipr2go[acc]
                    for go_id in go_ids:
                        matches2go[go_id] = go_info[go_id]
                        match["GOTERMS"] = matches2go
                except KeyError:
                    pass
    return matches_info


def add_pathways_info(matches_info, pathway_path: str):
    matches2pa = {}
    with open(pathway_path + ".ipr.json", "r") as ipr:
        with open(pathway_path + ".json", "r") as pa:
            ipr2pa = json.load(ipr)
            pa_info = json.load(pa)
            for match in matches_info:
                try:
                    acc = match["interpro_annotations_acc"]
                    pa_ids = ipr2pa[acc]
                    for pa_id in pa_ids:
                        matches2pa[pa_id] = pa_info[pa_id]
                        match["PATHWAY"] = matches2pa
                except KeyError:
                    pass
    return matches_info


def main():
    parser = argparse.ArgumentParser(
        description="complement match lookup with external files infos"
    )
    parser.add_argument("-matches", "--matches", type=str, help="Match lookup parsed")
    parser.add_argument(
        "-entries", "--entries", type=str, help="entries xref file path"
    )
    parser.add_argument("-go", "--goterms", type=str, default="", help="goterms xref file path")
    parser.add_argument("-pa", "--pathways", type=str, default="", help="pathways xref file path")
    args = parser.parse_args()

    matches_info = add_entries_info(args.matches, args.entries)

    if args.goterms:
        matches_info = add_goterms_info(matches_info, args.goterms)
    if args.pathways:
        matches_info = add_pathways_info(matches_info, args.pathways)

    print(matches_info)
    return matches_info


if __name__ == "__main__":
    main()
