import argparse
import json


def entries_info(matches_path: str, entries_path: str) -> dict:
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


def add_goterms_info(matches_info: dict, goterm_path: str):
    matches2go = {}
    with open(goterm_path + ".ipr.json", "r") as ipr:
        with open(goterm_path + ".json", "r") as go:
            ipr2go = json.load(ipr)
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


def add_pathways_info(matches_info: dict, pathway_path: str):
    matches2pa = {}
    with open(pathway_path + ".ipr.json", "r") as ipr:
        with open(pathway_path + ".json", "r") as pa:
            ipr2pa = json.load(ipr)
            pa_info = json.load(pa)
            for seq_id, matches in matches_info.items():
                for match in matches:
                    for domain in match["domains"]:
                        try:
                            acc = domain["interpro_annotations_acc"]
                            pa_ids = ipr2pa[acc]
                            for pa_id in pa_ids:
                                matches2pa[pa_id] = pa_info[pa_id]
                                domain["PATHWAY"] = matches2pa
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

    matches_info = entries_info(args.matches, args.entries)

    if args.goterms:
        matches_info = add_goterms_info(matches_info, args.goterms)
    if args.pathways:
        matches_info = add_pathways_info(matches_info, args.pathways)

    print(json.dumps(matches_info, indent=2))


if __name__ == "__main__":
    main()
