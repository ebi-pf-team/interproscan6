import json
import sys


def add_pathways(matches_path: str, pathway_path: str):
    matches2pa = {}
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(pathway_path + ".ipr.json", "r") as ipr:
        ipr2pa = json.load(ipr)
    with open(pathway_path + ".json", "r") as pa:
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
    args = sys.argv[1:]

    matches = args[0]
    pathways = args[1]

    matches2pathways = add_pathways(matches, pathways)

    print(json.dumps(matches2pathways, indent=2))


if __name__ == "__main__":
    main()
