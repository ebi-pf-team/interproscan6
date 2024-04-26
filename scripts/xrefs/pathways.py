import json
import sys


def add_pathways(matches_path: str, pathway_path: str):
    with open(matches_path, "r") as matches:
        matches_info = json.load(matches)
    with open(pathway_path + ".ipr.json", "r") as ipr:
        ipr2pa = json.load(ipr)
    with open(pathway_path + ".json", "r") as pa:
        pa_info = json.load(pa)

    for seq_id, match_info in matches_info.items():
        for match_key, data in match_info.items():
            if data["entry"]:
                ipr_id = data["entry"]["accession"]
                try:
                    pa_ids = ipr2pa[ipr_id]
                    for pa_id in pa_ids:
                        match_info[match_key]["entry"]["goXRefs"].append({pa_id: pa_info[pa_id]})
                except KeyError:
                    pass
        matches_info[seq_id] = match_info

    return matches_info


def main():
    args = sys.argv[1:]

    matches = args[0]
    pathways = args[1]

    matches2pathways = add_pathways(matches, pathways)

    print(json.dumps(matches2pathways, indent=2))


if __name__ == "__main__":
    main()
