import json
import sys
import re
# Parse prints output to standardised JSON format
# param prints_out: path to prints output file
# param hierarchy: path to prints hierarchy db
# param version: prints version number


def main():
    args = sys.argv[1:]
    hierarchy_map = parse_hierarchy(args[1])
    results = parse_prints(args[0], hierarchy_map, args[2])
    sorted_results = sort_results(results)
    selected_results = select_results(sorted_results, hierarchy_map)
    parsed_results = process_results(selected_results)
    print(json.dumps(parsed_results, indent=2))


def parse_hierarchy(hierarchy: str) -> dict:
    hierarchymap = {}
    with open(hierarchy, "r") as f:

        model_ids = []
        for line in f:
            if line.startswith(("/", "#")):
                continue
            row = line.split("|")
            if len(row) >= 3:
                model_id = row[0]
                model_ids.append(model_id)
                model_acc = row[1]
                evalue_cutoff = row[2]
                min_motif_count = row[3]
                hierarchical_rel = row[4]
                is_domain = not hierarchical_rel
                if not is_domain:
                    hierarchical_rel = hierarchical_rel.split(",")
                hierarchymap[model_id] = {
                    "model_acc": model_acc,
                    "evalue_cutoff": evalue_cutoff,
                    "min_motif_count": min_motif_count,
                    "sibling_list": hierarchical_rel,
                    "is_domain": is_domain
                }

    return hierarchymap


def parse_prints(prints_out: str, hierarchy_map: dict, version: str) -> dict:
    '''
Extracts fingerprint match info from prints output.
Sn line: protein_id
2TBN/2TBH lines: fingerprint match summary values
3TBN/3TBH lines: fingerprint motif match values
Other lines: blank or not required
    '''
    results = {}
    with open(prints_out) as f:
        for line in f:
            if line.startswith("Sn; "):
                idline = line.replace("Sn; ", "")
                idline = idline.strip("\n")
                protein_id = idline.split(" ")[0]
            if line.startswith(("2TBN", "2TBH")):
                fingerprint, nummotif, evalue, graphscan = process_2tb(line)
                # hierarchy map to get model ac
                if fingerprint in hierarchy_map:
                    model_acc = hierarchy_map[fingerprint]["model_acc"]
                match = {"accession": model_acc,
                         "name": fingerprint,
                         "member_db": "PRINTS",
                         "version": version,
                         "evalue": float(evalue),
                         "num_motif": nummotif,
                         "graphscan": graphscan,
                         "model-ac": model_acc,
                         "locations": []}
                if protein_id in results:
                    results[protein_id].append(match)
                else:
                    results[protein_id] = [match]

            if line.startswith(("3TBN", "3TBH")):
                motifname, motifnum, idscore, pvalue, pos, end = process_3tb(line)
                for match in results[protein_id]:
                    if motifname in match["name"]:
                        match["locations"].append({
                            "motifNumber": int(motifnum),
                            "pvalue": pvalue,
                            "score": idscore,
                            "start": pos,
                            "end": end,
                            "representative": "false",
                            "evalue": match["evalue"],
                            "model_id": match["accession"]})

    return results


def process_2tb(line):
    line = re.sub(r"\s+", "\t", line)
    line = line.split("\t")
    fingerprint = line[1]
    num_motifs = line[2] + line[3] + line[4]
    evalue = line[9]
    graphscan = line[10]
    return fingerprint, num_motifs, evalue, graphscan


def process_3tb(line):
    line = re.sub(r"\s+", "\t", line)
    line = line.split("\t")
    motifname = line[1]
    motifnum = int(line[2])
    idscore = line[5]
    pvalue = line[7]
    length = line[9]
    pos = line[11]
    end = int(pos) + int(length) - 1
    return motifname, motifnum, idscore, pvalue, pos, end


def sort_results(results: dict) -> dict:
    for protein_id in results:
        results[protein_id].sort(key=lambda d: (d["evalue"]))
        for model in results[protein_id]:
            model["locations"].sort(key=lambda d:
            (d["model_id"], d["motifNumber"], d["start"], d["end"]))
    return results


def select_results(sorted_res: dict, hierarchy: dict) -> dict:
    # from hierarchymap, get evalue for id
    pass_matches = {}
    for protein_id in sorted_res:
        for match in sorted_res[protein_id]:
            fingerprint = match["name"]
            if fingerprint in hierarchy.keys():
                # process sorted matches step
                # filter by min motif
                min_motif = hierarchy[fingerprint]["min_motif_count"]
                cutoff = float(hierarchy[fingerprint]["evalue_cutoff"])
                if match["evalue"] <= cutoff and match["num_motif"] > min_motif:
                    pass
                else:
                    continue

                if protein_id in pass_matches:
                    pass_matches[protein_id].append(match)
                else:
                    pass_matches[protein_id] = [match]

        if protein_id not in pass_matches:
            pass_matches[protein_id] = []

    return pass_matches


def process_results(prints_dict: dict) -> dict:
    rebuild = {}
    for protein in prints_dict:
        for match in prints_dict[protein]:
            rebuild[protein] = {}
            match.pop("num_motif")
            for location in match["locations"]:
                location.pop("evalue")
                location.pop("model_id")
                location["location-fragments"] = [
                    {"start": int(location["start"]),
                     "end": int(location["end"]),
                     "dc-status": "CONTINUOUS"}
                ]
            rebuild[protein] = {match["accession"]: match}

    return rebuild


if __name__ == "__main__":
    main()
