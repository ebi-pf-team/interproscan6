import json
import sys
import re


# postprocess prints
# param prints_out: path to prints output file
# param hierarchy: path to prints hierarchy db

def main():
    args = sys.argv[1:]
    parsed_results = postprocess(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def postprocess(prints_out: str, hierarchy: str) -> dict:
    results = {}
    hierarchy_map = parse_hierarchy(hierarchy)
    allprintsids = list(hierarchy_map.keys())
    results = parse_prints(prints_out, hierarchy_map)
    sorted_results = sort_results(results)


def sort_results(results: dict) -> dict:
    for protein in results:
        print(protein)
        # sort models by evalue
        # aka sort order of dictionaries in a list by one of the values in the dictionary
        results[protein].sort(key=lambda d: (d["evalue"]))
        for model in results[protein]:
            model["locations"].sort(key=lambda d: (d["model_id"], d["motifnumber"], d["start"], d["end"]))
            #comp = compare_prints_raw_match(model, model)
            #print(comp)
            print(model["evalue"])
            print(model)
        #    model["evalue"] = f"{model['evalue']:.1e}"


def parse_prints(prints_out: str, hierarchy_map: str) -> dict:
    results = {}
    with open(prints_out) as f:
        for line in f:
            if line.startswith("Sn; "):
                idline = line.replace("Sn; ", "")
                idline = idline.strip("\n")
                protein_id = idline.split(" ")[0]
            if line.startswith(("2TBN", "2TBH")):
                fingerprint, nummotif, sumid, aveid, profscore, ppvalue, evalue, graphscan = process_2tb(line)
                # hierarchy map to get model ac
                if fingerprint in hierarchy_map:
                    model_acc = hierarchy_map[fingerprint]["model_acc"]

                prints = {"model_acc": model_acc, "model_name": fingerprint, "evalue": float(evalue), "graphscan": graphscan,
                          "locations": []}
                if protein_id in results:
                    results[protein_id].append(prints)
                else:
                    results[protein_id] = [prints]

            if line.startswith(("3TBN", "3TBH")):
                motifname, nummotif, idscore, pfscore, pvalue, sequence, length, low, pos, high, end = process_3tb(line)
                for i in results[protein_id]:
                    if motifname in i["model_name"]:
                        i["locations"].append(
                            {"start": pos, "end": end, "representative": False, "pvalue": pvalue, "score": idscore,
                             "motifnumber": nummotif, "evalue": i["evalue"], "model_id": i["model_acc"]})
                        #print(i)
    return results


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
                hierarchical_rel = row[4].strip("\n")
                is_domain = False
                if hierarchical_rel == "":
                    is_domain = True
                else:
                    hierarchical_rel = list([hierarchical_rel.replace(",", ", ")])
                hierarchymap[model_id] = {
                    "model_acc": model_acc,
                    "evalue_cutoff": evalue_cutoff,
                    "min_motif_count": min_motif_count,
                    "sibling_list": hierarchical_rel,
                    "is_domain": is_domain
                }

    return hierarchymap


def process_2tb(line):
    line = re.sub(r"\s+", "\t", line)
    line = line.split("\t")
    fingerprint = line[1]
    num_motifs = line[2] + line[3] + line[4]
    sumid = line[5]
    aveid = line[6]
    profscore = line[7]
    ppvalue = line[8]
    evalue = line[9]
    graphscan = line[10]
    return fingerprint, num_motifs, sumid, aveid, profscore, ppvalue, evalue, graphscan


def process_3tb(line):
    line = re.sub(r"\s+", "\t", line)
    line = line.split("\t")
    motifname = line[1]
    nummotif = line[2]
    idscore = line[5]
    pfscore = line[6]
    pvalue = line[7]
    sequence = line[8]
    length = line[9]
    low = line[10]
    pos = line[11]
    high = line[12]
    end = int(pos) + int(length)
    return motifname, nummotif, idscore, pfscore, pvalue, sequence, length, low, pos, high, end


def parse_model(matches):
    match_info
    if line.startswith("Sn; "):
        idline = line.replace("Sn; ", "")
        idline = idline.strip("\n")
        protein_id = idline.split(" ")[0]
    if line.startswith("2TBN"):
        fingerprint, nummotif, sumid, aveid, profscore, ppvalue, evalue, graphscan = process_2tb(line)
        model_acc = ""
        match_info = {"protein_id": protein_id, "evalue": evalue, "graphscan": graphscan, "model_acc": model_acc,
                      "locations": []}
    return match_info


def filter(prints_out: str, hierarchy: str) -> dict:
    # from hierarchymap, get evalue for id
    if fingerprint in hierarchy_map.keys():
        # process sorted matches step
        # filter by min motif
        min_motif = hierarchy_map[fingerprint]["min_motif_count"]
        if num_motif < min_motif:
            modelpass = True
        if hierarchy_map[fingerprint]["is_domain"]:
            # no filtering required
            modelpass = True
            pass_matches.append({"protein_id": protein_id, "fingerprint": fingerprint})
        else:
            cutoff = hierarchy_map[fingerprint]["evalue_cutoff"]
            if evalue <= cutoff:
                modelpass = True
                pass_matches.append({"protein_id": protein_id, "fingerprint": fingerprint})
            # hierarchy limitations check
    else:
        print("model not found")

    return pass_matches


if __name__ == "__main__":
    main()
