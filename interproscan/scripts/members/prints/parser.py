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
    results = parse_prints(args[0], hierarchy_map)
    print(json.dumps(results, indent=2))


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


def parse_prints(prints_out: str, hierarchy_map: dict) -> dict:
    '''
    Extract fingerprint match info from prints output.

    Lines in the output file:
    Sn line: protein_id
    2TBN/2TBH lines: fingerprint match summary values
    3TBN/3TBH lines: fingerprint motif match values
    Other lines: blank or not required
    '''
    matches = {}
    protein_id = ""

    with open(prints_out) as f:
        for line in f:
            if line.startswith("Sn; "):
                protein_id = line.split()[1]
                protein_hits = {}
                matches[protein_id] = {}
            if line.startswith("1TBH"):
                #  keep track of hit motifs and motif names
                tbh = line.split()
                motif = tbh[1]
                desc = " ".join(tbh[3:len(tbh)-1])
                # key value store of hit motif id and motif name
                protein_hits[motif] = desc
            if line.startswith("2TBH"):
                fingerprint, nummotif, evalue, graphscan = process_2tb(line)
                # hierarchy map to get model ac
                if fingerprint in hierarchy_map:
                    model_acc = hierarchy_map[fingerprint]["model_acc"]
                    min_motif = hierarchy_map[fingerprint]["min_motif_count"]
                    cutoff = float(hierarchy_map[fingerprint]["evalue_cutoff"])
                else:
                    continue
                match = {"accession": model_acc,
                         "name": fingerprint,
                         "member_db": "PRINTS",
                         "version": prints_out.split("._.")[0],
                         "evalue": float(evalue),
                         "num_motif": nummotif,
                         "graphscan": graphscan,
                         "model-ac": model_acc,
                         "description": desc,
                         "locations": []}
                         # filter by min motif
                if match["evalue"] <= cutoff and match["num_motif"] > min_motif:
                    match.pop("num_motif")

                if model_acc not in matches[protein_id]:
                    matches[protein_id][model_acc] = match

            elif line.startswith("3TBH"):
                motifname, motifnum, idscore, pvalue, pos, end = process_3tb(line)
                # compile key value store for model id, model name from 1tb
                # use to assign location
                for hit in protein_hits:
                    if hit == motifname:
                        matches[protein_id][model_acc]["locations"].append({
                            "motifNumber": int(motifnum),
                            "pvalue": pvalue,
                            "score": idscore,
                            "start": pos,
                            "end": end,
                            "representative": "false",
                            "location-fragments": [{"start": pos,
                                 "end": end,
                                 "dc-status": "CONTINUOUS"}]})
                        matches[protein_id][model_acc]["description"] = protein_hits[hit]

    return matches


def process_2tb(line):
    line = line.split()
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
    sequence = line[8]
    length = line[9]
    pos = line[11]
    # corrects for pos merging with next column
    if len(pos) > 5:
        pos = pos[:6]
    end = int(pos) + int(length) - 1
    # corrects for motif overhanging start
    if int(pos) < 1:
        pos = 1
    # corrects for motif overhanging end
    if sequence.endswith("#"):
        motiflength = len(sequence)
        indexcheck = motiflength - 1
        while sequence[indexcheck] == "#":
            indexcheck -= 1
        end = end - (motiflength - indexcheck) + 1
    return motifname, motifnum, idscore, pvalue, pos, end


if __name__ == "__main__":
    main()

