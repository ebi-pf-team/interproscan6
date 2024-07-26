import json
import sys


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

        for line in f:
            if line.startswith(("/", "#")):
                continue
            row = line.split("|")
            if len(row) >= 3:
                model_id = row[0]
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
    1TBH line: fingerprint description
    2TBH lines: fingerprint match summary values
    3TBH lines: fingerprint motif match values
    Other lines: blank or not required
    '''
    matches = {}
    protein_id = ""
    version = prints_out.split("._.")[0]
    with open(prints_out) as f:
        for line in f:
            if line.startswith("Sn; "):
                protein_id = line.split(maxsplit=2)[1]
                protein_hits = {}
                matches[protein_id] = {}
            elif line.startswith("1TBH"):
                #  keep track of hit motifs and motif names
                #print(line)
                tbh = line.split()
                motif = tbh[1]
                protein_hits[motif] = {}
                desc = " ".join(tbh[3:len(tbh) - 1])
                motif_acc = tbh[len(tbh)-1]
                # key value store of hit motif id and motif name
                protein_hits[motif]["desc"] = desc
                protein_hits[motif]["acc"] = motif_acc
                #print(protein_hits[motif])

            elif line.startswith("2TBH"):
                # fingerprint = motif
                fingerprint, nummotif, evalue, graphscan = process_2tb(line)
                # hierarchy map to get model ac
                if fingerprint in hierarchy_map:
                    model_acc = hierarchy_map[fingerprint]["model_acc"]
                    min_motif = hierarchy_map[fingerprint]["min_motif_count"]
                    cutoff = float(hierarchy_map[fingerprint]["evalue_cutoff"])

                    if evalue <= cutoff and nummotif > min_motif:
                        match = {"accession": model_acc,
                                 "name": fingerprint,
                                 "member_db": "PRINTS",
                                 "version": version,
                                 "evalue": evalue,
                                 "graphscan": graphscan,
                                 "model-ac": model_acc,
                                 "description": protein_hits[fingerprint]["desc"],
                                 "locations": []}
                        if model_acc not in matches[protein_id]:
                            matches[protein_id][model_acc] = match

                else:
                    continue

            elif line.startswith("3TBH"):
                motifname, motifnum, idscore, pvalue, pos, end = process_3tb(
                    line)
                acc = protein_hits[motifname]["acc"]
                if matches[protein_id]:
                    #print(matches[protein_id])
                    #sprint(acc)
                    matches[protein_id][acc]["locations"].append({
                    "motifNumber": int(motifnum),
                    "pvalue": pvalue,
                    "score": idscore,
                    "start": pos,
                    "end": end,
                    "representative": "false",
                    "location-fragments": [{
                        "start": pos,
                        "end": end,
                        "dc-status": "CONTINUOUS"}]})

    return matches


def process_2tb(line):
    line = line.split()
    fingerprint = line[1]
    num_motifs = line[2]
    evalue = float(line[9])
    graphscan = line[10]
    return fingerprint, num_motifs, evalue, graphscan


def process_3tb(line):
    line = line.split()
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

