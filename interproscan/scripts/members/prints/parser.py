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
                motif, desc, motif_acc = process_1tb(line)
                # key value store of hit motif id and motif name
                protein_hits[motif] = {}
                protein_hits[motif]["desc"] = desc
                protein_hits[motif]["acc"] = motif_acc

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


def process_1tb(line):
    line_pattern = re.compile(r"^(\w+)\s+(\w+)\s+([\d+\.]*\d+e[+-]?\d+|[\d\.]+)\s+([A-Za-z0-9\s\-\/\(\)]+?)\s+(\w+)\s*$")
    # collects groups in order:
    # line, FingerPrint, Evalue, description, accession
    rematch = line_pattern.match(line)
    motif = rematch.group(2)
    desc = rematch.group(4)
    acc = rematch.group(5)
    return motif, desc, acc


def process_2tb(line):
    line_pattern = re.compile(r"^(\w+)\s+(\w+)\s+(\d+)\s+(of\s+\d+)\s+([\d\.]+)\s+([\d\.]+)\s+(\d+)\s+([\d+\.]*\d+e[+-]?\d+|[\d\.]+)\s+([\d\.]*\d+e[+-]?\d+|[\d\.])\s+([Ii.]+)\s*$")
    # collects groups in order:
    # line, FingerPrint, No.Motifs, SumId, AveId,
    # ProfScore, Ppvalue,Evalue, GraphScan
    rematch = line_pattern.match(line)
    fingerprint = rematch.group(2)
    num_motifs = rematch.group(3)
    evalue = float(rematch.group(9))
    graphscan = rematch.group(10)
    return fingerprint, num_motifs, evalue, graphscan


def process_3tb(line):
    line_pattern = re.compile(r"^(\w+)\s+(\w+)\s+(\d+)\s+(of\s+\d+)\s+([\d\.]+)\s+(\d+)\s+([\d+\.]*\d+e[+-]?\d+|[\d\.]+)\s+(#*[a-zA-Z]+#*)\s+(\d+)\s+(\d+)\s+(\d+)\s*(\d)\s*$")
    # collects groups in order:
    # line, MotifName, No.Motifs, of total number of motifs,
    # IdScore, PfScore, Pvalue, Sequence, Length, Low, Position, High
    rematch = line_pattern.match(line)
    motifname = rematch.group(2)
    motifnum = int(rematch.group(3))
    idscore = rematch.group(5)
    pvalue = float(rematch.group(7))
    sequence = rematch.group(8)
    length = rematch.group(9)
    pos = int(rematch.group(11))
    end = pos + int(length) - 1
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

