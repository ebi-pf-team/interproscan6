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
    pass_matches = []
    with open(prints_out) as f:
        for line in f:
            if line.startswith("Sn; "):
                idline = line.replace("Sn; ", "")
                idline = idline.strip("\n")
                print(idline)
                protein_id = idline.split(" ")[0]
                print(protein_id)
                # filtering based on evalue
                # see proteinidtofiltered match
            if line.startswith("2TBN"):
                # regex match spaces and replace
                line = re.sub(r"\s+", "\t", line)
                line = line.split("\t")
                num_motifs = line[2]+line[3]+line[4]
                evalue = line[9]
                fingerprint = line[1]

                # sorting matches with "sortedRawMatches"

                # from hierarchymap, get evalue for id
                if fingerprint in hierarchy_map.keys():
                    # process sorted matches step
                    cutoff = hierarchy_map[fingerprint]["evalue_cutoff"]
                    if evalue <= cutoff:
                        modelpass = True
                        pass_matches.append({"protein_id": protein_id, "fingerprint": fingerprint})
                else:
                    print("model not found")

    return pass_matches


def parse_hierarchy(hierarchy:str) -> dict:
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
                sibling_list = row[4].strip("\n")
                is_domain = False
                if sibling_list == "":
                    is_domain = True
                else:
                    sibling_list = list([sibling_list.replace(",", ", ")])
                hierarchymap[model_id] = {
                "model_acc": model_acc,
                "evalue_cutoff": evalue_cutoff,
                "min_motif_count": min_motif_count,
                "sibling_list": sibling_list,
                "is_domain": is_domain
                }
        #print(len(model_ids))

    return hierarchymap


if __name__ == "__main__":
    main()
    