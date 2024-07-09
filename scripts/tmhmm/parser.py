import json
import sys
# Parse tmhmm output to standardised JSON format
# param tmrgff: path to deep_tmhmm output file ('TMRs.gff3')
# param version: deep_tmhmm version number, e.g. '2.0c'

def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def parse(tmrgff: str, version: str) -> dict:
    results = {}

    with open(tmrgff, "r") as f:
        for line in f:
            if line.startswith(("/", "#")):
                continue
            protein_id = line.split("\t")[0]
            location_tag = line.split("\t")[1]
            start = line.split("\t")[2]
            end = line.split("\t")[3]
            location = {"location_tag": location_tag, "start": start, "end": end}
            if protein_id in results:
                results[protein_id]["transmembrane_prediction"]["locations"].append(location)
            else:
                results[protein_id] = {
                "transmembrane_prediction": {
                    "member_db": "DeepTMHMM",
                    "version": version,
                    "locations": [location]
                }
            }

    return results


if __name__ == "__main__":
    main()
