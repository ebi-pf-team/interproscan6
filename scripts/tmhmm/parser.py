import json
import sys

def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def parse(tmrgff: str, version: str) -> dict:
    results = {}

    with open(tmrgff, "r") as f:
        for line in f:
            ids = list(results.keys())
            if line.startswith(("/", "#")):
                continue
            protein_id = line.split("\t")[0]
            location = line.split("\t")[1]
            start = line.split("\t")[2]
            end = line.split("\t")[3]
            locations = {"location": location, "start": start, "end": end}
            if protein_id in results:
                results[protein_id] = {
                    "transmembrane_prediction": {
                        "member_db": "DeepTMHMM",
                        "version": version,
                        "locations": []
                    }
                }
            if protein_id in results:
            results[protein_id]["transmembrane_prediction"]["locations"].append(locations)

    return results


if __name__ == "__main__":
    main()
