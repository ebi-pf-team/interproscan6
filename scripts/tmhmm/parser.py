import json
import sys

def main():
    args = sys.argv[1:]
    #parsed_results = parse(args[0], args[1])
    parsed_results = parse(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def parse(resultsmd: str, version: str) -> dict:
#def parse(resultsmd: str, version: str) -> dict:
    results = {}

    with open(resultsmd, "r") as f:
        for line in f.readlines():
            ids = list(results.keys())
            if line.startswith("/"):
                continue
            if line.startswith("#"):
                continue
            id = line.split("\t")[0]
            location = line.split("\t")[1]
            start = line.split("\t")[2]
            end = line.split("\t")[3]
            locations = {"location": location, "start": start, "end": end}
            if id in ids:
                results[id]["transmembrane_prediction"]["locations"].append(locations)
            else:
                results[id] = {
                    "transmembrane_prediction": {
                        "member_db": "DeepTMHMM",
                        "version": version,
                        "locations": [locations]
                    }
                }

    return results


if __name__ == "__main__":
    main()
