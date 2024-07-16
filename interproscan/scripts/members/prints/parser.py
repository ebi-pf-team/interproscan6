import json
import sys
# Parse prints postprocessed output to standardised JSON format
# param prints_postproc: path to prints_postprocess output file
# param version: prints version number


def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0])
    print(json.dumps(parsed_results, indent=2))


def parse(prints: str) -> dict:
    rebuild = {}
    with open(prints, 'r') as file:
        prints_dict = json.load(file)
        for protein in prints_dict:
            for match in prints_dict[protein]:
                rebuild[protein] = {}
                match.pop("num_motif")
                for location in match["locations"]:
                    location.pop("evalue")
                    location.pop("model_id")
                    location["location-fragments"] = [
                        {"start": int(location["start"]), "end": int(location["end"]), "dc-status": "CONTINUOUS"}
                    ]
                rebuild[protein] = {match["accession"]: match}

    return rebuild


if __name__ == "__main__":
    main()
