import json
import sys


def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def parse(coils_out: str, version: float) -> dict:
    matches = {}
    with open(coils_out) as coils_file:
        for line in coils_file:
            if line.startswith(">"):
                seq_id = line.strip(">").split(maxsplit=1)[0]
                match = {"Coil": {"member_db": "Coils",
                                  "version": version, "name": "Coil",
                                  "accession":"Coil", "locations": []}}
            else:
                line = line.split()
                if len(line) > 1:
                    matches[seq_id] = match
                    location_start = int(line[0])
                    location_end = int(line[1])
                    location_fragments = {"start": location_start,
                                          "end": location_end,
                                          "dc-status": "CONTINUOUS"}
                    location = {"start": location_start,
                                "end": location_end,
                                "representative": "false",
                                "location-fragments": [location_fragments]}
                    matches[seq_id]["Coil"]["locations"].append(location)
    return matches


if __name__ == '__main__':
    main()
