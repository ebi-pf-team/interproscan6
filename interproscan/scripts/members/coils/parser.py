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
                seq_id = line.strip(">").strip("\n").split()[0]
                match = {"Coil": {"member_db": "Coils",
                                  "version": version, "name": "Coil",
                                  "accession":"Coil", "locations": []}}
                matches[seq_id] = match
            else:
                line = line.split()
                if len(line) > 1:
                    location_start = line[0]
                    location_end = line[1]
                    location_fragments = {"start": int(location_start),
                                          "end": int(location_end),
                                          "dc-status": "CONTINUOUS"}
                    location = {"start": location_start,
                                "end": location_end,
                                "representative": "false",
                                "location-fragments": [location_fragments]}
                    matches[seq_id]["Coil"]["locations"].append(location)
    return matches


if __name__ == '__main__':
    main()
