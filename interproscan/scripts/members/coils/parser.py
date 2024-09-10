import json
import sys


def main():
    """CL input:
    0. Str repr of path to the output from rpsblast
    1. Str repr of path to the output file
    """
    args = sys.argv[1:]
    parsed_results = parse(args[0])
    with open(args[1], "w") as fh:
        json.dump(parsed_results, fh)


def parse(coils_out: str) -> dict:
    matches = {}
    with open(coils_out) as coils_file:
        for line in coils_file:
            if line.startswith(">"):
                seq_id = line.strip(">").split(maxsplit=1)[0]
                match = {"Coil": {"member_db": "Coils",
                                  "name": "Coil",
                                  "accession": "Coil", "locations": []}}
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
                                "location-fragments": [location_fragments]}
                    matches[seq_id]["Coil"]["locations"].append(location)
    return matches


if __name__ == '__main__':
    main()
