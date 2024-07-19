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
                seq_id = line[1:].strip(">").strip("\n")
                print(seq_id)
                match = {"Coil": {"member_db": "Coils", "version": version,
                             "accession": "Coil", "name": "Coil",
                             "model-ac": "Coil",
                             "locations": []}}
                matches[seq_id] = match

            elif not line.startswith((">", "//", "\n")):
                line = line.split()
                if line:
                    location_start = line[0]
                    location_end = line[1]
                    location = {"start": location_start, "end": location_end}
                    matches[seq_id]["Coil"]["locations"].append(location)
    return matches


if __name__ == '__main__':
    main()
    