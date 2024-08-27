import json
import sys


def parse(input_file: str) -> dict:
    matches = {}
    with open(input_file, 'r') as reader:
        for line in reader:
            line_data = line.split()
            sequence_id = line_data[0]
            start = int(line_data[1])
            end = int(line_data[2])
            feature = line_data[3] if line_data[3] else ""
            try:
                matches[sequence_id]["mobidb"]["locations"].append({
                    "start": start,
                    "end": end,
                    "sequence-feature": feature
                })
            except KeyError:
                matches[sequence_id] = {
                    "mobidb": {
                        "member_db": "mobidb",
                        "accession": "mobidb",
                        "name": "disorder_prediction",
                        "description": "consensus disorder prediction",
                        "locations": [{
                            "start": start,
                            "end": end,
                            "sequence-feature": feature
                        }]
                    }
                }
    return matches


def main():
    """CL input:
    0. Str repr of path to mobidb-lite (idpred) output
    1. release number
    2. Str repr of path to output file"""
    args = sys.argv[1:]
    matches = parse(args[0])

    with open(args[2], "w") as fh:
        json.dump(matches, fh)


if __name__ == "__main__":
    main()
