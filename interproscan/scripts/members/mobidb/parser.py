import json
import sys


def parse(input_file, release):
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
                        "release": release,
                        "name": "disorder_prediction",
                        "description": "consensus disorder prediction",
                        "locations": []
                    }
                }
    return matches


def main():
    args = sys.argv[1:]
    matches = parse(args[0], args[1])

    print(json.dumps(matches, indent=4))


if __name__ == "__main__":
    main()
