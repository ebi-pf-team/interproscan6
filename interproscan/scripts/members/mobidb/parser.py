import json
import re
import sys


def parse(input_file, release):
    matches = {}
    with open(input_file, 'r') as reader:
        for line in reader:
            line_data = line.split()
            sequence_id = line_data[0]
            location_start = int(line_data[1])
            location_end = int(line_data[2])
            feature = line_data[3] if line_data[3] else ""
            info = {
                "member": "mobidb",
                "release": release,
                "location_start": location_start,
                "location_end": location_end,
                "feature": feature
            }
            try:
                matches[sequence_id].append(info)
            except KeyError:
                matches[sequence_id] = [info]

    return matches


def main():
    args = sys.argv[1:]
    matches = parse(args[0], args[1])

    print(json.dumps(matches, indent=4))


if __name__ == "__main__":
    main()
