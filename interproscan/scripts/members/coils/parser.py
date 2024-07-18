import json
import sys

def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def parse(coils_out: str, version:float) -> dict:
    matches = {}
    with open(coils_out) as coils_file:
        for line in coils_file:
            line = line.split()
            seq_id = line[5]
            coil_count = line[1]
            matches[seq_id] = {"count": coil_count}
    return matches

if __name__ == '__main__':
    main()