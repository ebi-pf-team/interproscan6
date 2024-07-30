import json
import sys


def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1])
    print(json.dumps(parsed_results, indent=2))


def parse(phobius_out: str, version: float) -> dict:
    matches = {}
    with open(phobius_out) as ph_file:
        for line in ph_file:
            if line.startswith("ID"):
                seq_id = line.strip("ID").split(maxsplit=1)[0]

            elif line.startswith("FT"):
                line = line.strip("FT").split()

    return matches


if __name__ == '__main__':
    main()
