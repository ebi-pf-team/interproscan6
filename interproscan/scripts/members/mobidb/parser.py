import args
import re

END_OF_RECORD_MARKER = "//"
PROTEIN_ID_LINE_START = '>'
DOMAIN_LINE_PATTERN = re.compile(r"^(\S+)\s+(\d+)\s+(\d+).*$")


def parse(input_file, library, release):
    match_data = {}
    raw_matches = parse_file_input(input_file, library, release)

    for raw_match in raw_matches:
        sequence_id = raw_match["sequence_identifier"]
        if sequence_id in match_data:
            match_data[sequence_id].append(raw_match)
        else:
            match_data[sequence_id] = [raw_match]

    return match_data

def parse_file_input(input_file, library, release):
    matches = []
    with open(input_file, 'r') as reader:
        for line in reader:
            line = line.strip()
            if line.startswith(PROTEIN_ID_LINE_START) or line == END_OF_RECORD_MARKER:
                continue

            match = DOMAIN_LINE_PATTERN.match(line)
            if match:
                sequence_identifier = match.group(1)
                location_start = int(match.group(2))
                location_end = int(match.group(3))
                feature = match.group(4).strip() if match.group(4) else ""
                matches.append({
                    "sequence_identifier": sequence_identifier,
                    "member": library,
                    "release": release,
                    "location_start": location_start,
                    "location_end": location_end,
                    "feature": feature
                })
    return matches


def main():
    args = sys.argv[1:]
    matches = parse(args[0], args[1], args[2])

    print(json.dumps(matches, indent=4))


if __name__ == "__main__":
    main()
