import json
import re
import sys
from collections import namedtuple

MATCH_REGIONS_PATTERN = re.compile(r"^\d+-\d+(,\d+-\d+)*$")
NAME_LINE = re.compile(r"^NAME\s+(.+)$")
DESC_LINE = re.compile(r"^DESC\s+(.+)$")
ACCESSION_PATTERN = re.compile(r"^ACC\s+([A-Z0-9]+)\.?.*$")
LENGTH_LINE = re.compile(r"^LENG\s+([0-9]+)$")


def parse(ass3_out_path: str, hmmlib_info: dict, member_db: str) -> dict:
    data = {}

    with open(ass3_out_path, 'r') as reader:
        for line in reader:
            match = _parse_line(line.strip())
            if match:
                acc_id = hmmlib_info[match['model_id']].acc_id
                name = hmmlib_info[match['model_id']].description
                hmm_length = hmmlib_info[match['model_id']].length

                if match['sequence_id'] not in data:
                    data[match['sequence_id']] = {}
                if acc_id not in data[match['sequence_id']]:
                    data[match['sequence_id']][acc_id] = {
                        'accession': acc_id,
                        'name': name,
                        'hmm_length': hmm_length,
                        'member_db': member_db,
                        'evalue': match['evalue'],
                        'model-ac': match['model_id'],
                        'locations': []
                    }

                match_regions = match['match_regions'].split(',')
                fragments = []
                for region in match_regions:
                    start, end = map(int, region.split('-'))
                    fragments.append({'start': start, 'end': end})

                if len(fragments) > 1:
                    fragments.sort(key=lambda x: x['start'])
                    fragments[0]['dc-status'] = "C_TERMINAL_DISC"
                    fragments[-1]['dc-status'] = "N_TERMINAL_DISC"
                    for i in range(1, len(fragments) - 1):
                        fragments[i]['dc-status'] = "NC_TERMINAL_DISC"
                    new_location = {
                        'start': fragments[0]['start'],
                        'end': fragments[-1]['end'],
                        'location-fragments': fragments,
                        'representative': '',
                        'evalue': match['evalue'],
                        'alignment': match['alignment_to_model'],
                        'family_evalue': match['family_evalue'],
                        'scop_domain_id': match['scop_domain_id'],
                        'scop_family_id': match['scop_family_id']
                    }
                    data[match['sequence_id']][acc_id]['locations'].append(new_location)
                else:
                    location = {
                        'start': fragments[0]['start'],
                        'end': fragments[0]['end'],
                        'representative': '',
                        'evalue': match['evalue'],
                        'alignment': match['alignment_to_model'],
                        'family_evalue': match['family_evalue'],
                        'scop_domain_id': match['scop_domain_id'],
                        'scop_family_id': match['scop_family_id']
                    }
                    data[match['sequence_id']][acc_id]['locations'].append(location)
    return data


def _parse_line(line: str) -> list:
    if not line or line.startswith("-"):
        return None

    values = line.split()
    if len(values) != 9 or not MATCH_REGIONS_PATTERN.match(values[2]):
        return None
    try:
        match = {
            'sequence_id': values[0],
            'model_id': values[1],
            'match_regions': values[2],
            'evalue': float(values[3]),
            'model_match_start_pos': int(values[4]),
            'alignment_to_model': values[5],
            'family_evalue': float(values[6]),
            'scop_domain_id': int(values[7]),
            'scop_family_id': int(values[8])
        }
    except ValueError:
        return None

    return match


def parse_hmmlib(hmmlib_path: str) -> dict:
    model_info = {}
    with open(hmmlib_path, 'r') as file:
        ModelInfo = namedtuple('ModelInfo', ['acc_id', 'description', 'length'])
        acc_id, model_ac, description, length = None, None, None, None
        for line in file:
            line = line.strip()
            if line and line.startswith('//'):
                if model_ac:
                    model_info[model_ac] = ModelInfo(acc_id, description, length)
                acc_id, model_ac, description, length = None, None, None, None
            elif line:
                if line.startswith('A') and not acc_id:
                    match = ACCESSION_PATTERN.match(line)
                    if match:
                        acc_id = "SSF" + match.group(1)
                elif line.startswith('D') and not description:
                    match = DESC_LINE.match(line)
                    if match:
                        description = match.group(1)
                elif line.startswith('N') and not model_ac:
                    match = NAME_LINE.match(line)
                    if match:
                        model_ac = match.group(1)
                elif line.startswith('L') and not length:
                    match = LENGTH_LINE.match(line)
                    if match:
                        length = int(match.group(1))

        if model_ac:
            model_info[model_ac] = ModelInfo(acc_id, description, length)

    return model_info


def main():
    """CL input:
    0. Str repr of path to the hmm lib
    1. Str repr of path to the ass3 output file
    2. member database
    3. Str repr of path to write the output file
    """
    args = sys.argv[1:]
    hmmlib_info = parse_hmmlib(args[0])
    parsed_result = parse(args[1], hmmlib_info, args[1])
    with open(args[3], "w") as fh:
        json.dump(parsed_result, fh)


if __name__ == "__main__":
    main()
