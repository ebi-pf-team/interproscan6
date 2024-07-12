import json
import re
import sys

MATCH_REGIONS_PATTERN = re.compile(r"^\d+-\d+(,\d+-\d+)*$")


def parse(ass3_out_path: str, hmmlib_info: dict) -> dict:
    data = {}
    path_segments = ass3_out_path.split("/")[-1].split("._.")
    version = path_segments[0]
    member_db = path_segments[1]
    with open(ass3_out_path, 'r') as reader:
        for line in reader:
            raw_matches = _parse_line(line.strip())
            for match in raw_matches:
                sequence_id = match['sequence_id']
                model_ac = match['model_id']
                acc_id = hmmlib_info[model_ac][0]
                name = hmmlib_info[model_ac][1]
                hmm_length = hmmlib_info[model_ac][2]

                if sequence_id not in data:
                    data[sequence_id] = {}
                if acc_id not in data[sequence_id]:
                    data[sequence_id][acc_id] = {
                        'accession': acc_id,
                        'name': name,
                        'hmm_length': hmm_length,
                        'member_db': member_db,
                        'version': version,
                        'evalue': match['evalue'],
                        'model-ac': model_ac,
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
                    data[sequence_id][acc_id]['locations'].append(new_location)
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
                    data[sequence_id][acc_id]['locations'].append(location)
    return data


def _parse_line(line: str) -> list:
    if not line or line.startswith("-"):
        return []
    values = line.split()
    if len(values) != 9:
        return []
    try:
        sequence_id = values[0]
        model_id = values[1]
        match_regions = values[2]
        evalue = float(values[3])
        model_match_start_pos = int(values[4])
        alignment_to_model = values[5]
        family_evalue = float(values[6])
        scop_domain_id = int(values[7])
        scop_family_id = int(values[8])
    except ValueError:
        return []

    matches = []
    if MATCH_REGIONS_PATTERN.match(match_regions):
        match = {
            'sequence_id': sequence_id,
            'model_id': model_id,
            'match_regions': match_regions,
            'evalue': evalue,
            'model_match_start_pos': model_match_start_pos,
            'alignment_to_model': alignment_to_model,
            'family_evalue': family_evalue,
            'scop_domain_id': scop_domain_id,
            'scop_family_id': scop_family_id
        }
        matches.append(match)
    return matches


def parse_hmmlib(hmmlib_path: str) -> dict:
    NAME_LINE = re.compile(r"^NAME\s+(.+)$")
    DESC_LINE = re.compile(r"^DESC\s+(.+)$")
    ACCESSION_PATTERN = re.compile(r"^ACC\s+([A-Z0-9]+)\.?.*$")
    LENGTH_LINE = re.compile(r"^LENG\s+([0-9]+)$")

    model_info = {}

    with open(hmmlib_path, 'r') as file:
        acc_id = None
        model_ac = None
        description = None
        length = None

        for line in file:
            line = line.strip()
            if line and line.startswith('//'):
                if model_ac is not None:
                    model_info[model_ac] = (acc_id, description, length)
                acc_id = None
                model_ac = None
                description = None
                length = None
            elif line:
                if line.startswith('A') and acc_id is None:
                    match = ACCESSION_PATTERN.match(line)
                    if match:
                        acc_id = "SSF" + match.group(1)
                elif line.startswith('D') and description is None:
                    match = DESC_LINE.match(line)
                    if match:
                        description = match.group(1)
                elif line.startswith('N') and model_ac is None:
                    match = NAME_LINE.match(line)
                    if match:
                        model_ac = match.group(1)
                elif line.startswith('L') and length is None:
                    match = LENGTH_LINE.match(line)
                    if match:
                        length = int(match.group(1))

        if model_ac is not None:
            model_info[model_ac] = (acc_id, description, length)

    return model_info


def main():
    args = sys.argv[1:]
    hmmlib_info = parse_hmmlib(args[0])
    parsed_result = parse(args[1], hmmlib_info)
    print(json.dumps(parsed_result, indent=4))


if __name__ == "__main__":
    main()
