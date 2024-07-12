import json
import re
import sys

MATCH_REGIONS_PATTERN = re.compile(r"^\d+-\d+(,\d+-\d+)*$")


def parse(ass3_out_path: str, modelac2acc: dict) -> dict:
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
                acc_id = f"SSF{modelac2acc[model_ac]}"

                if sequence_id not in data:
                    data[sequence_id] = {}
                if acc_id not in data[sequence_id]:
                    data[sequence_id][acc_id] = {
                        'accession': acc_id,
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
                    dc_statuses = ["C_TERMINAL_DISC"] + ["INTERNAL"] * (len(fragments) - 2) + ["N_TERMINAL_DISC"]
                    fragments_with_status = [
                        {**frag, 'dc-status': dc_status}
                        for frag, dc_status in zip(fragments, dc_statuses)
                    ]
                    new_location = {
                        'start': fragments[0]['start'],
                        'end': fragments[-1]['end'],
                        'location-fragments': fragments_with_status,
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


def modelac2acc(model_tab_path: str):
    modelac2acc = {}
    with open(model_tab_path, 'r') as reader:
        for line in reader:
            values = line.split()
            if len(values) < 2:
                continue
            modelac2acc[values[0]] = values[1]
    return modelac2acc


def main():
    args = sys.argv[1:]
    modelac2acc_dict = modelac2acc(args[0])
    parsed_result = parse(args[1], modelac2acc_dict)
    print(json.dumps(parsed_result, indent=4))


if __name__ == "__main__":
    main()
