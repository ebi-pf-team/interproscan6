import argparse
import json
import stockholm_parser


def post_process(hmm_matches: str, model2clans: dict) -> dict:
    with open(hmm_matches, "r") as matches:
        matches_info = json.load(matches)

    sorted_matches_info = {}
    for protein_id, domains in matches_info.items():
        sorted_domains = sorted(
            domains.items(),
            key=lambda item: (float(item[1]['evalue']), -float(item[1]['score']))
        )
        sorted_matches_info[protein_id] = {domain_id: domain_details for domain_id, domain_details in sorted_domains}

    filtered_matches = {}

    for protein_id, domains in sorted_matches_info.items():
        filtered_matches[protein_id] = {}
        for domain_id, domain_details in domains.items():
            keep = True
            candidate_match_info = model2clans.get(domain_id, {})
            for filtered_domain_id, filtered_domain_details in filtered_matches[protein_id].items():
                filtered_match_info = model2clans.get(filtered_domain_id, {})
                candidate_clan = candidate_match_info.get("clan", None)
                filtered_clan = filtered_match_info.get("clan", None)
                # same clan?
                if candidate_clan and (candidate_clan == filtered_clan):
                    not_overlapping_locations = []
                    # overlap?
                    for candidate_location in domain_details["locations"]:
                        overlapping = False
                        for filtered_location in filtered_domain_details["locations"]:
                            if _matches_overlap(candidate_location, filtered_location):
                                # if overlapping, check if they are NOT nested
                                candidate_nested = candidate_match_info.get('nested', [])
                                filtered_nested = filtered_match_info.get('nested', [])
                                if not _matches_are_nested(candidate_nested, filtered_nested):
                                    overlapping = True
                                    break
                        if not overlapping:
                            not_overlapping_locations.append(candidate_location)

                    if not not_overlapping_locations:
                        keep = False
                    else:
                        domain_details["locations"] = not_overlapping_locations
            if keep:
                filtered_matches[protein_id][domain_id] = domain_details

    return filtered_matches


def build_fragments(filtered_matches: dict, dat_parsed: dict, min_length: int) -> list:
    processed_matches = []

    for protein_id, domains in filtered_matches.items():
        for domain_id, pfam_match in domains.items():
            if not pfam_match["locations"]:
                # complete sequence match but no domain matches
                continue
            model_id = pfam_match['accession']
            nested_models = dat_parsed.get(model_id, [])

            if nested_models:
                location_fragments = []
                for other_domain_id, other_match in domains.items():
                    for location in other_match["locations"]:
                        if (other_match['accession'] in nested_models) \
                                and _matches_overlap(location, pfam_match["locations"][0]):
                            location_fragments.append({
                                'start': location['start'],
                                'end': location['end']
                            })
                location_fragments.sort(key=lambda x: (x['start'], x['end']))

                fragment_dc_status = "CONTINUOUS"
                raw_discontinuous_matches = [pfam_match]
                for fragment in location_fragments:
                    two_actual_regions = False
                    (
                        new_location_start,
                        new_location_end,
                        final_location_end,
                        fragment_dc_status,
                        two_actual_regions
                    ) = _create_fragment(
                        raw_discontinuous_matches[0],
                        fragment,
                        fragment_dc_status,
                        two_actual_regions
                    )

                    for loc in raw_discontinuous_matches[0]["locations"]:  # it's always only one location here?
                        if 'location-fragments' not in loc:
                            loc['location-fragments'] = []

                        loc['location-fragments'].append({
                            'start': new_location_start,
                            'end': new_location_end,
                            'dc-status': fragment_dc_status
                        })

                        if two_actual_regions:
                            loc['location-fragments'].append({
                                'start': int(fragment['end']) + 1,
                                'end': final_location_end,
                                'dc-status': "N_TERMINAL_DISC"
                            })

                for raw_discontinuous_match in raw_discontinuous_matches:
                    match_length = int(raw_discontinuous_match["locations"][0]['end']) - \
                                   int(raw_discontinuous_match["locations"][0]['start']) + 1
                    if match_length >= min_length:
                        processed_matches.append((protein_id, domain_id, raw_discontinuous_match))
            else:
                if len(pfam_match["locations"]) > 0:
                    match_length = int(pfam_match["locations"][0]['end']) - int(pfam_match["locations"][0]['start']) + 1
                    if match_length >= min_length:
                        processed_matches.append((protein_id, domain_id, pfam_match))

    return processed_matches


def _create_fragment(raw_match: dict, fragment: dict, fragment_dc_status: str, two_actual_regions: bool):
    new_location_start = int(raw_match["locations"][0]['start'])
    new_location_end = int(raw_match["locations"][0]['end'])
    final_location_end = int(raw_match["locations"][0]['end'])

    if int(fragment['start']) <= new_location_start and int(fragment['end']) >= new_location_end:
        fragment_dc_status = "NC_TERMINAL_DISC"
    elif fragment_dc_status == "CONTINUOUS":
        fragment_dc_status = None

    if int(fragment['start']) <= new_location_start:
        new_location_start = int(fragment['end']) + 1
        fragment_dc_status = "N_TERMINAL_DISC"
    elif int(fragment['end']) >= new_location_end:
        new_location_end = int(fragment['start']) - 1
        fragment_dc_status = "C_TERMINAL_DISC"
    elif int(fragment['start']) > new_location_start and int(fragment['end']) < new_location_end:
        new_location_end = int(fragment['start']) - 1
        two_actual_regions = True
        fragment_dc_status = "C_TERMINAL_DISC"

    return new_location_start, new_location_end, final_location_end, fragment_dc_status, two_actual_regions


def _matches_overlap(one: dict, two: dict) -> bool:
    return max(int(one['start']), int(two['start'])) <= min(int(one['end']), int(two['end']))


def _matches_are_nested(one: dict, two: dict):
    if one is None or two is None:
        return False
    return not set(one).isdisjoint(set(two))


def build_result(filtered_matches):
    result = {}
    for protein_id, domain_id, domain_details in filtered_matches:
        if protein_id not in result:
            result[protein_id] = {}
        result[protein_id][domain_id] = domain_details
    return result


def main():
    parser = argparse.ArgumentParser(description='Postprocess pfam matches')
    parser.add_argument('--hmm_parsed', required=True, help='hmm parsed file')
    parser.add_argument('--min_length', required=True, help='minimum length')
    parser.add_argument('--seed', required=True, help='seed file')
    parser.add_argument('--clans', required=True, help='clans file')
    parser.add_argument('--dat', required=True, help='dat file')
    args = parser.parse_args()

    hmm_parsed = args.hmm_parsed
    seed = args.seed.split("=")[1]
    clans = args.clans.split("=")[1]
    dat = args.dat.split("=")[1]
    min_length = int(args.min_length)

    # Need to return Pfam clans AND nesting relationships between models
    seed_nesting = stockholm_parser.parser_seed_nesting(seed)
    clans_parsed = stockholm_parser.parser_clans(clans, seed_nesting)
    dat_parsed = stockholm_parser.get_pfam_a_dat(dat)

    filtered_matches = post_process(hmm_parsed, clans_parsed)
    filtered_fragments = build_fragments(filtered_matches, dat_parsed, min_length)
    result = build_result(filtered_fragments)
    print(json.dumps(result, indent=4))


if __name__ == '__main__':
    main()
