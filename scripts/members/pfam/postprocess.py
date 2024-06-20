import argparse
import json
import stockholm_parser


def post_process(hmm_matches, model2clans):
    with open(hmm_matches, "r") as matches:
        matches_info = json.load(matches)

    matches_list = []
    for protein_id, domains in matches_info.items():
        for domain_id, domain_details in domains.items():
            matches_list.append((protein_id, domain_id, domain_details))
    # evalue ASC e score DESC to keep the best matches in the optimistic algorithm below
    matches_sorted = sorted(matches_list, key=lambda match: (float(match[2]['evalue']), -float(match[2]['score'])))

    filtered_matches = []
    for candidate_match in matches_sorted:
        protein_id, domain_id, domain_details = candidate_match
        keep = True
        try:
            candidate_match_info = model2clans[domain_id]
            for filtered_match in filtered_matches:
                filtred_match_info = model2clans.get(filtered_match[1], {})
                # same clan?
                if candidate_match_info["clan"] == filtred_match_info["clan"]:
                    # overlap?
                    not_overlapping_locations = []
                    for candidate_location in domain_details["locations"]:
                        overlapping = False
                        for filtered_location in filtered_match[2]["locations"]:
                            if _matches_overlap(candidate_location, filtered_location):
                                # if overlapping, check if they are NOT nested
                                candidate_nested = candidate_match_info.get('nested', [])
                                filtered_nested = filtred_match_info.get('nested', [])
                                if not _matches_are_nested(candidate_nested, filtered_nested):
                                    overlapping = True
                                    break
                        if not overlapping:
                            not_overlapping_locations.append(candidate_location)

                    if not not_overlapping_locations:
                        keep = False
                        break
                    else:
                        domain_details["locations"] = not_overlapping_locations
        except KeyError:
            pass

        if keep:
            filtered_matches.append(candidate_match)

    return filtered_matches


def build_fragments(filtered_matches, dat_parsed, min_length):
    processed_matches = []
    for info_pfam in filtered_matches:
        protein_id, domain_id, pfam_match = info_pfam
        model_id = pfam_match['accession']
        nested_models = dat_parsed.get(model_id, [])

        if nested_models:
            location_fragments = []
            for info_raw in filtered_matches:
                raw_protein_id, raw_domain_id, raw_match = info_raw
                for location in raw_match["locations"]:
                    if (raw_match['accession'] in nested_models) and _matches_overlap(location, pfam_match["locations"][0]):
                        location_fragments.append({
                            'start': location['start'],
                            'end': location['end']
                        })
            location_fragments.sort(key=lambda x: (x['start'], x['end']))

            fragment_dc_status = "CONTINUOUS"
            raw_discontinuous_matches = [pfam_match]
            for fragment in location_fragments:
                two_actual_regions = False
                new_location_start, new_location_end, final_location_end, fragment_dc_status, two_actual_regions = _create_fragment(
                    raw_discontinuous_matches[0], fragment, fragment_dc_status, two_actual_regions)

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
                match_length = int(raw_discontinuous_match["locations"][0]['end']) - int(raw_discontinuous_match["locations"][0]['start']) + 1
                if match_length >= min_length:
                    processed_matches.append((protein_id, domain_id, raw_discontinuous_match))
        else:
            match_length = int(pfam_match["locations"][0]['end']) - int(pfam_match["locations"][0]['start']) + 1
            if match_length >= min_length:
                processed_matches.append((protein_id, domain_id, pfam_match))

    return processed_matches


def _create_fragment(raw_match, fragment, fragment_dc_status, two_actual_regions):
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


def _matches_overlap(one, two) -> bool:
    return max(one['start'], two['start']) <= min(one['end'], two['end'])


def _matches_are_nested(one, two):
    if one is None or two is None:
        return False
    return not set(one).isdisjoint(set(two))


def _regions_overlap(start_region_one, end_region_one, start_region_two, end_region_two):
    return max(start_region_one, start_region_two) <= min(end_region_one, end_region_two)


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

    clans = stockholm_parser.get_clans(seed, clans)
    dat_parsed = stockholm_parser.get_pfam_a_dat(dat)

    filtered_matches = post_process(hmm_parsed, clans)
    filtered_fragments = build_fragments(filtered_matches, dat_parsed, min_length)
    result = build_result(filtered_fragments)
    print(json.dumps(result, indent=4))


if __name__ == '__main__':
    main()
