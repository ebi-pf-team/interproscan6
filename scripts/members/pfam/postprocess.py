import json
import stockholm_parser


#  return a Map of proteinIds to a List of filtered matches.
def postprocess(hmm_matches, model2clans, dat_parsed, min_length):
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
                filtred_match_info = model2clans[filtered_match[1]]
                # same clan?
                if candidate_match_info["clan"] == filtred_match_info["clan"]:
                    # overlap?
                    not_overlapping_locations = set()
                    for candidate_location in candidate_match[2]["locations"]:
                        for filtered_location in filtered_match[2]["locations"]:
                            if _matches_overlap(candidate_location, filtered_location):
                                # if overlapping, check if they are NOT nested
                                candidate_nested = candidate_match_info['nested'] if 'nested' in candidate_match_info else []
                                filtred_nested = filtred_match_info['nested'] if 'nested' in filtred_match_info else []
                                nested = _matches_are_nested(candidate_nested, filtred_nested)
                                if not nested:
                                    pass
                                else:
                                    not_overlapping_locations.add(candidate_location)
                            else:
                                not_overlapping_locations.add(candidate_location)
                    if len(not_overlapping_locations) != 0:
                        filtered_match[2]["locations"] = list(not_overlapping_locations)
                    else:
                        keep = False
                        break

        except KeyError:
            pass

        if keep:
            filtered_matches.append(candidate_match)

    return filtered_matches


def _matches_overlap(one, two) -> bool:
    return max(one['start'], two['start']) <= min(one['end'], two['end'])


def _matches_are_nested(one, two):
    if one is None or two is None:
        return False
    return not set(one).isdisjoint(set(two))


def _regions_overlap(start_region_one, end_region_one, start_region_two, end_region_two):
    return max(start_region_one, start_region_two) <= min(end_region_one, end_region_two)


def main():
    # args = sys.argv[1:]
    #
    # hmm_parsed = args[0]
    # min_length = args[1]
    # pfam_a_seed_file = args[2]
    # pfam_clans_file = args[3]
    # pfam_a_dat_file = args[4]

    clans = stockholm_parser.get_clan("/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.seed",
                                                     "/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_clans")

    dat_parsed = "" # stockholm_parser.parse_pfam_a_dat("/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.dat")

    hmmer_out_parsed = "/Users/lcf/PycharmProjects/interproscan6/work/de/5d75df573812164ace5cc776c2d06c/hmmer_parsed_37.0_pfam_a.hmm.out.json"
    min_length = 8
    pp_pfam_result = postprocess(hmmer_out_parsed, clans, dat_parsed, min_length)
    print(json.dumps(pp_pfam_result, indent=4))



if __name__ == '__main__':
    main()
