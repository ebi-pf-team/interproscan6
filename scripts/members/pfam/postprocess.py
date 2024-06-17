import argparse
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
                    not_overlapping_locations = []
                    for candidate_location in domain_details["locations"]:
                        for filtered_location in filtered_match[2]["locations"]:
                            if _matches_overlap(candidate_location, filtered_location):
                                # if overlapping, check if they are NOT nested
                                candidate_nested = candidate_match_info.get('nested', [])
                                filtred_nested = filtred_match_info.get('nested', [])
                                nested = _matches_are_nested(candidate_nested, filtred_nested)
                                if not nested:
                                    pass
                                else:
                                    not_overlapping_locations.append(candidate_location)
                            else:
                                not_overlapping_locations.append(candidate_location)
                        domain_details["locations"] = not_overlapping_locations
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


def filter_fragments(filtered_matches):
    return filtered_matches
    # itera nos matches que ja foram filtrados
    #      pega o model id do match da vez
    #      cria uma variavel nested models pegando os nested desse id ou inicializa um set vazio
    #      IF tem algum nested relacionado a esse id E (o match ta no nested_models OU tem algum overlap:
    #         cria aquele uui e atribui ao grupo
    #         location_fragments = []
    #         nested_fragments = 0
    #         Para cada filtered_match
    #           IF o model id ta em nested_models E tem overlap:
    #                adiciona em location_fragments o location start e end do match
    #                incrementa nested_fragments
    #           ordena os locations fragments
    #           fragment_dc_status = 'CONTINUOUS'
    #           cria uma variavel raw_discontinuous_matches
    #           adiciona esse match a raw_discontinuous_matches

            # PARA cada fragment_start, fragment_end no location_fragments:
            #     inicializa new_matches_from_fragment

                # Para cada raw_discontinuous_matches:
                #     new_location_start = start desse match
                #     new_location_end = end dessa match

                    # se nao tem overlap, so adiciona no new_matches_from_fragment
                    # se tem overlap, verifica se o fragmento esta no meio do match
                    #     fragment_dc_status = 'NC_TERMINAL_DISC'
                    #     adiciona o match com o novo DC_status ao new_matches_from_fragment

                    # ta construindo uma nova location???
                    # if fragment_dc_status == 'CONTINUOUS':
                    #     fragment_dc_status = None
                    # if fragment_start <= new_location_start:
                    #     new_location_start = fragment_end + 1
                    #     fragment_dc_status = 'N_TERMINAL_DISC'
                    # elif fragment_end >= new_location_end:
                    #     new_location_end = fragment_start - 1
                    #     fragment_dc_status = fragment_dc_status ou 'C_TERMINAL_DISC', se nao tiver
                    # elif fragment_start > new_location_start and fragment_end < new_location_end:
                    #     new_location_end = fragment_start - 1
                    #     fragment_dc_status = fragment_dc_status ou 'C_TERMINAL_DISC', se nao tiver
                    # new_match = match com new_location_start, new_location_end, fragment_dc_status

                    # new_match.setLocFragmentDCStatus(fragment_dc_status)
                    # new_matches_from_fragment.append(new_match)

                    # new_location_start = fragment_end + 1

                    # SE fragment_start > new_location_start and fragment_end < new_location_end:
                    #     fragment_dc_status = 'N_TERMINAL_DISC'
                    #     final_location_end = end do match
                    #     final_match = o match com a new_location_start, final_location_end e fragment_dc_status
                    #     adiciona o final_match ao new_matches_from_fragment

                # raw_discontinuous_matches = new_matches_from_fragment

             # Para cada raw_discontinuous_matches:
             #    match_length = end - start + 1
             #    IF match_length >= o min match length:
             #        filtered_raw_protein adicionar o raw_discontinuous_match
        # Caso contrario,
    #       SE o end - start + 1 >= min_match_length
    #           adiciona o match ao filtered_raw_protein


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

    seed = args.seed.split("=")[1]
    clans = args.clans.split("=")[1]
    dat = args.dat.split("=")[1]

    clans = stockholm_parser.get_clans(seed, clans)
    dat_parsed = stockholm_parser.get_pfam_a_dat(dat)  # used to fragments block

    filtered_matches = postprocess(args.hmm_parsed, clans, dat_parsed, args.min_length)
    filtered_fragments = filter_fragments(filtered_matches)
    result = build_result(filtered_fragments)
    print(json.dumps(result, indent=4))


if __name__ == '__main__':
    main()
