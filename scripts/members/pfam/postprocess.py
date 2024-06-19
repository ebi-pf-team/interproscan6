import argparse
import json
import stockholm_parser


def postprocess(hmm_matches, model2clans):
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


def build_fragments(filtered_matches, dat_parsed, min_length):
    filtered_raw_protein = []

    for match in filtered_matches:
        match_seq_id, match_acc, match_info = match
        try:
            nested_models = dat_parsed[match_acc]
        except KeyError:
            nested_models = []

        for location in match_info["locations"]:
            if nested_models:
                location_fragments = []
                for raw_match in filtered_matches:
                    raw_seq_id, raw_acc, raw_info = raw_match
                    for raw_location in raw_info["locations"]:
                        if raw_acc in nested_models and _matches_overlap(raw_location, location):
                            location_fragments.append((int(raw_location["start"]), int(raw_location["end"])))

                location_fragments.sort()

                fragment_dc_status = "CONTINUOUS"
                raw_discontinuous_matches = [location]
                location["location-fragments"] = []

                for fragment_start, fragment_end in location_fragments:
                    new_matches_from_fragment = []

                    for raw_discontinuous_match in raw_discontinuous_matches:
                        new_location_start = int(raw_discontinuous_match["start"])
                        new_location_end = int(raw_discontinuous_match["end"])

                        if not _regions_overlap(new_location_start, new_location_end, fragment_start, fragment_end):
                            new_matches_from_fragment.append(raw_discontinuous_match)
                            continue

                        if fragment_start <= new_location_start and fragment_end >= new_location_end:
                            fragment_dc_status = "NC_TERMINAL_DISC"
                            raw_discontinuous_match["loc_fragment_dc_status"] = fragment_dc_status
                            new_matches_from_fragment.append(raw_discontinuous_match)
                            continue

                        if fragment_dc_status == "CONTINUOUS":
                            fragment_dc_status = None
                        if fragment_start <= new_location_start:
                            new_location_start = fragment_end + 1
                            fragment_dc_status = "N_TERMINAL_DISC"
                        elif fragment_end >= new_location_end:
                            new_location_end = fragment_start - 1
                            fragment_dc_status = "C_TERMINAL_DISC"
                        elif fragment_start > new_location_start and fragment_end < new_location_end:
                            new_location_end = fragment_start - 1
                            fragment_dc_status = "C_TERMINAL_DISC"

                        new_match = {
                            "start": new_location_start,
                            "end": new_location_end,
                            "dc-status": fragment_dc_status
                        }
                        new_matches_from_fragment.append(new_match)
                        location["location-fragments"].append(new_match)

                    raw_discontinuous_matches = new_matches_from_fragment

                for raw_discontinuous_match in raw_discontinuous_matches:
                    match_length = int(raw_discontinuous_match["end"]) - int(raw_discontinuous_match["start"]) + 1
                    if match_length >= min_length:
                        filtered_raw_protein.append({
                            "model_id": match_acc,
                            "location_start": raw_discontinuous_match["start"],
                            "location_end": raw_discontinuous_match["end"],
                            "loc_fragment_dc_status": raw_discontinuous_match.get("loc_fragment_dc_status", "")
                        })
            else:
                match_length = int(location["end"]) - int(location["start"]) + 1
                if match_length >= min_length:
                    filtered_raw_protein.append({
                        "model_id": match_acc,
                        "location_start": int(location["start"]),
                        "location_end": int(location["end"]),
                        "loc_fragment_dc_status": location.get("loc_fragment_dc_status", "")
                    })

    return filtered_raw_protein


def _matches_overlap(one, two) -> bool:
    return max(one['start'], two['start']) <= min(one['end'], two['end'])


def _matches_are_nested(one, two):
    if one is None or two is None:
        return False
    return not set(one).isdisjoint(set(two))


def _regions_overlap(start_region_one, end_region_one, start_region_two, end_region_two):
    return max(start_region_one, start_region_two) <= min(end_region_one, end_region_two)


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

                    # constroi a nova location com o dc-status

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
    # parser = argparse.ArgumentParser(description='Postprocess pfam matches')
    # parser.add_argument('--hmm_parsed', required=True, help='hmm parsed file')
    # parser.add_argument('--min_length', required=True, help='minimum length')
    # parser.add_argument('--seed', required=True, help='seed file')
    # parser.add_argument('--clans', required=True, help='clans file')
    # parser.add_argument('--dat', required=True, help='dat file')
    # args = parser.parse_args()
    #
    # hmm_parsed = args.hmm_parsed
    # seed = args.seed.split("=")[1]
    # clans = args.clans.split("=")[1]
    # dat = args.dat.split("=")[1]
    # min_length = int(args.min_length)

    hmm_parsed = "/Users/lcf/PycharmProjects/interproscan6/work/ad/0ee28305825498f95f502350293ff0/hmmer_parsed_37.0_pfam_a.hmm.out.json"
    seed = '/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.seed'
    clans = '/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_clans'
    dat = '/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.dat'
    min_length = 8

    clans = stockholm_parser.get_clans(seed, clans)
    dat_parsed = stockholm_parser.get_pfam_a_dat(dat)

    filtered_matches = postprocess(hmm_parsed, clans)
    # print(filtered_matches)
    filtered_fragments = build_fragments(filtered_matches, dat_parsed, min_length)
    # print(json.dumps(filtered_fragments, indent=4))
    # result = build_result(filtered_fragments)
    # print(json.dumps(result, indent=4))


if __name__ == '__main__':
    main()
