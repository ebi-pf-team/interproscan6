import json
import os
import re

ACCESSION_LINE = "#=GF AC"
CLAN_MEMBER_MODEL_AC_LINE = "#=GF MB"
NESTING_LINE = "#=GF NE"
ACCESSION_EXTRACTOR_PATTERN = re.compile("^\\#=GF\\s+[A-Z]{2}\\s+([A-Z0-9]+).*$")


def get_pfamA_seed(pfam_seed_file):
    # Maybe We shouldnt parse the seed file as it is not being used at the moment (???)
    models = {}
    record = {"model_ac": None, "nested_domains": []}

    with open(pfam_seed_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("//"):
                if record["model_ac"]:
                    models[record["model_ac"]] = record["nested_domains"]
                record = {"model_ac": None, "nested_domains": []}
            elif line.startswith(ACCESSION_LINE):
                match = ACCESSION_EXTRACTOR_PATTERN.match(line)
                if match:
                    record["model_ac"] = match.group(1)
            elif line.startswith(NESTING_LINE):
                match = ACCESSION_EXTRACTOR_PATTERN.match(line)
                if match:
                    record["nested_domains"].append(match.group(1))
                else:
                    raise ValueError(
                        f"Line: {line} appears to be a nesting line, but does not match the regex to retrieve the accession number.")

    if record["model_ac"] or record["nested_domains"]:
        if record["model_ac"]:
            models[record["model_ac"]] = record["nested_domains"]

    return models


def get_pfam_clans(pfam_clan_file):
    clans = {}
    current_clan = None

    with open(pfam_clan_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("//"):
                current_clan = None
            elif line.startswith(ACCESSION_LINE):
                match = ACCESSION_EXTRACTOR_PATTERN.match(line)
                if match:
                    accession = match.group(1)
                    current_clan = accession
                    if accession not in clans:
                        clans[accession] = {"members": []}
            elif line.startswith(CLAN_MEMBER_MODEL_AC_LINE):
                if current_clan is None:
                    raise ValueError(f"Found a clan member line without a current clan: {line}")
                match = ACCESSION_EXTRACTOR_PATTERN.match(line)
                if match:
                    member_accession = match.group(1)
                    clans[current_clan]["members"].append(member_accession)
                else:
                    raise ValueError(f"Cannot parse the clan member accession: {line}")

        return clans


def build_model(record, model_accession_nests_model_accession):
    if record["model_ac"]:
        if record["model_ac"] not in model_accession_nests_model_accession:
            model_accession_nests_model_accession[record["model_ac"]] = []
        model_accession_nests_model_accession[record["model_ac"]].extend(record["nested_domains"])


def add_nesting_information(model_accession_nests_model_accession, clan_data):
    for nesting_model_ac, nested_model_acs in model_accession_nests_model_accession.items():
        nesting_model = clan_data.get(nesting_model_ac)
        if nesting_model is None:
            raise ValueError(f"Cannot find nesting model with accession {nesting_model_ac}")

        for nested_model_ac in nested_model_acs:
            nested_model = clan_data.get(nested_model_ac)
            if nested_model is None:
                raise ValueError(f"Cannot find nested model with accession {nested_model_ac}")
            if "nested_in" not in nested_model:
                nested_model["nested_in"] = []
            nested_model["nested_in"].append(nesting_model_ac)


def get_pfam_a_dat(pfam_a_dat_file):
    domain_name_to_accession = {}
    pfam_hmm_data = {}

    with open(pfam_a_dat_file, 'r') as file:
        accession = None
        nested_domains = set()

        for line in file:
            line = line.strip()

            if line.startswith("#=GF ID"):
                accession = line.split()[-1]
            elif line.startswith("#=GF NE"):
                domain_name = line.split()[-1]
                nested_domains.add(domain_name)
            elif line.startswith("//"):
                if accession:
                    domain_name_to_accession[accession] = accession
                    pfam_hmm_data[accession] = nested_domains
                accession = None
                nested_domains = set()

    return pfam_hmm_data


def main():
    # args = sys.argv[1:]
    #
    # # parsePfamASeed(); //think of an easier way to get this info as the rest of the alignment info is not used (???)
    # seed_parsed = get_pfamA_seed(args[0])
    # clan_parsed = get_pfam_clans(args[1])
    # dat_parsed = get_pfamA_dat(args[2])
    # add_nesting_information(model_accession_nests_model_accession, clan_parsed)

    # parse_result = get_pfamA_dat("/Users/lcf/PycharmProjects/interproscan6/scripts/members/pfam/test_dat")

    # it needs to return all the details of Pfam clans AND nesting relationships between models.
    # print(json.dumps(parse_result, indent=2))


    pfam_hmm_data_path = '/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.dat'
    alt_pfam_hmm_data = get_pfam_a_dat(pfam_hmm_data_path)

    for key, value in alt_pfam_hmm_data.items():
        if len(value) > 0:
            print(key, value)


if __name__ == '__main__':
    main()

