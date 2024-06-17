import json
import os
import re

ACCESSION_LINE = "#=GF AC"
CLAN_MEMBER_MODEL_AC_LINE = "#=GF MB"
NESTING_LINE = "#=GF NE"
ACCESSION_EXTRACTOR_PATTERN = re.compile("^\\#=GF\\s+[A-Z]{2}\\s+([A-Z0-9]+).*$")


# it needs to return all the details of Pfam clans AND nesting relationships between models.
def get_clans(pfam_a_seed_file: str, pfam_clans_file: str):
    seed_nesting = parser_seed_nesting(pfam_a_seed_file)
    clans_parsed = parser_clans(pfam_clans_file, seed_nesting)
    return clans_parsed


def parser_seed_nesting(pfam_a_seed_file: str) -> dict[str, str]:
    nesting_info = {}
    with open(pfam_a_seed_file, 'r') as file:
        for line in file:
            if line.startswith(ACCESSION_LINE):
                model_accession = re.match(ACCESSION_EXTRACTOR_PATTERN, line)
                if model_accession:
                    accession = model_accession.group(1)
                else:
                    raise ValueError(f"Cannot parser the model accession: {line}")
            elif line.startswith(NESTING_LINE):
                nested_accession = re.match(ACCESSION_EXTRACTOR_PATTERN, line)
                if nested_accession:
                    try:
                        nesting_info[accession]["nested"].append(nested_accession.group(1))
                    except KeyError:
                        nesting_info[accession] = {"nested": [nested_accession.group(1)]}
    return nesting_info


def parser_clans(pfam_clans_file: str, parsed_seed: dict[str, str]) -> dict[str, str]:
    with open(pfam_clans_file, 'r', errors="replace") as file:
        for line in file:
            if line.startswith(ACCESSION_LINE):
                clan_accession = re.match(ACCESSION_EXTRACTOR_PATTERN, line)
                if clan_accession:
                    accession = clan_accession.group(1)
                else:
                    raise ValueError(f"Cannot parser the clan accession: {line}")
            elif line.startswith(CLAN_MEMBER_MODEL_AC_LINE):
                model_accession = re.match(ACCESSION_EXTRACTOR_PATTERN, line)
                if model_accession:
                    try:
                        parsed_seed[model_accession.group(1)]["clan"] = accession
                    except KeyError:
                        parsed_seed[model_accession.group(1)] = {"clan": accession}
    return parsed_seed


def parse_seed_alignment(file_path):
    seed_alignments = {}
    model_ac = None
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("//"):
                model_ac = None
            elif line.startswith("#=GF AC"):
                match = re.match(r'^#=GF\s+AC\s+([A-Z0-9]+)\..*$', line)
                if match:
                    model_ac = match.group(1)
            elif line.startswith("#=GS"):
                match = re.match(r'^#=GS\s+.+?/(\d+)-(\d+)\s+AC\s+([A-Z0-9]+)\.(\d+)$', line)
                if match and model_ac:
                    start_coordinate = int(match.group(1))
                    stop_coordinate = int(match.group(2))
                    uniprot_ac = match.group(3)
                    version_number = int(match.group(4))
                    md5 = "" #retrieve_md5_from_uniprot(uniprot_ac, version_number)
                    if uniprot_ac not in seed_alignments:
                        seed_alignments[uniprot_ac] = []
                    seed_alignment = (model_ac, md5, start_coordinate, stop_coordinate)
                    seed_alignments[uniprot_ac].append(seed_alignment)

    return seed_alignments


def get_pfam_a_dat(pfam_a_dat_file):
    alt_pfam_hmm_data = {}
    domain_name_to_accession = {}
    pfam_hmm_data = {}
    accession = None
    name = None
    clan = None
    nested_domains = set()

    with open(pfam_a_dat_file, 'r') as reader:
        for line_number, line in enumerate(reader, 1):
            if line.startswith("#=GF "):
                parts = line.split(maxsplit=2)
                if len(parts) >= 3:
                    tag = parts[1]
                    value = parts[2]
                    if tag == "ID":
                        if name is None:
                            name = value
                    elif tag == "AC":
                        if accession is None:
                            accession = value.split(".")[0]
                    elif tag == "NE":
                        domain_name = value
                        nested_domains.add(domain_name)
                    elif tag == "CL":
                        if clan is None:
                            clan = value
                elif line.startswith("//"):
                    if accession is not None:
                        domain_name_to_accession[name] = accession
                        pfam_hmm_data[accession] = nested_domains
                        accession = None
                        name = None
                        clan = None
                        nested_domains = set()

    for accession, nested_domains in pfam_hmm_data.items():
        domain_accessions = {domain_name_to_accession[domain_name] for domain_name in nested_domains}
        alt_pfam_hmm_data[accession] = domain_accessions

    return alt_pfam_hmm_data


# if __name__ == '__main__':
    # alignments = parse_seed_alignment("/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.seed")
    # 'W5KKI6': [('PF03813', '', 169, 309), ('PF17403', '', 312, 452), ('PF17404', '', 455, 617), ('PF17405', '', 629, 835), ('PF17406', '', 837, 992), ('PF17407', '', 994, 1128)], 'F1Q4S6': [('PF03813', '', 183, 325), ('PF17403', '', 328, 468), ('PF17404', '', 471, 633), ('PF17405', '', 645, 851), ('PF17406', '', 853, 1008), ('PF17407', '', 1010, 1144)],
    # clans = get_clans("/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.seed", "/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_clans")
    # {'PF00001': {'nested': ['PF01498'], 'clan': 'CL0192'}, 'PF00005': {'nested': ['PF00385'], 'clan': 'CL0023'}, 'PF00006': {'nested': ['PF05203', 'PF05204'], 'clan': 'CL0023'}, 'PF00026': {'nested': ['PF03489', 'PF05184'], 'clan': 'CL0129'}, 'PF00060': {'nested': ['PF00497'], 'clan': 'CL0030'}, 'PF00082': {'nested': ['PF02225', 'PF07676']}, 'PF00133': {'nested': ['PF08264'], 'clan': 'CL0039'}, 'PF00152': {'nested': ['PF02938', 'PF00454'], 'clan': 'CL0040'}, 'PF00176': {'nested': ['PF02037', 'PF00628', 'PF00538', 'PF07496', 'PF00439'], 'clan': 'CL0023'}, 'PF00270': {'nested': ['PF00622'], 'clan': 'CL0023'}, 'PF00278': {'nested': ['PF02784'], 'clan': 'CL0105'},
    # aninhamentos = parser_seed_nesting("/Users/lcf/PycharmProjects/interproscan6/data/pfam/37.0/pfam_a.seed")
    # print(aninhamentos)
