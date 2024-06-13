import json
import os
import re

ACCESSION_LINE = "#=GF AC"
CLAN_MEMBER_MODEL_AC_LINE = "#=GF MB"
NESTING_LINE = "#=GF NE"
ACCESSION_EXTRACTOR_PATTERN = re.compile("^\\#=GF\\s+[A-Z]{2}\\s+([A-Z0-9]+).*$")


# it needs to return all the details of Pfam clans AND nesting relationships between models.
def get_clan(pfam_a_seed_file: str, pfam_clans_file: str):
    seed_parsed = parser_seed(pfam_a_seed_file)
    clan_parsed = parser_clans(pfam_clans_file, seed_parsed)
    return clan_parsed


def parser_seed(pfam_a_seed_file: str) -> dict[str, str]:
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



def get_pfam_a_dat(pfam_a_dat_file):
    alt_pfam_hmm_data = {}
    domain_name_to_accession = {}
    pfam_hmm_data = {}
    accession = None
    name = None
    clan = None
    nested_domains = set()

    with open(pfam_hmm_data_path, 'r') as reader:
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
