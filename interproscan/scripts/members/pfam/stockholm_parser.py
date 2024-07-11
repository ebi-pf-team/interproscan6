import re

ACCESSION_LINE = "#=GF AC"
CLAN_MEMBER_MODEL_AC_LINE = "#=GF MB"
NESTING_LINE = "#=GF NE"
ACCESSION_EXTRACTOR_PATTERN = re.compile("^\\#=GF\\s+[A-Z]{2}\\s+([A-Z0-9]+).*$")


def parser_seed_nesting(pfam_a_seed_file: str) -> dict:
    nesting_info = {}
    with open(pfam_a_seed_file, 'rb') as file:
        for line in map(_decode, file):
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
    with open(pfam_clans_file, 'rb') as file:
        for line in map(_decode, file):
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
    name2acc = {}
    acc2nested = {}

    with open(pfam_a_dat_file, 'rb') as reader:
        for line in map(_decode, reader):
            line = line.strip()
            if line.startswith("#=GF ID"):
                name = line.split()[2]
            elif line.startswith("#=GF AC"):
                accession = line.split()[2].split(".")[0]
                name2acc[name] = accession
            elif line.startswith("#=GF NE"):
                try:
                    acc2nested[accession].add(line.split()[2])
                except KeyError:
                    acc2nested[accession] = {line.split()[2]}
            # elif line.startswith("#=GF CL"):
            #     clan = line.split()[2]

    parsed_dat = {}
    for accession, nested in acc2nested.items():
        for name in nested:
            try:
                parsed_dat[accession].append(name2acc[name])
            except KeyError:
                parsed_dat[accession] = [name2acc[name]]
    return parsed_dat


def _decode(b: bytes) -> str:
    """Decode bytes to a string using UTF-8 or Latin-1

    :param b: bytes to decode
    :return: a string, with any trailing whitespace trimmed
    """
    try:
        s = b.decode("utf-8")
    except UnicodeDecodeError:
        pass
    else:
        return s.rstrip()

    return b.decode("latin-1").rstrip()
