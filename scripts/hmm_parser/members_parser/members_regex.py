import re


def get_accession_regex(appl):
    if appl == "antifam":
        accession_regex = re.compile("^Accession:\\s+(ANF\\d{5})\\s*$")
        return accession_regex
    if appl == "ncbifam":
        accession_regex = re.compile("^Accession:\\s+((TIGR|NF)\\d+)\\.\\d+$")
        return accession_regex
