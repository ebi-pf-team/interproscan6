import re


def get_accession_regex(appl):
    if appl == "antifam":
        return re.compile("^Accession:\\s+(ANF\\d{5})\\s*$")
    if appl == "ncbifam":
        return re.compile("^Accession:\\s+((TIGR|NF)\\d+)\\.\\d+$")
