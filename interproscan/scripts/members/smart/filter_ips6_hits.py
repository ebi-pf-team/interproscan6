import json
import re
import sys


SMART_SER_THR_KINASE_METHOD = "SM00220"
SMART_SER_THR_REGEX = re.compile(".*D[LIVM]K\\w\\wN.*")

SMART_TYR_KINASE_METHOD = "SM00219"
SMART_TYR_REGEX = re.compile(".*HRD[LIV][AR]\\w\\wN.*")


def load_protein_seqs(fasta: str) -> dict:
    """Load protein sequences into memory for kinase filtering.

    Keyed by protein ID, valued by str repr of protein seq
    :param fasta: str repr of path to fasta file of query prot seqs
    """
    protein_seqs = {}
    current_id = None
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                current_id = line.strip(">").split(maxsplit=1)[0]
                protein_seqs[current_id] = ""
            else:
                protein_seqs[current_id] += line.strip()
    return protein_seqs


def filter_matches(matches: dict, protein_seqs: dict) -> dict:
    """
    Additional filtering to differentiate between serine/threonine kiases and 
    tyrosine kinases when both types are matches in a protein sequence.

    Someone decided in 2012 that if both serine-theronine (SM00220) and tyrosine kinase 
    (SM00219) matches are present we only keep both matches if the protein sequences
    matches the respective REs for the domains.

    Using the unlicsened distibution (which will be most users),
    keep all matches and only apply additional checks for matches that contain
    both a serine-theronine (SM00220) and tyrosine kinase (SM00219) matches.
    """
    for protein_id in matches:
        if SMART_SER_THR_KINASE_METHOD in matches[protein_id] and SMART_TYR_KINASE_METHOD in matches[protein_id]:
            protein_seq = protein_seqs[protein_id]
            if not protein_seq.match(SMART_SER_THR_REGEX):
                del matches[protein_id][SMART_SER_THR_KINASE_METHOD]
            if not protein_seq.match(SMART_TYR_REGEX):
                del matches[protein_id][SMART_TYR_KINASE_METHOD]

    return matches


def main():
    """System arguments:
    1. str repr of the path the internal IPS6 JSON file from the hmmpfam parser
    2. str repr of the path to the fasta file of the query protein seqs
    3. str repr of path to write out the final results
    """
    with open(sys.argv[1], "r") as fh:
        matches = json.load(fh)

    protein_seqs = load_protein_seqs(sys.argv[2])
    with open("seqs.dict.json", "w") as fh:
        json.dump(protein_seqs, fh, indent=2)

    parsed_matches = filter_matches(matches, protein_seqs)

    with open(sys.argv[3], "w") as fh:
        json.dump(parsed_matches, fh, indent=2)


if __name__ == "__main__":
    main()
