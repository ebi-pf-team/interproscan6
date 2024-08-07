import json
import sys
import re

FT_PATTERN = re.compile(
    r"^(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\s+([\w+\-\s*]+\.)?(\n)$"
)

FEATUREDICT = {
    ("SIGNAL", None): {
        "acc": "SIGNAL_PEPTIDE",
        "name": "Signal Peptide",
        "desc": "Signal peptide region",
    },
    ("DOMAIN", "CYTOPLASMIC."): {
        "acc": "CYTOPLASMIC_DOMAIN",
        "name": "Cytoplasmic domain",
        "desc": "Region of a membrane-bound protein predicted to be outside the membrane, in the cytoplasm.",
    },
    ("DOMAIN", "NON CYTOPLASMIC."): {
        "acc": "NON_CYTOPLASMIC_DOMAIN",
        "name": "Non cytoplasmic domain",
        "desc": "Region of a membrane-bound protein predicted to be outside the membrane, in the extracellular region.",
    },
    ("TRANSMEM", None): {
        "acc": "TRANSMEMBRANE",
        "name": "Transmembrane region",
        "desc": "Region of a membrane-bound protein predicted to be embedded in the membrane.",
    },
    ("DOMAIN", "N-REGION."): {
        "acc": "SIGNAL_PEPTIDE_N_REGION",
        "name": "Signal peptide N-region",
        "desc": "N-terminal region of a signal peptide.",
    },
    ("DOMAIN", "H-REGION."): {
        "acc": "SIGNAL_PEPTIDE_H_REGION",
        "name": "Signal peptide H-region",
        "desc": "Hydrophobic region of a signal peptide.",
    },
    ("DOMAIN", "C-REGION."): {
        "acc": "SIGNAL_PEPTIDE_C_REGION",
        "name": "Signal peptide C-region",
        "desc": "C-terminal region of a signal peptide.",
    },
}


def main():
    # args 0 = fasta file parsed by phobius (to get seq lens)
    # args 1 = output form phobius
    args = sys.argv[1:]
    load_seqs = args[0]
    parsed_results = parse(args[0], load_seqs)
    print(json.dumps(parsed_results, indent=2))


def load_seqs(fasta: str) -> dict:
    seqs = {}
    current_seq = ""
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                current_seq = line.strip(">")
                seqs[current_seq] = 0
            else:
                seqs[current_seq] += len(line.strip())
    return seqs


def parse(phobius_out: str, seqs_dict: dict) -> dict:
    """Example Phobius hit from the output file. Data is grouped by the query protein seq ID:
     //
    ID   sp|O16846|1AK_HETMG
    FT   SIGNAL        1     22       
    FT   DOMAIN        1      5       N-REGION.
    FT   DOMAIN        6     17       H-REGION.
    FT   DOMAIN       18     22       C-REGION.
    FT   DOMAIN       23     74       NON CYTOPLASMIC.
    //
    """
    matches = {}
    version = phobius_out.split("._.")[0]
    seq_id = None
    with open(phobius_out) as ph_file:
        for line in ph_file:
            if line.startswith("ID"):
                if seq_id:
                    # drop previous protein if any location covers the seq entirely
                    # filtering step brought over from i5
                    for location in matches[seq_id][]
                seq_id = line.strip("ID").split(maxsplit=1)[0]
                matches[seq_id] = {}
            elif line.startswith("FT"):
                ftmatch = FT_PATTERN.match(line)
                if ftmatch:
                    feature = (ftmatch.group(2), ftmatch.group(5) if ftmatch.group(6) else None)
                    start = int(ftmatch.group(3))
                    end = int(ftmatch.group(4))
                else:
                    raise Exception("Unrecognised line formatting:", line)
                acc, name, desc = FEATUREDICT[feature].values()
                match = {
                    "member_db": "Phobius",
                    "version": version,
                    "name": name,
                    "accession": acc,
                    "description": desc,
                    "locations": [],
                }
                location = {
                    "start": start,
                    "end": end,
                    "representative": "false",
                    "location-fragments": [{
                        "start": start,
                        "end": end,
                        "dc-status": "CONTINUOUS",
                    }],
                }

                if acc not in matches[seq_id]:
                    matches[seq_id][acc] = match
                matches[seq_id][acc]["locations"].append(location)

    return matches


if __name__ == "__main__":
    main()
