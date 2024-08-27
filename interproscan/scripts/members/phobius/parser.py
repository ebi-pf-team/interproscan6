import json
import sys
import re

FT_PATTERN = re.compile(
    r"^\w{2}\s+(\w+)\s+(\d+)\s+(\d+)\s*(|\w{11}\.|\w{3}\s\w{11}\.|\w\-\w{6}\.)$"
)
FEATUREDICT = {
    "SIGNAL": {
        "acc": "SIGNAL_PEPTIDE",
        "name": "Signal Peptide",
        "desc": "Signal peptide region",
    },
    "C-REGION.": {
        "acc": "SIGNAL_PEPTIDE_C_REGION",
        "name": "Signal peptide C-region",
        "desc": "C-terminal region of a signal peptide.",
    },
    "H-REGION.": {
        "acc": "SIGNAL_PEPTIDE_H_REGION",
        "name": "Signal peptide H-region",
        "desc": "Hydrophobic region of a signal peptide.",
    },
    "N-REGION.": {
        "acc": "SIGNAL_PEPTIDE_N_REGION",
        "name": "Signal peptide N-region",
        "desc": "N-terminal region of a signal peptide.",
    },
    "TRANSMEM": {
        "acc": "TRANSMEMBRANE",
        "name": "Transmembrane region",
        "desc": (
            "Region of a membrane-bound protein predicted to be "
            "embedded in the membrane."
        ),
    },
    "CYTOPLASMIC.": {
        "acc": "CYTOPLASMIC_DOMAIN",
        "name": "Cytoplasmic domain",
        "desc": (
            "Region of a membrane-bound protein predicted to be "
            "outside the membrane, in the cytoplasm."
        ),
    },
    "NON CYTOPLASMIC.": {
        "acc": "NON_CYTOPLASMIC_DOMAIN",
        "name": "Non cytoplasmic domain",
        "desc": (
            "Region of a membrane-bound protein predicted to be "
            "outside the membrane, in the extracellular region."
        ),
    },
}


class PhobiusHit:
    """Represent protein and associated hits in the Phobius output file.

    We store like domain hits (signal peptide, transmembrane and
    (non-)cytoplasmic domain)) together as one hit with multiple locations.
    In the final output each domain hit is written as a separate
    'signature' hit - this approach ensures all locations are
    retrieved. Each domain type can be found in multiple locations.
    """
    def __init__(self):
        self.seq_id = None  # query protein seq id
        self.signal_peptides = {}  # 4 potential keys: SIGNAL_PEPTIDE, SIGNAL_PEPTIDE_N/H/C_REGION
        self.transmembrane_domains = {}
        self.other_domains = {}

    def get_protein_id(self, line):
        self.seq_id = line.strip().split(maxsplit=1)[1]

    def add_domain(
        self, self_domains,
        acc: str, name: str, desc: str,
        start: str, end: str
    ):
        if acc not in self_domains:
            self_domains[acc] = {
                "member_db": "Phobius",
                "name": name,
                "accession": acc,
                "description": desc,
                "locations": [],
            }
        self_domains[acc]["locations"].append({
            "start": start,
            "end": end,
            "location-fragments": [{
                "start": start,
                "end": end,
                "dc-status": "CONTINUOUS",
            }],
        })

    def parse_domain(self, line: str):
        ftmatch = FT_PATTERN.match(line)
        if ftmatch:
            feature = ftmatch.group(4) if ftmatch.group(4) else ftmatch.group(1)
            start = int(ftmatch.group(2))
            end = int(ftmatch.group(3))
        else:
            raise Exception("Unrecognised line formatting:", line)

        acc, name, desc = FEATUREDICT[feature].values()
        if acc.startswith("SIGNAL_PEPTIDE"):
            self.add_domain(
                self.signal_peptides,
                acc, name, desc, start, end
            )
        elif acc == "TRANSMEMBRANE":
            self.add_domain(
                self.transmembrane_domains,
                acc, name, desc, start, end
            )
        else:
            self.add_domain(
                self.other_domains,
                acc, name, desc, start, end
            )


def main():
    """CL input:
    0. Str repr of path to the output from Phobius
    1. Str repr of path for the output file
    """
    # args 
    args = sys.argv[1:]
    parsed_results = parse(args[0])
    with open(args[1], "w") as fh:
        json.dump(parsed_results, fh)


def parse(phobius_out: str) -> dict:
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
    protein = PhobiusHit()
    with open(phobius_out) as ph_file:
        for line in ph_file:

            if line.startswith("ID"):
                if protein.seq_id:
                    # Only store details of proteins with at least one signal peptide or
                    # transmembrane domain region.
                    # Proteins with only "CYTOPLASMIC" or "NON-CYTOPLASMIC" domains
                    # are junk according to the Phobius domcumentation
                    # -- brought over from interproscan5
                    if protein.signal_peptides or protein.transmembrane_domains:
                        if protein.seq_id not in matches:
                            matches[protein.seq_id] = {}
                        for acc, sp_matches in protein.signal_peptides.items():
                            matches[protein.seq_id][acc] = sp_matches
                        if protein.transmembrane_domains:
                            matches[
                                protein.seq_id
                            ]['TRANSMEMBRANE'] = protein.transmembrane_domains['TRANSMEMBRANE']
                        for acc, other_matches in protein.other_domains.items():
                            matches[protein.seq_id][acc] = other_matches

                protein = PhobiusHit()
                protein.get_protein_id(line)

            elif line.startswith("FT"):
                protein.parse_domain(line.strip())

    return matches


if __name__ == "__main__":
    main()
