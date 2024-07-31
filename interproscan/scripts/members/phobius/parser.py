import json
import sys
import re


def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0])
    print(json.dumps(parsed_results, indent=2))


def parse(phobius_out: str) -> dict:
    matches = {}

    featuredict = {
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

    with open(phobius_out) as ph_file:
        for line in ph_file:
            if line.startswith("ID"):
                seq_id = line.strip("ID").split(maxsplit=1)[0]
                locations = []
                matches[seq_id] = {}
            elif line.startswith("FT"):
                line = line.strip("\n")
                ft_pattern = re.compile(
                    r"^(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\s+([\w+\-\s?\.]+)?$"
                )
                ftmatch = ft_pattern.match(line)
                feature = ftmatch.group(2)
                featuretype = None
                start = int(ftmatch.group(3))
                end = int(ftmatch.group(4))
                if ftmatch.group(5):
                    featuretype = ftmatch.group(5)
                featkey = (feature, featuretype)
                acc, name, desc = featuredict[featkey].values()
                match = {
                        "member_db": "Phobius",
                        "version": "",
                        "name": name,
                        "accession": acc,
                        "desc": desc,
                        "locations": locations,
                    }
                location = {
                    "start": start,
                    "end": end,
                    "representative": "false",
                    "location-fragments": {
                        "start": start,
                        "end": end,
                        "dc-status": "CONTINUOUS",
                    },
                }
                if acc not in matches[seq_id]:
                    matches[seq_id][acc] = match
                if matches[seq_id][acc]:
                        matches[seq_id][acc]["locations"].append(location)

    return matches


if __name__ == "__main__":
    main()