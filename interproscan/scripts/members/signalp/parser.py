import json
import sys

COMMENT_LINE = "#"


def main():
    """
    Input params for SignalP parser:
    0. SignalP output gff3 file path
    1. SignalP output csv file path
    2. SignalP p-value threshold
    3. SignalP version
    4. SignalP organism
    5. Str repr of path for the output file
    """
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1], float(args[2]), args[3], args[4])
    with open(args[5], "w") as fh:
        json.dump(parsed_results, fh)


def parse(signalp_out: str, signalp_cs: str, threshold: float, signalp_version: str, signalp_db: str):
    """Parse signalP output into JSON file standardised for InterProScan

    :param signalp_out: path to signalP signal peptide location output file
        ('output.gff3')
    :param signalp_cs: path to signalP cleavage site location output CSV file
        ('prediction_results.txt')
    :param threshold: p-value threshold that must be met or exceeded
    :param signalp_version: str, version num of signalP e.g. '6.0h'
    :param signalp_db: str, organism arg, used to assign member database
    """
    sequence_matches = {}

    member_db = "SignalP" if signalp_db == "other" else "SignalP_EUK"
    org = "Other" if signalp_db == "other" else "Eukarya"

    with open(signalp_out, "r") as f:
        for line in f:
            if line.startswith(COMMENT_LINE):
                continue
            seq_identifer = line.split("\t")[0]
            acc = seq_identifer.split(" ")[0].strip()
            pvalue = line.split("\t")[5]

            if len(pvalue) > 1:
                if float(pvalue) >= threshold:
                    start = line.split("\t")[3]
                    end = line.split("\t")[4]
                    pvalue = float(pvalue)
                else:
                    start, end, pvalue = None, None, None
            # checks if there is only one signal peptide prediction per protein
            if acc in sequence_matches:
                raise Exception(f"Protein {acc} has more than one SignalP match")
            # reports signal peptide start and end location
            if start:
                sequence_matches[acc] = {
                    "signal_peptide": {
                        "accession": "SignalP",
                        "name": "SignalP",
                        "member_db": member_db,
                        "version": signalp_version,
                        "orgType": org,
                        "model-ac": "SignalP",
                        "locations": [{
                            "start": start,
                            "end": end,
                            "pvalue": pvalue,
                            "cleavage_start": "",
                            "cleavage_end": "",
                            "representative": "false"
                        }]
                    }
                }

    matches = get_cleavage_site(signalp_cs, sequence_matches)

    return matches


def get_cleavage_site(signalp_cs: str, seq_matches: dict):
    acc_list = []
    with open(signalp_cs, "r") as fh:
        for line in fh:
            if line.startswith(COMMENT_LINE):
                continue
            seq_identifer = line.split("\t")[0]
            acc = seq_identifer.split(" ")[0].strip()
            cs_start = None
            cs_end = None

            cs_prediction = line.split("\t")[-1]
            if len(cs_prediction) > 1:
                cs_start = int(
                    cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[
                        0].strip())
                cs_end = int(
                    cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[
                        1].strip())

            if cs_start and acc in seq_matches:
                acc_list.append(acc)
                for location in seq_matches[acc]["signal_peptide"]["locations"]:
                    location["cleavage_start"] = int(cs_start)
                    location["cleavage_end"] = int(cs_end)

            elif acc in seq_matches and not cs_start:
                raise Exception(f"{acc} is missing cleavage start in prediction_results.txt")
            # if acc not in seq matches, signalp prediction did not pass p-value threshold

    return seq_matches


if __name__ == "__main__":
    main()
