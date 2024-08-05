import json
import sys


COMMENT_LINE = "#"


def main():
    """
    Input params for SignalP parser:
    1. SignalP output gff3 file path
    2. SignalP output csv file path
    3. SignalP p-value threshold
    4. SignalP version
    5. SignalP organism
"""
    args = sys.argv[1:]
    parsed_results = parse(args[0], args[1], float(args[2]), args[3], args[4])
    print(json.dumps(parsed_results, indent=2))


def parse(signalp_out: str, signalp_cs: str, threshold: float, signalp_version: str, signalp_db = str):
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

    with open(signalp_out, "r") as f:
        for line in f:
            if line.startswith(COMMENT_LINE):
                continue
            seq_identifer = line.split("\t")[0]
            acc = seq_identifer.split(" ")[0].strip()
            start = None
            end = None
            pvalue = None
            predict_pvalue = line.split("\t")[5]

            if len(predict_pvalue) > 1:
                if float(predict_pvalue) >= threshold:
                    start = line.split("\t")[3]
                    end = line.split("\t")[4]
                    pvalue = float(predict_pvalue)

            if acc in sequence_matches:
                raise Exception(f"Protein {acc} had more than one SignalP match")

            if start:
                sequence_matches[acc] = {
                    "signal_peptide": {
                        "member_db": member_db,
                        "version": signalp_version,
                        "locations": [{
                            "start": start,
                            "end": end,
                            "pvalue": pvalue,
                        }]
                    }
                }


    with open(signalp_cs, "r") as fh:
        for line in fh:
            if line.startswith(COMMENT_LINE):
                continue
            seq_identifer = line.split("\t")[0]
            acc = seq_identifer.split(" ")[0].strip()
            cs_start = None
            cs_end = None
            pvalue = None

            cs_prediction = line.split("\t")[-1]
            if len(cs_prediction) > 1:
                # reports start and end of cleavage site
                # to do: report signal peptide location
                if float(cs_prediction.split("Pr:")[-1].strip()) >= threshold:
                    cs_start = int(cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[0].strip())
                    cs_end = int(cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[1].strip())
                    pvalue = float(cs_prediction.split("Pr:")[-1].strip())
            # add check if multiple signal peptides predicted, raise error
            if cs_start:
                for location in sequence_matches[acc]["signal_peptide"]["locations"]:
                    if location["pvalue"] == pvalue:
                        location["cleavage_start"] = cs_start
                        location["cleavage_end"] = cs_end


    return sequence_matches


if __name__ == "__main__":
    main()
