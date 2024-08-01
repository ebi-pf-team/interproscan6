import json
import sys


COMMENT_LINE = "#"


def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0], float(args[1]), args[2], args[3])
    print(json.dumps(parsed_results, indent=2))


def parse(signalp_out: str, threshold: float, signalp_version: str, signalp_db = str):
    """Parse signalP output into JSON file standardised for InterProScan

    :param signalp_out: path to signalP output CSV file
        ('prediction_results.txt')
    :param threshold: p-value threshold that must be met or exceeded
    :param signalp_version: str, version num of signalP e.g. '6.0h'
    """
    sequence_matches = {}

    if signalp_db == "other":
        member_db =  "SignalP"
    elif signalp_db == "euk":
        member_db = "SignalP_EUK"

    with open(signalp_out, "r") as fh:
        for line in fh:
            if line.startswith(COMMENT_LINE):
                continue
            seq_identifer = line.split("\t")[0]
            acc = seq_identifer.split(" ")[0].strip()
            start_location = None
            end_location = None
            pvalue = None

            cs_prediction = line.split("\t")[-1]
            if len(cs_prediction) > 1:
                if float(cs_prediction.split("Pr:")[-1].strip()) >= threshold:
                    start_location = int(cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[0].strip())
                    end_location = int(cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[1].strip())
                    pvalue = float(cs_prediction.split("Pr:")[-1].strip())

            if start_location:
                sequence_matches[acc] = {
                    "signal_peptide": {
                        "member_db": member_db,
                        "version": signalp_version,
                        "locations": [{
                            "start": start_location,
                            "end": end_location,
                            "pvalue": pvalue
                        }]
                    }
                }

    return sequence_matches


if __name__ == "__main__":
    main()
