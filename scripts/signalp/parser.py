import json
import re
import sys


COMMENT_LINE = "#"
THRESHOLD = 0.95


def main():
    args = sys.argv[1:]
    parsed_results = parse(args[0])
    print(json.dumps(parsed_results, indent=2))


def parse(signalp_out: str):
    """Parse signalP output into JSON file standardised for InterProScan

    :param signalp_out: path to signalP output CSV file
        ('prediction_results.txt')
    """
    sequence_matches = {}
    with open(signalp_out, "r") as fh:
        for line in fh:
            if line.startswith(COMMENT_LINE):
                continue
            seq_identifer = line.split("\t")[0]
            acc = seq_identifer.split(" ")[0].strip()
            if acc.startswith("sp|"):
                acc = acc.split("|")[1]

            cs_prediction = line.split("\t")[-1]
            if len(cs_prediction) == 1:
                start_location = None
                end_location = None
            elif float(cs_prediction.split("Pr:")[-1].strip()) > THRESHOLD:
                start_location = int(cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[0].strip())
                end_location = int(cs_prediction.split(". ")[0].strip("CS pos: ").split("-")[1].strip())

            sequence_matches[acc] = {"start": start_location, "end": end_location}

    return sequence_matches


if __name__ == "__main__":
    main()
