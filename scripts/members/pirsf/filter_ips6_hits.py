import argparse
import json
import re

from pathlib import Path


CHILD_DAT_LINE = re.compile(r">\D+\d+\schild:\s(.*)$")


class DatEntry:
    def __init__(self):
        self.model_id = None
        self.name = None
        self.mean_l = None
        self.std_l = None
        self.min_s = None
        self.mean_s = None
        self.std_s = None
        self.blast = None
        self.children = []

    def add_model_id(self, value: str):
        self.model_id = value.split()[0].replace(">", "")
        line_match = CHILD_DAT_LINE.match(value)
        if line_match:
            self.children = line_match.group(1).split()

    def add_name(self, value: str):
        self.name = value

    def add_data_values(self, value: str):
        self.mean_l, self.std_l, self.min_s, self.mean_s, self.std_s = value.split()

    def add_blast(self, value: str):
        self.blast = 0 if value.split()[1] == "NO" else 1


def build_parser() -> argparse.ArgumentParser:
    """Build cmd-line argument parser"""
    parser = argparse.ArgumentParser(
        prog="pirsf_post_processing",
        description="Reimplementation of pl script and mod from i5 to post-processing PIRSF hits",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "ips6",
        type=Path,
        help="Path to internal IPS6 JSON containing HMMER matches"
    )
    parser.add_argument(
        "dat",
        type=Path,
        help="Path to PIRSF.dat file"
    )
    parser.add_argument(
        "out_json",
        type=Path,
        help="Path to write output JSON"
    )
    parser.add_argument('-verbose', action='store_true', help="Report No matches")
    parser.add_argument('-outfmt', choices=['pirsf', 'ips6'], default='pirsf', help="Output format")
    parser.add_argument('-cpus', type=int, default=1, help="Number of CPUs to use")
    parser.add_argument('-tmpdir', default='tmp', help="Directory for HMMER to use")

    return parser


def load_dat(dat_path: str) -> tuple[dict, dict]:
    """Load data from pirsf.dat file

    :param dat_path: str repr of path to pirsf.dat file

    Return dict with all data entries
    and dict with {subfamily: parent family}
    """
    stage = ''
    entry = None
    dat_entries = {}  # model id: datEntry
    children = {}  # subfamily: parent family
    with open(dat_path, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if entry:
                    dat_entries[entry.model_id] = entry
                    for child in entry.children:
                        children[child] = entry.model_id
                entry = DatEntry()
                entry.add_model_id(line)
                stage = 'GET_NAME'
            elif stage == 'GET_NAME':
                entry.add_name(line)
                stage = 'GET_VALUES'
            elif stage == 'GET_VALUES':
                entry.add_data_values(line)
                stage = 'GET_BLAST'
            elif stage == 'GET_BLAST':
                entry.add_blast(line)
                stage = 'GET_MODEL_ID'

    return dat_entries, children


def get_location_data(match: dict):
    seq_len, seq_start, seq_end, hmm_start, hmm_end, score = [0]*6
    for location in match["locations"]:
        seq_len = location["hmmLength"]
        seq_start = location["start"]
        seq_end = location["end"]
        hmm_start = location["hmmStart"]
        hmm_end = location["hmmEnd"]
        score += location["score"]

        # update to match georgetown 2017 script
        if location["start"] < seq_start and location["hmmStart"] < hmm_start:
            seq_start = location["start"]
            hmm_start = location["hmmStart"]
        if location["end"] < seq_end and location["hmmEnd"] < hmm_end:
            seq_end = location["end"]
            hmm_end = location["hmmEnd"]

    return seq_len, seq_start, seq_end, hmm_start, hmm_end, score


def store_match(
    match: dict,
    matches_dict: dict,
    protein_id: str,
    model_id: str,
    seq_len: int,
    seq_start: int,
    seq_end: int,
    hmm_start: int,
    hmm_end: int,
    score: float
) -> dict:
    """Add a match that passes the filtering criteria in the processed_matches dict"""
    if protein_id not in matches_dict:
        matches_dict[protein_id] = {}

    if model_id not in matches_dict[protein_id]:
        matches_dict[protein_id][model_id] = match
        matches_dict[protein_id][model_id]["locations"] = [{
            "start": seq_start,
            "end": seq_end,
            "hmmStart": hmm_start,
            "hmmEnd": hmm_end,
            "hmmLength": seq_len,
            "score": score,
            "postProcessed": "false"
        }]

    return matches_dict


def process_results(ips6_json: Path, pirsf_dat: dict, children: dict) -> dict:
    """
    :param ips6_json: path to internal IPS6 JSON file with hmmer hits
    :param pirsf_dat: dict containing data from the pirsf.dat file
    :param children: dict containing {subfam: parent fam}
    """
    with open(ips6_json, "r") as fh:
        matches = json.load(fh)
    processed_matches = {}

    for protein_id in matches:
        for model_id in matches[protein_id]:
            seq_len, seq_start, seq_end, hmm_start, hmm_end, score = get_location_data(
                matches[protein_id][model_id]
            )

            # calculate ratios and deviations
            ovl = (seq_end - seq_start + 1) / seq_len
            ratio = (hmm_end - hmm_start + 1) / (seq_end - seq_start + 1)
            ld = abs(seq_len - pirsf_dat[model_id].mean_l)

            # check if fam (has children) or subfamily (does not have children)
            if model_id in children:
                # If a sub-family, process slightly differently. Only consider the score.
                if ratio > 0.67 and score >= pirsf_dat[model_id].min_s:
                    # subfamily passes the criteria
                    processed_matches = store_match(
                        matches[protein_id][model_id],
                        processed_matches,
                        protein_id, model_id,
                        seq_len, seq_start, seq_end, hmm_start, hmm_end, score
                    )

                    # if we have a subfamily match we should consider the parent to also be a match
                    parent = children[model_id]  # parent = parent model id
                    if parent not in processed_matches[protein_id]:
                        seq_len, seq_start, seq_end, hmm_start, hmm_end, score = get_location_data(
                            matches[protein_id][parent]
                        )
                        processed_matches = store_match(
                            matches[protein_id][parent],
                            processed_matches,
                            protein_id, model_id,
                            seq_len, seq_start, seq_end, hmm_start, hmm_end, score
                        )

            elif (ratio > 0.67 and ovl >= 0.8 and score >= pirsf_dat[model_id].min_s and \
                  (ld < 3.5 * pirsf_dat[model_id].std_l or ld < 50)):
                # everything passes the threshold of length, score and standard deviations of length
                processed_matches = store_match(
                    matches[protein_id][model_id],
                    processed_matches,
                    protein_id, model_id,
                    seq_len, seq_start, seq_end, hmm_start, hmm_end, score
                )

    return processed_matches


def main():
    parser = build_parser()
    args = parser.parse_args()

    pirsf_dat, children = load_dat(args.dat)
    procossed_results = process_results(args.ips6, pirsf_dat, children)

    with open(args.out_json, "w") as fh:
        json.dump(procossed_results, fh, indent=2)


if __name__ == "__main__":
    main()
