import json

"""
Use this script to extract just the match data that is used
as the input for testing the parsing of matches when writing
the JSON and XML outputs.
"""

PREFIX = "tests/unit_tests/test_inputs/format_writer/tsv"
SEQ_MAT_PATH = "tests/unit_tests/test_inputs/format_writer/protein/seq_match_dict.json"
MEMBERS = [
    "ANTIFAM",
    "CDD", "COILS",
    "FUNFAM", "GENE3D",
    "HAMAP",
    "MOBIDB", "NCBIFAM",
    "PANTHER", "PFAM", "PHOBIUS",
    "PIRSF", "PIRSR", "PRINTS",
    "PROSITE_PATTERNS", "PROSITE_PROFILES",
    "SFLD", "SIGNALP", "SMART", "SUPERFAMILY"
]


def load_seq_matches() -> dict:
    """Load a JSON file that contains the data found in the 'seqs_matches' variable,
    which is an input argument for build_json_output_protein() and 
    build_json_output_nucleic()"""
    with open(SEQ_MAT_PATH, "r") as fh:
        seq_matches = json.load(fh)
    return seq_matches


def get_match_data(seq_matches: dict) -> None:
    """For each member database, create the seqs_matches dict used as input
    to the TSV and TSV-PRO format writer functions.
    """
    for current_member in MEMBERS:
        member_seq_matches_dict = {}
        # use the data for the seq id that has the most hits and locations
        for seq_id, data in seq_matches.items():
            test_matches = {}  # only retain matches for the member db
            for sig_acc, match_data in data['matches'].items():
                if match_data['member_db'].upper() == current_member:
                    test_matches[sig_acc] = match_data
            if test_matches:
                data_copy = seq_matches[seq_id].copy()
                data_copy['matches'] = test_matches
                member_seq_matches_dict[seq_id] = data_copy

        if member_seq_matches_dict:
            _path = f"{PREFIX}/{current_member}.seq-match-dict.json"
            with open(_path, "w") as fh:
                json.dump(member_seq_matches_dict, fh)
        else:
            print(f"No matches for {current_member}")


if __name__ == "__main__":
    get_match_data(load_seq_matches())
