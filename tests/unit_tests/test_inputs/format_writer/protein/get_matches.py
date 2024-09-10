import json

"""
Use this script to extract just the match data that is used
as the input for testing the parsing of matches when writing
the JSON and XML outputs.
"""

PREFIX = "tests/unit_tests/test_inputs/format_writer/matches"
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
    """For each member database, find the sequence that has the most signature hits and 
    locations per signature. For this sequence, only retain the signature hits relevant
    to the current member being processed, and write out the results to a file.
    This JSON file contains the data that will be used for the input argument
    'data' in the functions get_matches() in the JSON and XML format_writers.
    """
    for current_member in MEMBERS:
        counts = {}
        # use the data for the seq id that has the most hits and locations
        for seq_id, data in seq_matches.items():
            signatures = {}  # match structure of data['matches'] but only with sigs for the current db
            locations = 0
            # data is the var passed by the get_matches() func in the json+xml format writers
            # count number of hits and number of locations per seq_id
            for sig_acc, match_data in data['matches'].items():
                if match_data['member_db'].upper() == current_member:
                    signatures[sig_acc] = match_data
                    locations += len(match_data['locations'])
            if signatures:
                counts[seq_id] = {'signatures': signatures, 'locations': locations}

        if counts:
            # sort by most signatures and locations 
            sorted_seqids = dict(sorted(
                counts.items(),
                key=lambda item: (len(item[1]['signatures']), item[1]['locations']),
                reverse=True
            ))

            # 'data' var passed by the get_matches() func in the json+xml format writers
            # containing data for only one member db
            seq_id = list(sorted_seqids.keys())[0]
            test_matches = sorted_seqids[seq_id]['signatures']
            data_copy = seq_matches[seq_id].copy()
            data_copy['matches'] = test_matches
            _path = f"{PREFIX}/{current_member}.match-data.json"
            with open(_path, "w") as fh:
                json.dump(data_copy, fh)

        else:
            print(f"No matches for {current_member}")


if __name__ == "__main__":
    get_match_data(load_seq_matches())
