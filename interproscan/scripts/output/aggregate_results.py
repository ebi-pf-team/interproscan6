import json
import os
import sys


def aggregate_results(aggreated_results_path: str, results_file: list) -> dict:
    """Parsers an internal IPS6 JSON, adding the results to the growing file 
    of all aggregated IPS6 hits.

    Each  internal IPS6 JSON file is keyed by the protein sequence ID, and valued by
    a dict, keyed by the signature accession and valued by information on the 
    hit (name, evalue, member_db, locations, etc.)

    :param aggreated_results_path: str repr of path to the growing IPS6 JSON object of all hits
    :param results_file: str repr of path to single IPS6 JSON to be added to aggreated_results
    """
    with open(aggreated_results_path, "r") as fh:
        all_results = json.load(fh)

    with open(results_file.strip(), 'r') as file:
        if os.path.getsize(results_file) > 0:
            try:
                data = json.load(file)
            except json.JSONDecodeError:
                pass

    for seq_id, match_info in data.items():
        # where match_info is a dict
        # keyed by signature accs and valued by hit data
        try:
            all_results[seq_id].update(match_info)
        except KeyError:
            all_results[seq_id] = match_info

    return all_results


def main():
    """
    args:
    0. str representation of list of internal IPS6 JSON files
    1. str repr of path to results_aggregated.json output file
    """
    args = sys.argv[1:]
    aggreated_results_path = args[0]
    results_file = args[1].strip().strip('[').strip(']')
    all_results = aggregate_results(aggreated_results_path, results_file)
    with open(aggreated_results_path, "w") as fh:
        json.dump(all_results, fh)


if __name__ == "__main__":
    main()
