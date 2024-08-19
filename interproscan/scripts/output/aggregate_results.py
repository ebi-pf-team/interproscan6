import json
import os
import sys


def aggregate_results(aggreated_results_path: str, result_files: list) -> dict:
    """Parsers a list of strs, each representing a path to a JSON file
    containing the output from an application's analysis.

    Each JSON file is keyed by the protein sequence ID, and valued by
    a dict, keyed by the signature accession and valued by information on the 
    hit (name, evalue, member_db, locations, etc.)

    The content of these JSON files is combined into a single dict
    """
    with open(aggreated_results_path, "r") as fh:
        all_results = json.load(fh)
    for file_path in result_files:
        if file_path:
            with open(file_path, 'r') as file:
                if os.path.getsize(file_path) > 0:
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
    result_files = args[0].strip('[').strip(']').strip().replace(" ", "").split(',')
    aggreated_results_path = args[1]
    all_results = aggregate_results(aggreated_results_path, result_files)
    with open(aggreated_results_path, "w") as fh:
        json.dump(all_results, fh)


if __name__ == "__main__":
    main()
