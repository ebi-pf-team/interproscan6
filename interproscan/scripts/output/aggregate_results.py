import json
import sys

from pathlib import Path


class FileFormatError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


def aggregate_results(aggreated_results_path: str, results_file: Path) -> dict:
    """Parsers an internal IPS6 JSON, adding the results to the growing file 
    of all aggregated IPS6 hits.

    Each  internal IPS6 JSON file is keyed by the protein sequence ID, and valued by
    a dict, keyed by the signature accession and valued by information on the 
    hit (name, evalue, member_db, locations, etc.)

    :param aggreated_results_path: str repr of path to the growing IPS6 JSON object of all hits
    :param results_file: str repr of path to single IPS6 JSON to be added to aggreated_results
    """
    def _multiple_lines(file_handle):
        lines = 0
        for line in file_handle.readlines():
            lines += 1
            if lines > 10:
                continue
        return True if lines > 1 else False

    data = {}

    with open(aggreated_results_path, "r") as fh:
        all_results = json.load(fh)

    with open(results_file, 'r') as file:
        try:
            data = json.load(file)
        except json.JSONDecodeError as exc:
            # are we dealing with an empty file/json obj then pass
            # or is the file incorrectly formatted/corrupted - then raise an error
            with open(results_file, 'r') as fh:
                if not _multiple_lines(fh):
                    if fh.readline().strip() == "{}" or fh.readline().strip() == "":
                        return all_results
                raise FileFormatError(
                    (f"Internal results file '{results_file}' is incorrectly formatted "
                     "or corrupted")
                ) from exc

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
    aggreated_results_path = Path(args[0])
    if not aggreated_results_path.exists():
        raise FileNotFoundError(
            f"Could not find aggregated results file {str(aggreated_results_path)} during\n"
            "the AGGREGATE_RESULTS process"
        )

    results_file = Path(args[1].strip().strip('[').strip(']').strip())
    if results_file.exists() and results_file.stat().st_size > 0:
        all_results = aggregate_results(aggreated_results_path, results_file)
    else:
        raise FileNotFoundError(
            f"Could not find internal results file {str(results_file)} during\n"
            "the AGGREGATE_RESULTS process"
        )

    with open(aggreated_results_path, "w") as fh:
        json.dump(all_results, fh)


if __name__ == "__main__":
    main()
