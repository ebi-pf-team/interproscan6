import json
import sys

from format_writer.tsv_output import tsv_pro_output
from format_writer.tsv_output import tsv_output
from format_writer.json_output import json_output
from format_writer.xml_output import xml_output


def write_results(
    sequences_path: str,
    matches_path: str, 
    output_format: list,
    output_path: str,
    version: str
):
    seq_matches = {}

    all_sequences = {}
    with open(matches_path, 'r') as match_data:
        all_matches = json.load(match_data)
    with open(sequences_path, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
            all_sequences.update(sequence)

    for key, value in all_sequences.items():
        seq_matches[key] = {
            'sequences': value,
            'matches': all_matches.get(key, {})
        }

    if "TSV" in output_format:
        tsv_output(seq_matches, output_path)
    if "TSV-PRO" in output_format:
        tsv_pro_output(seq_matches, output_path)
    if "JSON" in output_format:
        json_output(seq_matches, output_path, version)
    if "XML" in output_format:
        xml_output(seq_matches, output_path, version)


def main():
    args = sys.argv[1:]

    sequences = args[0]
    matches = args[1]
    formats_str = args[2]
    output_path = args[3]
    version = args[4]

    formats = formats_str.upper().split(',')
    write_results(sequences, matches, formats, output_path, version)


if __name__ == "__main__":
    main()
