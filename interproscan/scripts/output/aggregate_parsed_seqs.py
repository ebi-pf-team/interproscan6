import json
import sys


def aggregate_parsed_seqs(parsed_seq_paths: list) -> dict:
    parsed_seqs = {}
    for file_path in parsed_seq_paths:
        with open(file_path, "r") as fh:
            data = json.load(fh)
            for seq_id, seq_data in data.items():
                parsed_seqs[seq_id] = seq_data

    return parsed_seqs


def main():
    """CL input
    0. Str repr of a list of paths to parsed_sequences.py output JSON files
    1. Str repr of the path for the output file"""
    args = sys.argv[1:]
    parsed_seq_files = args[0].strip('[]').replace(" ", "").split(',')
    all_seqs = aggregate_parsed_seqs(parsed_seq_files)
    with open(args[1], "w") as fh:
        json.dump(all_seqs, fh)


if __name__ == "__main__":
    main()
