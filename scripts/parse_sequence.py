import argparse
import hashlib
import json


def get_sequences(fasta_file: str) -> dict:
    sequences = {}
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                seq_key = line[1:]
                sequences[seq_key] = ""
            else:
                sequences[seq_key] += line
    return sequences


def parse(sequences: dict):
    results = {}
    for key, sequence in sequences.items():
        sequence_info = []
        acc = key.split(" ", 1)[0]
        sequence_info.append(key)
        sequence_info.append(sequence)
        sequence_info.append(hashlib.md5(sequence.encode()).hexdigest().upper())
        sequence_info.append(len(sequence))

        results[acc] = sequence_info
    return results


def main():
    parser = argparse.ArgumentParser(
        description="sequences parser"
    )
    parser.add_argument(
        "-fasta", "--fasta_file", type=str, help="fasta sequences to be parsed"
    )
    args = parser.parse_args()

    sequences = get_sequences(args.fasta_file)
    sequence_parsed = parse(sequences)
    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
