import hashlib
import json
import sys


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
    args = sys.argv[1:]
    sequences = get_sequences(*args)
    sequence_parsed = parse(sequences)
    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
