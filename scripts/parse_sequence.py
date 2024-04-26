import hashlib
import json
import sys


def is_fasta_check(fasta_file: str):
    with open(fasta_file, "r") as f:
        line = f.readline().strip()
        if not line.startswith(">"):
            raise ValueError(f"{fasta_file} is not in FASTA format")


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
        sequence_info.append(hashlib.md5(sequence.encode()).hexdigest())
        sequence_info.append(len(sequence))
        results[acc] = sequence_info
    return results


def main():
    args = sys.argv[1:]
    is_fasta_check(args[0])
    sequences = get_sequences(args[0])
    sequence_parsed = parse(sequences)
    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
