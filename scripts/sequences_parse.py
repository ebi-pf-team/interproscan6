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


def parse(sequences):
    results = []
    for key, sequence in sequences.items():
        sequence_info = {}
        id_desc = key.split(" ", 1)
        sequence_info["id"] = id_desc[0]
        sequence_info["id_desc"] = key
        sequence_info["sequence"] = sequence
        sequence_info["md5"] = hashlib.md5(sequence.encode()).hexdigest()
        sequence_info["length"] = len(sequence)

        results.append(sequence_info)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="sequences parser"
    )
    parser.add_argument(
        "-seq", "--sequences", type=str, help="fasta file with sequences"
    )
    args = parser.parse_args()

    sequences = get_sequences(args.sequences)
    sequence_info = parse(sequences)
    print(json.dumps(sequence_info))


if __name__ == "__main__":
    main()
