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


def reverse_parse(md5: dict, seq_info: str):
    with open(seq_info, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)

    no_matches_md5 = md5["no_matches"]

    md52seqinfo = {}
    for seq_id, seq_info in sequence.items():
        md52seqinfo[seq_info[-2]] = f"{seq_info[0]} {seq_info[1]}"

    seq_fasta = ""
    for md5 in no_matches_md5:
        seq_fasta += f"{md52seqinfo[md5]}\n"
    return seq_fasta


def main():
    parser = argparse.ArgumentParser(
        description="sequences parser"
    )
    parser.add_argument(
        "-file", "--input_file", type=str, help="file to process parse/reverse parse"
    )
    parser.add_argument(
        "-reverse", "--reverse", type=bool, help="flag to normal or reverse parse process"
    )
    args = parser.parse_args()

    if args.reverse:
        sequence_info = reverse_parse(args.input_file)
    else:
        sequences = get_sequences(args.input_file)
        sequence_info = parse(sequences)
    print(json.dumps(sequence_info))


if __name__ == "__main__":
    main()
