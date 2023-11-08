import argparse
import hashlib
import json
import ast


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


def reverse_parse(md5: str, seq_info: str):
    with open(seq_info, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)

    with open(md5, 'r') as md5_data:
        md5_info = md5_data.read()
    checked_seq_md5 = ast.literal_eval(md5_info)
    md5_no_matches = checked_seq_md5['no_matches']

    md52seqinfo = {}
    for seq_id, seq_info in sequence.items():
        md52seqinfo[seq_info[-2]] = f"{seq_info[0]} {seq_info[1]}"

    seq_fasta = ""
    for hash_key in md5_no_matches:
        seq_fasta += f"{md52seqinfo[hash_key]}\n"
    return seq_fasta


def main():
    parser = argparse.ArgumentParser(
        description="sequences parser"
    )
    parser.add_argument(
        "-seq_info", "--sequence_info", type=str, help="sequences to be parsed"
    )
    parser.add_argument(
        "-md5", "--md5",  type=str, required=False, help="md5 list to reverse parse"
    )
    args = parser.parse_args()

    if args.md5:
        md5, seq_info = args.sequence_info
        sequence_parsed = reverse_parse(md5, seq_info)
    else:
        sequences = get_sequences(args.sequence_info)
        sequence_parsed = parse(sequences)
    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
