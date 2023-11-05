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
        "-file", "--input_file", type=str, help="lookup checked md5 dict"
    )
    parser.add_argument(
        "-seq_info", "--seq_info", type=str, required=False, default="", help="sequences parsed"
    )
    args = parser.parse_args()

    print(f"args input file: {args.input_file}")
    print(f"args seq info: {args.seq_info}")
    if args.seq_info:
        sequence_info = reverse_parse(args.input_file, args.seq_info)
        print(f"sequence info do reverse: {sequence_info}")
    else:
        sequences = get_sequences(args.input_file)
        sequence_info = parse(sequences)
        print(f"sequences do get sequences: {sequences}")
        print(f"sequence info do parse normal: {sequence_info}")
    print(json.dumps(sequence_info))


if __name__ == "__main__":
    main()
