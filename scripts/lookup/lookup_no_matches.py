import argparse
import json


def md52fasta(md5: set, sequence: dict):
    md52seqinfo = {}
    for seq_id, seq_info in sequence.items():
        md52seqinfo[seq_info[-2]] = f"{seq_info[0]}\n{seq_info[1]}"

    seq_fasta = ""
    for hash_key in md5:
        seq_fasta += f">{md52seqinfo[hash_key]}\n"
    return seq_fasta


def main():
    parser = argparse.ArgumentParser(
        description="Parse to no match lookup"
    )
    parser.add_argument(
        "-checked", "--checked_lookup", type=str, help="dict with md5 lookup checked"
    )
    args = parser.parse_args()

    with open(args.checked_lookup, 'r') as md5_data:
        checked_data = json.load(md5_data)
    no_matches = checked_data["no_matches"]
    sequences_data = checked_data["sequences_info"]

    fasta_to_analyse = ""
    if no_matches:
        fasta_to_analyse = md52fasta(no_matches, sequences_data)

    print(fasta_to_analyse)


if __name__ == "__main__":
    main()
