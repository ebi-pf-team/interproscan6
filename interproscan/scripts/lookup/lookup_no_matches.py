import json
import sys


def md52fasta(md5: set, sequence: dict):
    md52seqinfo = {}
    for seq_id, seq_info in sequence.items():
        md52seqinfo[seq_info['md5']] = f"{seq_info['seq_id']}\n{seq_info['sequence']}"

    seq_fasta = ""
    for hash_key in md5:
        try:
            seq_fasta += f">{md52seqinfo[hash_key]}\n"
        except KeyError:
            seq_fasta += f">{md52seqinfo[hash_key.lower()]}\n"
    return seq_fasta


def main():
    args = sys.argv[1:]

    with open(args[0], 'r') as md5_data:
        checked_data = json.load(md5_data)
    no_matches = checked_data["no_matches"]
    sequences_data = checked_data["sequences_info"]

    fasta_to_analyse = ""
    if no_matches:
        fasta_to_analyse = md52fasta(no_matches, sequences_data)

    print(fasta_to_analyse)


if __name__ == "__main__":
    main()
