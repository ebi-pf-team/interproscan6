import argparse
import json
import requests


def check_precalc(md5: list, url: str) -> list:
    sequences_md5 = ', '.join(md5)
    checkout = requests.get(f"{url}?md5={sequences_md5}")
    is_precalc = checkout.text
    precalc = is_precalc.strip().split("\n")
    return precalc


def md52fasta(md5: set, sequence: dict):
    md52seqinfo = {}
    for seq_id, seq_info in sequence.items():
        md52seqinfo[seq_info[-2]] = f"{seq_info[0]} {seq_info[1]}"

    seq_fasta = ""
    for hash_key in md5:
        seq_fasta += f"{md52seqinfo[hash_key]}\n"
    return seq_fasta


def main():
    parser = argparse.ArgumentParser(
        description="Check if sequence is pre calculated"
    )
    parser.add_argument(
        "-seq", "--sequences", type=str, help="sequences hash"
    )
    parser.add_argument("-url", "--url", type=str, help="url to check precalc lookup")
    args = parser.parse_args()

    sequences = args.sequences
    url = args.url
    seq_md5 = []
    fasta_to_analyse = None
    with open(sequences, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
    for seq_id, seq_info in sequence.items():
        seq_md5.append(seq_info[-2])
    md5_checked_matches = check_precalc(seq_md5, url)
    no_matches_md5 = set(seq_md5) - set(md5_checked_matches)
    if no_matches_md5:
        fasta_to_analyse = md52fasta(no_matches_md5, sequence)
    checked_result = {"matches": md5_checked_matches, "no_matches": fasta_to_analyse}
    print(checked_result)


if __name__ == "__main__":
    main()
