import argparse
import json
import requests


def check_precalc(md5: list, url: str) -> list:
    sequences_md5 = ', '.join(md5)
    checkout = requests.get(f"{url}?md5={sequences_md5}")
    is_precalc = checkout.text
    precalc = is_precalc.strip().split("\n")
    return precalc


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
    with open(sequences, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
    for seq_id, seq_info in sequence.items():
        seq_md5.append(seq_info[-2])
    md5_checked_matches = check_precalc(seq_md5, url)
    no_matches_md5 = set(seq_md5) - set(md5_checked_matches)
    check_result = {"matches": md5_checked_matches, "no_matches": no_matches_md5}
    print(check_result)


if __name__ == "__main__":
    main()
