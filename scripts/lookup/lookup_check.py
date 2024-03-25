import json
import requests
import sys


def check_precalc(md5: list, url: str) -> list:
    sequences_md5 = ', '.join(md5)
    checkout = requests.get(f"{url}?md5={sequences_md5}")
    is_precalc = checkout.text
    precalc = is_precalc.strip().split("\n")
    return precalc


def main():
    args = sys.argv[1:]

    sequences = args[0]
    url = args[1]
    seq_md5 = []
    with open(sequences, 'r') as seq_data:
        sequences_data = json.load(seq_data)
    for seq_id, seq_info in sequences_data.items():
        seq_md5.append(seq_info[-2])
    md5_checked_matches = check_precalc(seq_md5, url)
    no_matches_md5 = set(seq_md5) - set(md5_checked_matches)
    checked_result = {"matches": md5_checked_matches,
                      "no_matches": list(no_matches_md5),
                      "sequences_info": sequences_data}
    print(json.dumps(checked_result))


if __name__ == "__main__":
    main()
