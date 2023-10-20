import argparse

import requests


def check_precalc(md5: list, url: str) -> tuple[list, list]:
    precalc = []
    no_precalc = []
    for i in range(0, len(md5), 4):
        sequences_md5 = ', '.join(md5[i:i + 4])
        checkout = requests.get(f"{url}?md5={sequences_md5}")
        is_precalc = checkout.text
        if is_precalc:
            precalc.append(md5)
        else:
            no_precalc.append(md5)
    return precalc, no_precalc


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
    md5 = []
    for seq_id, seq_info in sequences.items():
        md5.append(seq_info[-2])
    md5_upper = [item.upper() for item in md5]
    results = check_precalc(md5_upper, url)
    return results


if __name__ == "__main__":
    main()
