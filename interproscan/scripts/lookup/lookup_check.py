import json
import logging
import sys
import urllib.request

from retry_conn_decorator import lookup_retry_decorator

"""
Checks for pre-calculated matches from any of the member dbs/applications.
This differentiates between cases where the previous calculations found not matches and when a previous calculation has not been performed during LOOKUP_MATCH.
"""


@lookup_retry_decorator
def check_precalc(md5: list, url: str, **kwargs) -> list:
    sequences_md5 = ','.join(md5)
    checkout = urllib.request.urlopen(f"{url}?md5={sequences_md5}")
    is_precalc = checkout.read().decode('utf-8')
    precalc = is_precalc.strip().split("\n")
    return precalc


def main():
    """CL input:
    0. Str repr of path to the JSON file of the input sequences
    1. URL for the MLS
    2. Num of times to retry a failed connection
    3. Str repr of the path for the output file"""
    args = sys.argv[1:]

    sequences = args[0]
    url = args[1]
    retries = int(args[2])
    seq_md5 = []
    with open(sequences, 'r') as seq_data:
        sequences_data = json.load(seq_data)
    for seq_id, seq_info in sequences_data.items():
        seq_md5.append(seq_info[-2].upper())

    md5_checked_matches, err = check_precalc(seq_md5, url, retries=retries)

    if err:
        logging.error(err)
        #return all matches to be calculated locally
        checked_result = {"matches": [],
                          "no_matches": list(seq_md5),
                          "sequences_info": sequences_data}
    else:
        no_matches_md5 = set(seq_md5) - set(md5_checked_matches)
        checked_result = {"matches": md5_checked_matches,
                          "no_matches": list(no_matches_md5),
                          "sequences_info": sequences_data}

    with open(args[3], "w") as fh:
        json.dump(checked_result, fh)


if __name__ == "__main__":
    main()
