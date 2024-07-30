import json
import sys
import urllib.request

"""
Checks for pre-calculated matches from any of the member dbs/applications.
This differentiates between cases where the previous calculations found not matches and when a previous calculation has not been performed during LOOKUP_MATCH.
"""


def lookup_decorator(func):
    """Decorator to re-invoke the wrapped function up to 'args.retries' times."""

    def wrapper(*args, **kwargs):
        logger = logging.getLogger(__name__)
        tries, success, err = 0, False, None

        while not success and (tries < kwargs['max_tries']):
            # reset storing error messsage
            err_message = None

            try:
                func(*args, **kwargs)

            except (
                    IOError,
                    HTTPError,
                    URLError,
                    timeout,
                    ConnectionError,
                    OSError,
                    MissingSchema,
                    RequestError,
            ) as err_message:
                success = False
                err = err_message

            if err is None:
                success = True

            tries += 1

            if (not success) and (tries < kwargs['max_tries']):
                logger.warning(
                    f'Failed to connect to {url} on try {tries}/{kwargs["max_tries"]}\n'
                    f'Error raised: {err}\n'
                    'Retrying connection in 10s'
                )
                time.sleep(10)

        if success is False:
            logger.warning(
                f'Failed to connect to {url} after {kwargs["max_tries"]} tries\n'
                f'Error raised: {err}\n'
            )
            return err
        else:
            return None

    return wrapper


@lookup_decorator
def check_precalc(md5: list, url: str) -> list:
    sequences_md5 = ','.join(md5)
    checkout = urllib.request.urlopen(f"{url}?md5={sequences_md5}")
    is_precalc = checkout.read().decode('utf-8')
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
        seq_md5.append(seq_info[-2].upper())
    md5_checked_matches = check_precalc(seq_md5, url)
    no_matches_md5 = set(seq_md5) - set(md5_checked_matches)
    checked_result = {"matches": md5_checked_matches,
                      "no_matches": list(no_matches_md5),
                      "sequences_info": sequences_data}
    print(json.dumps(checked_result))


if __name__ == "__main__":
    main()
