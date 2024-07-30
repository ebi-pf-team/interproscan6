import logging
import time

from socket import timeout
from urllib.error import HTTPError, URLError


def lookup_retry_decorator(func):
    """Decorator to re-invoke the wrapped function up to 'args.retries' times."""
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(__name__)
        tries, success, err = 0, False, None

        while not success and (tries < kwargs["retries"]):
            err_message = None
            try:
                result = func(*args, **kwargs)
                success = True
            except (IOError, URLError, timeout, OSError) as err_message:
                success = False
                err = err_message

            tries += 1
            if (not success) and (tries < kwargs["retries"]):
                time.sleep(10)

        if not success:
            logger.error(
                f'Failed to connect to lookup url after {kwargs["retries"]} tries\n'
                f'Error raised: {err}\n'
            )
            raise err
        else:
            return result

    return wrapper
