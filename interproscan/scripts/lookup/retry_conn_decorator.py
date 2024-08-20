import time

from socket import timeout
from urllib.error import URLError


def lookup_retry_decorator(func):
    """Decorator to re-invoke the wrapped function up to 'args.retries' times."""
    def wrapper(*args, **kwargs):
        tries, success, err = 0, False, None
        result = []

        while not success and (tries < kwargs["retries"]):
            try:
                result = func(*args, **kwargs)
                success = True
            except (IOError, URLError, timeout, OSError) as err_message:
                success = False
                err = err_message

            tries += 1
            if (not success) and (tries < kwargs["retries"]):
                time.sleep(10)

        return result, err

    return wrapper
