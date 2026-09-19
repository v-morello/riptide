import logging
import time
from functools import wraps


def timing(func, *args, **kwargs):
    """Decorate a function to log its runtime."""

    @wraps(func)
    def wrapped(*args, **kwargs):
        log = logging.getLogger("riptide.timing")
        t0 = time.time()
        output = func(*args, **kwargs)
        t1 = time.time()
        dt = t1 - t0
        log.debug(f"{func.__name__!r} runtime: {dt * 1000.0:.2f} ms")
        return output

    return wrapped
