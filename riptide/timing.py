import time
import logging
from functools import wraps


def timing(func, *args, **kwargs):
    @wraps(func)
    def wrapped(*args, **kwargs):
        log = logging.getLogger('riptide.timing')
        t0 = time.time()
        output = func(*args, **kwargs)
        t1 = time.time()
        dt = t1 - t0
        log.debug("{!r} runtime: {:.2f} ms".format(func.__name__, dt * 1000.0))
        return output
    return wrapped