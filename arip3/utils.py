import functools
from time import time
from .typing import *

round0 = lambda x, prec=3: x.round(prec)
round1 = lambda x, prec=3: round(x, prec)
rnd = round0


def timer(disable_print=False):
    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            s = time()
            r = fn(*args, **kwargs)
            t = time()
            if not disable_print:
                print(f'[Timer]: {fn.__name__} took {t - s:.3f}s')
            return r
        return wrapper
    return decorator


def timer_s(fn:Callable[[Any], Tuple[float, ...]]):
    def wrapper(*args, **kwargs):
        s = time()
        r = fn(*args, **kwargs)
        t = time()
        if r: return (t - s,) + tuple(r)
        else: return (t - s)
    return wrapper


def die(info:str, errcode:int=-1):
    print(info)
    exit(errcode)
