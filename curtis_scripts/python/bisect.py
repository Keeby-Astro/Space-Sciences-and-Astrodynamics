# Compatibility for the original bisect module name.
# The root-finding algorithm is contained in bisect_root.py.

from _bisect import bisect_left, bisect_right, insort_left, insort_right
from bisect_root import bisect as bisect_root

insort = insort_right

def bisect(fun, xl, xu=None, *args, **kwargs):
    """Preserve root-finding calls and standard sequence-bisection imports."""
    if callable(fun):
        return bisect_root(fun, xl, xu)
    if xu is None:
        return bisect_right(fun, xl, *args, **kwargs)
    return bisect_right(fun, xl, xu, *args, **kwargs)