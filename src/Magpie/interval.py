
from interval import interval, imath
import numpy as np

"""

Interval arithmetic for GP singularity detection, using the pyinterval library.

An interval [a, b] represents the set of all real values between a and b.
Arithmetic over intervals produces a result interval that is guaranteed to
contain the true result for any point-value inputs within the input intervals.

The purpose here (following Keijzer 2003) is to detect dangerous expressions
that could produce inf or nan on the variable bounds we expect at runtime.
If evaluating the expression with interval inputs produces an interval
containing inf, we raise SingularityException to reject the individual.

"""


def _bounds(x):
    """Return (lo, hi) of a pyinterval interval or scalar."""
    if isinstance(x, interval):
        return min(c[0] for c in x), max(c[1] for c in x)
    v = float(x)
    return v, v

def interval_abs(x):
    lo, hi = _bounds(x)
    if lo >= 0:
        return interval([lo, hi])
    elif hi <= 0:
        return interval([-hi, -lo])
    else:
        return interval([0.0, max(-lo, hi)])

def interval_max(x, y):
    xl, xh = _bounds(x)
    yl, yh = _bounds(y)
    return interval([max(xl, yl), max(xh, yh)])

def interval_min(x, y):
    xl, xh = _bounds(x)
    yl, yh = _bounds(y)
    return interval([min(xl, yl), min(xh, yh)])

iv_fn_mappings = {
    "sin": imath.sin,
    "cos": imath.cos,
    "log": imath.log,
    "sqrt": imath.sqrt,
    "exp": imath.exp,
    "abs": interval_abs,
    "max": interval_max,
    "min": interval_min,
}


def generate_bounds(X, moe=0.0):
    # moe = margin of error, eg 0.05 for a 5% margin of error
    bounds = []
    for lb, ub in zip(np.min(X, axis=0), np.max(X, axis=0)):
        bounds.append(interval([lb * (1 - moe), ub * (1 + moe)]))
    return bounds  # plain list, not np.array: numpy would flatten interval objects into a float array
