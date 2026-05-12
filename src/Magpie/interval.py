
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


iv_fn_mappings = {
    "sin": imath.sin,
    "cos": imath.cos,
    "log": imath.log,
    "sqrt": imath.sqrt,
    "exp": imath.exp,
}


def generate_bounds(X, moe=0.0):
    # moe = margin of error, eg 0.05 for a 5% margin of error
    bounds = []
    for lb, ub in zip(np.min(X, axis=0), np.max(X, axis=0)):
        bounds.append(interval([lb * (1 - moe), ub * (1 + moe)]))
    return np.array(bounds)
