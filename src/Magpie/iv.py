import math
import numpy as np

"""
Interval arithmetic for GP singularity detection — pure Python/numpy,
no external C dependencies (replaces pyinterval + crlibm).

An interval [lo, hi] represents the set of all real values between lo and hi.
Arithmetic over intervals produces a result interval that contains the true
result for any point-value input within the input interval. We accept a small
floating-point rounding error in the bounds in exchange for not requiring crlibm.

The purpose here (following Keijzer 2003) is to detect dangerous expressions
that could produce inf or nan on the variable bounds we expect at runtime.
If evaluating the expression with interval inputs produces an interval
that is not finite, we raise SingularityException to reject the individual.

Return value contract
---------------------
| Case                               | Return             | Detected as singularity? |
|------------------------------------|--------------------|--------------------------|
| Valid, bounded                     | (lo, hi) finite    | No  (correct)            |
| Unbounded (e.g. log([0, 2]))        | (-inf, hi)         | Yes                      |
| Entirely out of domain (e.g. log([-2,-1])) | (nan, nan) | Yes                      |
| Denominator spans zero             | (-inf, +inf)       | Yes                      |
| Indeterminate (e.g. inf - inf)     | (nan, nan)         | Yes                      |

Note: fitness.py's check_intervals must use isinstance(result, Iv) rather than
the pyinterval interval class, and check math.isfinite(result.lo) and
math.isfinite(result.hi) (which already catches nan since isfinite(nan) is False).
"""

_nan = float('nan')
_inf = float('inf')


class Iv:
    """A closed interval [lo, hi] over the extended reals."""

    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi):
        self.lo = float(lo)
        self.hi = float(hi)

    def __repr__(self):
        return f'Iv[{self.lo}, {self.hi}]'

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __abs__(self):
        return interval_abs(self)

    def __add__(self, other):
        if not isinstance(other, Iv):
            other = Iv(other, other)
        if math.isnan(self.lo) or math.isnan(self.hi) or math.isnan(other.lo) or math.isnan(other.hi):
            return Iv(_nan, _nan)
        return Iv(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __sub__(self, other):
        if not isinstance(other, Iv):
            other = Iv(other, other)
        if math.isnan(self.lo) or math.isnan(self.hi) or math.isnan(other.lo) or math.isnan(other.hi):
            return Iv(_nan, _nan)
        return Iv(self.lo - other.hi, self.hi - other.lo)

    def __rsub__(self, other):
        return Iv(other, other) - self

    def __mul__(self, other):
        if not isinstance(other, Iv):
            other = Iv(other, other)
        a, b, c, d = self.lo, self.hi, other.lo, other.hi
        products = [a*c, a*d, b*c, b*d]
        if any(math.isnan(p) for p in products):
            return Iv(_nan, _nan)
        return Iv(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Iv):
            other = Iv(other, other)
        if math.isnan(other.lo) or math.isnan(other.hi):
            return Iv(_nan, _nan)
        if other.lo <= 0.0 <= other.hi:
            return Iv(-_inf, _inf)  # denominator spans zero
        return self * Iv(1.0 / other.hi, 1.0 / other.lo)

    def __rtruediv__(self, other):
        return Iv(other, other) / self

    def __pow__(self, other):
        if not isinstance(other, Iv):
            n = float(other)
            lo, hi = self.lo, self.hi
            if math.isnan(lo) or math.isnan(hi):
                return Iv(_nan, _nan)
            if n == 0.0:
                return Iv(1.0, 1.0)
            ni = int(n)
            if n == ni and ni > 0:
                if ni % 2 == 0:
                    if lo >= 0:
                        return Iv(lo**ni, hi**ni)
                    elif hi <= 0:
                        return Iv((-hi)**ni, (-lo)**ni)
                    else:
                        return Iv(0.0, max((-lo)**ni, hi**ni))
                else:
                    return Iv(lo**ni, hi**ni)  # odd: monotone
            # non-integer or negative exponent: require positive base
            if lo <= 0:
                return Iv(_nan, _nan)
            return Iv(lo**n, hi**n)
        # interval exponent: x^[a,b] via exp(b * log(x)), requires x > 0
        if math.isnan(self.lo) or math.isnan(other.lo):
            return Iv(_nan, _nan)
        if self.lo <= 0:
            return Iv(_nan, _nan)
        return iv_exp(other * iv_log(self))

    def __rpow__(self, base):
        return Iv(base, base).__pow__(self)


def _safe_exp(x):
    try:
        return math.exp(x)
    except OverflowError:
        return _inf


def _bounds(x):
    """Return (lo, hi) of an Iv or scalar."""
    if isinstance(x, Iv):
        return x.lo, x.hi
    v = float(x)
    return v, v


def iv_sin(x):
    if not isinstance(x, Iv):
        v = float(x)
        return Iv(math.sin(v), math.sin(v))
    lo, hi = x.lo, x.hi
    if math.isnan(lo) or math.isnan(hi):
        return Iv(_nan, _nan)
    if not math.isfinite(lo) or not math.isfinite(hi):
        return Iv(-1.0, 1.0)
    if hi - lo >= 2 * math.pi:
        return Iv(-1.0, 1.0)
    result_lo = min(math.sin(lo), math.sin(hi))
    result_hi = max(math.sin(lo), math.sin(hi))
    # check if a +1 peak (pi/2 + 2k*pi) falls inside [lo, hi]
    k = math.ceil((lo - math.pi / 2) / (2 * math.pi))
    if lo <= math.pi / 2 + k * 2 * math.pi <= hi:
        result_hi = 1.0
    # check if a -1 trough (3pi/2 + 2k*pi) falls inside [lo, hi]
    k = math.ceil((lo - 3 * math.pi / 2) / (2 * math.pi))
    if lo <= 3 * math.pi / 2 + k * 2 * math.pi <= hi:
        result_lo = -1.0
    return Iv(result_lo, result_hi)


def iv_cos(x):
    if not isinstance(x, Iv):
        v = float(x)
        return Iv(math.cos(v), math.cos(v))
    lo, hi = x.lo, x.hi
    if math.isnan(lo) or math.isnan(hi):
        return Iv(_nan, _nan)
    # cos(x) = sin(x + pi/2)
    return iv_sin(Iv(lo + math.pi / 2, hi + math.pi / 2))


def iv_exp(x):
    if not isinstance(x, Iv):
        v = float(x)
        return Iv(_safe_exp(v), _safe_exp(v))
    lo, hi = x.lo, x.hi
    if math.isnan(lo) or math.isnan(hi):
        return Iv(_nan, _nan)
    return Iv(_safe_exp(lo), _safe_exp(hi))


def iv_log(x):
    if not isinstance(x, Iv):
        v = float(x)
        if v < 0:
            return Iv(_nan, _nan)
        if v == 0:
            return Iv(-_inf, -_inf)
        return Iv(math.log(v), math.log(v))
    lo, hi = x.lo, x.hi
    if math.isnan(lo) or math.isnan(hi):
        return Iv(_nan, _nan)
    if hi < 0:
        return Iv(_nan, _nan)    # entirely out of domain
    if hi == 0:
        return Iv(-_inf, -_inf)  # log(0) = -inf
    hi_log = math.log(hi)
    if lo < 0:
        return Iv(_nan, _nan)    # domain contains negative values: undefined
    if lo == 0:
        return Iv(-_inf, hi_log) # natural limit: log(0+) → -inf
    return Iv(math.log(lo), hi_log)


def iv_sqrt(x):
    if not isinstance(x, Iv):
        v = float(x)
        if v < 0:
            return Iv(_nan, _nan)
        return Iv(math.sqrt(v), math.sqrt(v))
    lo, hi = x.lo, x.hi
    if math.isnan(lo) or math.isnan(hi):
        return Iv(_nan, _nan)
    if hi < 0:
        return Iv(_nan, _nan)  # entirely out of domain
    if lo < 0:
        return Iv(_nan, _nan)  # input range contains negatives: signal singularity
    return Iv(math.sqrt(lo), math.sqrt(hi))


def interval_abs(x):
    lo, hi = _bounds(x)
    if math.isnan(lo) or math.isnan(hi):
        return Iv(_nan, _nan)
    if lo >= 0:
        return Iv(lo, hi)
    elif hi <= 0:
        return Iv(-hi, -lo)
    else:
        return Iv(0.0, max(-lo, hi))


def interval_max(x, y):
    xl, xh = _bounds(x)
    yl, yh = _bounds(y)
    if math.isnan(xl) or math.isnan(xh) or math.isnan(yl) or math.isnan(yh):
        return Iv(_nan, _nan)
    return Iv(max(xl, yl), max(xh, yh))


def interval_min(x, y):
    xl, xh = _bounds(x)
    yl, yh = _bounds(y)
    if math.isnan(xl) or math.isnan(xh) or math.isnan(yl) or math.isnan(yh):
        return Iv(_nan, _nan)
    return Iv(min(xl, yl), min(xh, yh))


iv_fn_mappings = {
    "sin": iv_sin,
    "cos": iv_cos,
    "log": iv_log,
    "sqrt": iv_sqrt,
    "exp": iv_exp,
    "abs": interval_abs,
    "max": interval_max,
    "min": interval_min,
}


def generate_bounds(X, moe=0.0):
    # moe = margin of error, eg 0.05 for a 5% margin of error
    bounds = []
    for lb, ub in zip(np.min(X, axis=0), np.max(X, axis=0)):
        margin = moe * (ub - lb)
        bounds.append(Iv(lb - margin, ub + margin))
    return bounds  # plain list, not np.array: numpy would flatten Iv objects
