
from .exceptions import SingularityException
import math

"""

A simple and incomplete implementation of interval arithmetic.

Interval arithmetic is arithmetic over intervals, instead of over
single values.

An interval is a segment of the real line, [a, b], eg [0, 1].

For example, if x is in [0, 1] and y is in [-1, 0], then x + y is in
[-1, 1] (eg x = 0, y = -1 gives x + y = -1; and eg x = 1, y = 0 gives
x + y = 1).

x * y is in [-1, 0]

1 / x is in [1, inf]

(and so on).

The point of interval arithmetic in GP is to identify dangerous
expressions, which could give values such as inf, -inf, nan. We don't
want those in our programs. Keijzer (2003?) used this idea.

To use this idea, we need to have bounds on our variables x0, x1,
etc. For example, if x0 is measured in kg we know it will never be
negative. Another approach is just to guess the bounds based on the
training data.

"""


class Interval:
    def __init__(self, a, b):
        if (not math.isfinite(a)) or (not math.isfinite(b)):
            raise SingularityException
        if a > b:
            a, b = b, a       
        self.a, self.b = a, b
    def __repr__(self):
        return f"Interval([{self.a}, {self.b}])"
    def __neg__(self):
        return Interval(-self.b, -self.a)
    def __add__(self, c):
        if type(c) == type(self):
            return Interval(self.a + c.a, self.b + c.b)
        else:
            return Interval(self.a + c, self.b + c)
    def __radd__(self, c):
        return self + c
    def __sub__(self, c):
        if type(c) == type(self):
            return Interval(self.a - c.b, self.b - c.a)
        else:
            return Interval(self.a - c, self.b - c)
    def __rsub__(self, c):
        return c + -self
    def __mul__(self, c):
        if type(c) == type(self):
            # lazy
            vs = self.a * c.a, self.a * c.b, self.b * c.a, self.b * c.b
            return Interval(min(vs), max(vs))
        else:
            return Interval(self.a * c, self.b * c)
    def __rmul__(self, c):
        return self * c
    def __truediv__(self, c):
        if type(c) == type(self):
            if c.a <= 0 <= c.b: raise SingularityException
            # lazy
            vs = self.a / c.a, self.a / c.b, self.b / c.a, self.b / c.b
            return Interval(min(vs), max(vs))
        else:
            try:
                return Interval(self.a / c, self.b / c)
            except ZeroDivisionError:
                raise SingularityException
    def __rtruediv__(self, c):
        return Interval(c, c) / self
    def __pow__(self, c):
        if type(c) == type(self):
            if self.a < 0 and abs(c.a - c.b) > 0.000000001:
                # eg [-3, 3] ** [0, 1]: this includes singularity at -3 ** 0.5
                # (note np.power(-3, 0.5) -> math domain error, but -3 ** 0.5 -> ok)
                raise SingularityException
            try:
                # lazy
                vs = self.a ** c.a, self.a ** c.b, self.b ** c.a, self.b ** c.b
                return Interval(min(vs), max(vs))
            except (ValueError, ZeroDivisionError): # math domain error: pow(x, y) when x <0 and y non-integer
                # FIXME check this error
                #   File code/interval.py, line 91, in __pow__
                # vs = self.a ** c.a, self.a ** c.b, self.b ** c.a, self.b ** c.b
                # ZeroDivisionError: 0.0 cannot be raised to a negative power
                #
                raise SingularityException
        else:
            try:
                return Interval(self.a ** c, self.b ** c)
            except ValueError:
                raise SingularityException
    def __rpow__(self, c):
        return Interval(c, c) ** self

    # FIXME sin and cos are very conservative
    def sin(self):
        return Interval(-1, 1)
    def cos(self):
        return Interval(-1, 1)
    def log(self):
        if self.a <= 0 or self.b <= 0: raise SingularityException
        return Interval(math.log(self.a), math.log(self.b))
    def sqrt(self):
        if self.a < 0 or self.b < 0: raise SingularityException
        return Interval(math.sqrt(self.a), math.sqrt(self.b))
    def exp(self):
        return Interval(math.exp(self.a), math.exp(self.b))
    


## += (iadd) and other augmented operators are provided "for free" by python

## not done: floordiv, mod

def generate_bounds(X, moe=0.0):
    # moe = margin of error, eg 0.05 for a 5% margin of error
    import numpy as np
    bounds = []

    for lb, ub in zip(np.min(X, axis=1), np.max(X, axis=1)):

        bounds.append(Interval(lb * (1-moe), ub * (1+moe)))

    bounds = np.array(bounds)

    return bounds

