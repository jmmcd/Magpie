import numpy as np
from scipy.optimize import curve_fit
from .exceptions import *
from .interval import Interval
# import sympy


def rmse(a, b): # RMSE. cost function to be minimised.
    return np.linalg.norm(a - b)

def one_m_r2(a, b): # 1-R^2. cost function to be minimised.
    if np.isscalar(a):
        return 1.0
    return 1 - np.corrcoef(a, b)[0, 1]

def replaceall(s, d):
    for k in d:
        s = s.replace(k, d[k])
    return s

def readable_eqn(ps, varnames):
    ps = replaceall(ps, {f'X[{i}]': varnames[i] for i in range(len(varnames))})
    return ps

def latex_eqn(ps, varnames):
    import sympy
    ps = readable_eqn(ps, varnames)
    # symbols = sympy.symbols(varnames)
    ps = sympy.sympify(ps)
    s = sympy.latex(ps)
    return s

def optimise(ps, X, y):
    ps = replaceall(ps, {"sin": "np.sin", "cos": "np.cos",
                         "log": "np.log", "sqrt": "np.sqrt"})
    p = eval("lambda X, *C: " + ps)
    C_init = np.ones(len(X)) # TODO should this be zeros?
    try:
        popt, pcov = curve_fit(p, X, y, p0=C_init)
    except RuntimeError:
        raise FailedOptimisationError
    return popt


def check_intervals(p, bounds):
    # calculate p(bounds). we can return the result but it doesn't
    # matter too much. the point is that if there are any singularities,
    # an exception will be raised during this
    return p(bounds)

def evaluate(ps, X_train, y_train, X_test=None, y_test=None, X_bounds=None):
    # X is in sklearn format, but we need it in transposed format,
    # because we will be using Interval and (later, Sympy
    # simplification) where it is natural to have a 1d array of
    # intervals. So, transpose here.
    X = X_train
    y = y_train
    X = X.T 
    n_vars = len(X)
    
    # p is a function of X and C.
    # X can be 2d numpy array, or an array of intervals, or Sympy vars

    newc = optimise(ps, X, y)
    fn_mappings = {"sin": np.sin, "log": np.log,
                   "cos": np.cos, "exp": np.exp} # TODO there could be more to add here
    p = eval("lambda X, C: " + ps, fn_mappings)
    newp = lambda X: p(X, newc)
    # newp: p with consts built in
    # cost = rmse(newp(X), y)

    replace_dict = {f"C[{i}]": f"{newc[i]:.4f}" for i in range(len(newc))}
    psc = replaceall(ps, replace_dict)

    p_transpose = lambda X: newp(X.T) # p_transpose accepts sklearn format X
    newpX = p_transpose(X_train)
    if np.isscalar(newpX):
        newpX = np.ones_like(y) * newpX
    cost = one_m_r2(newpX, y)
    # could raise a SingularityException
    if X_bounds is not None:
        check_intervals(newp, X_bounds)
    # p_simp = simplify(newp, n_vars)

    if X_test is not None:
        newpX_test = p_transpose(X_test)
        if np.isscalar(newpX_test):
            newpX_test = np.ones_like(y_test) * newpX_test
        cost_test = one_m_r2(newpX_test, y_test)
    else:
        cost_test = None
    
    # use a.atoms() to get a count of all symbols
    return cost, cost_test, ps, psc, p, newc, newp, p_transpose







