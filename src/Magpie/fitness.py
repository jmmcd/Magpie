import numpy as np
from scipy.optimize import curve_fit
from .exceptions import *
from .interval import iv_fn_mappings
# import sympy


def rmse(a, b): # RMSE. cost function to be minimised.
    return np.linalg.norm(a - b)

def one_m_r2(a, b): # 1-R^2. cost function to be minimised.
    if np.isscalar(a):
        return 1.0
    ss_res = np.sum((b - a) ** 2)
    ss_tot = np.sum((b - np.mean(b)) ** 2)
    if ss_tot == 0:
        return 0.0 if np.allclose(a, b) else 1.0
    return ss_res / ss_tot  # this is 1 - R²


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

def simplify(p, n_vars):
    X = sympy.symbols(f"X:{n_vars}")
    C = sympy.symbols(f"C:{n_vars}")
    s = p(X, C)
    zs = s.simplify()
    return str(zs)

def optimise(ps, X, y):
    """ Optimise the constants in ps to fit X, y.
    ps is a string, e.g. "C[0] * X[0]"
    X is a 2d numpy array, with shape (n_vars, n_samples)
    y is a 1d numpy array, with shape (n_samples,)
    Returns the optimal constants and the modified ps string with
    the constants renumbered to be C[0], C[1], ...
    """

    # eg ps = "((C[2] * X[2]) ** (X[0] * np.log(C[4])))"

    # count the number of constants in ps and remap the indices
    # to be 0, 1, 2, ... n_consts-1
    # eg C[2], C[4] -> C[0], C[1]
    replace_dict = {}
    
    # detect all constants using a regex on C[int]
    import re
    pattern = r'C\[(\d+)\]'
    matches = re.findall(pattern, ps)
    # make them unique and sort them
    unique_consts = sorted(set(int(m) for m in matches))
    n_consts = len(unique_consts)

    if n_consts == 0:
        return np.array([]), ps # no constants to optimise
    
    # create a replace dict
    for new_i, old_i in enumerate(unique_consts):
        replace_dict[f"C[{old_i}]"] = f"C[{new_i}]"
    
    ps = replaceall(ps, replace_dict)
    ps_np = replaceall(ps, {"sin": "np.sin", "cos": "np.cos",
                         "log": "np.log", "sqrt": "np.sqrt",
                         "exp": "np.exp", "abs": "np.abs",
                         "max": "np.maximum", "min": "np.minimum"})
    p = eval("lambda X, *C: " + ps_np) # we need the np-version for optimisation, but return the version without np. prefix
    
    # attempt to optimise constants
    # try several random starting points from normal dist
    # keep the best
    best_cost = np.inf
    best_popt = None
    for _ in range(5):
        try:
            C_init = np.random.randn(n_consts)
            popt, _ = curve_fit(p, X, y, p0=C_init)
            cost = np.sum((p(X, *popt) - y) ** 2)
            if cost < best_cost:
                best_cost = cost
                best_popt = popt
        except RuntimeError:
            continue
    if best_popt is None:
        raise FailedOptimisationError
    return best_popt, ps


def check_intervals(p, bounds):
    from interval import interval
    import math
    result = p(bounds)
    if isinstance(result, interval):
        if not all(math.isfinite(c[0]) and math.isfinite(c[1]) for c in result):
            raise SingularityException
    elif not math.isfinite(float(result)):
        raise SingularityException
    return result

def evaluate(ps, X_train, y_train, X_test=None, y_test=None, X_bounds=None):
    # X is in sklearn format, but we need it in transposed format,
    # because we will be using interval arithmetic and (later, Sympy
    # simplification) where it is natural to have a 1d array of
    # intervals. So, transpose here.
    X = X_train
    y = y_train
    X = X.T 
    n_vars = len(X)
    
    # p is a function of X and C.
    # X can be 2d numpy array, or an array of intervals, or Sympy vars

    newc, ps = optimise(ps, X, y)
    fn_mappings = {"sin": np.sin, "log": np.log, "cos": np.cos,
                   "exp": np.exp, "sqrt": np.sqrt, "abs": np.abs,
                   "max": np.maximum, "min": np.minimum}
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
        p_iv = eval("lambda X, C: " + ps, dict(iv_fn_mappings))
        newp_iv = lambda X: p_iv(X, newc)
        check_intervals(newp_iv, X_bounds)
    # p_simp = simplify(newp, n_vars) # not used - needs work TODO

    if X_test is not None:
        newpX_test = p_transpose(X_test)
        if np.isscalar(newpX_test):
            newpX_test = np.ones_like(y_test) * newpX_test
        cost_test = one_m_r2(newpX_test, y_test)
    else:
        cost_test = None
    
    # use a.atoms() to get a count of all symbols
    return cost, cost_test, ps, psc, p, newc, newp, p_transpose







