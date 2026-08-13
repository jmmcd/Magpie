import math
import numpy as np
import pytest
import SBGP.iv as m

Iv           = m.Iv
_bounds      = m._bounds
interval_abs = m.interval_abs
interval_max = m.interval_max
interval_min = m.interval_min
iv_fn_mappings = m.iv_fn_mappings
generate_bounds = m.generate_bounds
iv_sin  = m.iv_sin
iv_cos  = m.iv_cos
iv_log  = m.iv_log
iv_sqrt = m.iv_sqrt
iv_exp  = m.iv_exp


# --- helpers ---

def iv(lo, hi):
    return Iv(lo, hi)

def lo(x):
    return _bounds(x)[0]

def hi(x):
    return _bounds(x)[1]

def assert_contains(iv_result, value, tol=1e-12):
    l, h = _bounds(iv_result)
    assert l - tol <= value <= h + tol, \
        f"interval [{l}, {h}] does not contain {value}"

def assert_containment(iv_fn, scalar_fn, iv_input, n_samples=200, seed=42):
    l, h = _bounds(iv_input)
    result = iv_fn(iv_input)
    rng = np.random.default_rng(seed)
    for p in rng.uniform(l, h, n_samples):
        assert_contains(result, scalar_fn(p))

def is_nan_iv(x):
    l, h = _bounds(x)
    return math.isnan(l) and math.isnan(h)


# --- _bounds ---

class TestBounds:
    def test_scalar_positive(self):
        assert _bounds(3.5) == (3.5, 3.5)

    def test_scalar_negative(self):
        assert _bounds(-2.0) == (-2.0, -2.0)

    def test_simple_interval(self):
        l, h = _bounds(iv(1.0, 4.0))
        assert l == pytest.approx(1.0)
        assert h == pytest.approx(4.0)


# --- interval_abs ---

class TestIntervalAbs:
    def test_purely_positive(self):
        result = interval_abs(iv(2.0, 5.0))
        assert lo(result) == pytest.approx(2.0)
        assert hi(result) == pytest.approx(5.0)

    def test_purely_negative(self):
        result = interval_abs(iv(-5.0, -2.0))
        assert lo(result) == pytest.approx(2.0)
        assert hi(result) == pytest.approx(5.0)

    def test_spanning_zero_positive_side_larger(self):
        result = interval_abs(iv(-3.0, 4.0))
        assert lo(result) == pytest.approx(0.0)
        assert hi(result) == pytest.approx(4.0)

    def test_spanning_zero_negative_side_larger(self):
        result = interval_abs(iv(-5.0, 2.0))
        assert lo(result) == pytest.approx(0.0)
        assert hi(result) == pytest.approx(5.0)

    def test_containment(self):
        assert_containment(interval_abs, abs, iv(-3.0, 4.0))

    def test_nan_input(self):
        assert is_nan_iv(interval_abs(iv(float('nan'), 1.0)))


# --- interval_max ---

class TestIntervalMax:
    def test_second_dominates(self):
        result = interval_max(iv(1.0, 2.0), iv(3.0, 4.0))
        assert lo(result) == pytest.approx(3.0)
        assert hi(result) == pytest.approx(4.0)

    def test_first_dominates(self):
        result = interval_max(iv(5.0, 6.0), iv(1.0, 2.0))
        assert lo(result) == pytest.approx(5.0)
        assert hi(result) == pytest.approx(6.0)

    def test_overlapping(self):
        result = interval_max(iv(1.0, 4.0), iv(2.0, 5.0))
        assert lo(result) == pytest.approx(2.0)
        assert hi(result) == pytest.approx(5.0)

    def test_containment(self):
        rng = np.random.default_rng(42)
        result = interval_max(iv(1.0, 4.0), iv(2.0, 5.0))
        for px, py in zip(rng.uniform(1, 4, 100), rng.uniform(2, 5, 100)):
            assert_contains(result, max(px, py))

    def test_nan_propagates(self):
        assert is_nan_iv(interval_max(iv(float('nan'), 1.0), iv(1.0, 2.0)))


# --- interval_min ---

class TestIntervalMin:
    def test_second_dominates(self):
        result = interval_min(iv(3.0, 4.0), iv(1.0, 2.0))
        assert lo(result) == pytest.approx(1.0)
        assert hi(result) == pytest.approx(2.0)

    def test_first_dominates(self):
        result = interval_min(iv(1.0, 2.0), iv(3.0, 4.0))
        assert lo(result) == pytest.approx(1.0)
        assert hi(result) == pytest.approx(2.0)

    def test_overlapping(self):
        result = interval_min(iv(1.0, 4.0), iv(2.0, 5.0))
        assert lo(result) == pytest.approx(1.0)
        assert hi(result) == pytest.approx(4.0)

    def test_containment(self):
        rng = np.random.default_rng(42)
        result = interval_min(iv(1.0, 4.0), iv(2.0, 5.0))
        for px, py in zip(rng.uniform(1, 4, 100), rng.uniform(2, 5, 100)):
            assert_contains(result, min(px, py))


# --- iv_sin ---

class TestSin:
    def test_containment_positive(self):
        assert_containment(iv_sin, math.sin, iv(0.0, math.pi))

    def test_containment_spanning_zero(self):
        assert_containment(iv_sin, math.sin, iv(-math.pi / 2, math.pi / 2))

    def test_peak_included(self):
        result = iv_sin(iv(0.0, math.pi))
        assert hi(result) == pytest.approx(1.0, abs=1e-9)

    def test_trough_included(self):
        result = iv_sin(iv(math.pi, 2 * math.pi))
        assert lo(result) == pytest.approx(-1.0, abs=1e-9)

    def test_full_period_covers_range(self):
        result = iv_sin(iv(-math.pi, math.pi))
        assert lo(result) == pytest.approx(-1.0, abs=1e-9)
        assert hi(result) == pytest.approx(1.0, abs=1e-9)

    def test_narrow_interval_tighter_than_minus1_1(self):
        # sin on [pi/4, 3pi/4] should give [sqrt(2)/2, 1], not [-1, 1]
        result = iv_sin(iv(math.pi / 4, 3 * math.pi / 4))
        assert lo(result) == pytest.approx(math.sqrt(2) / 2, abs=1e-9)
        assert hi(result) == pytest.approx(1.0, abs=1e-9)

    def test_nan_input(self):
        assert is_nan_iv(iv_sin(iv(float('nan'), 1.0)))

    def test_infinite_input(self):
        result = iv_sin(iv(float('-inf'), float('inf')))
        assert lo(result) == pytest.approx(-1.0)
        assert hi(result) == pytest.approx(1.0)


# --- iv_cos ---

class TestCos:
    def test_containment(self):
        assert_containment(iv_cos, math.cos, iv(0.0, math.pi))

    def test_cos_at_zero(self):
        result = iv_cos(iv(0.0, 0.0))
        assert lo(result) == pytest.approx(1.0, abs=1e-10)
        assert hi(result) == pytest.approx(1.0, abs=1e-10)

    def test_full_period_covers_range(self):
        result = iv_cos(iv(0.0, 2 * math.pi))
        assert lo(result) == pytest.approx(-1.0, abs=1e-9)
        assert hi(result) == pytest.approx(1.0, abs=1e-9)

    def test_narrow_interval_not_minus1_1(self):
        # cos on [0, pi/2]: range is [0, 1]
        result = iv_cos(iv(0.0, math.pi / 2))
        assert lo(result) == pytest.approx(0.0, abs=1e-9)
        assert hi(result) == pytest.approx(1.0, abs=1e-9)

    def test_nan_input(self):
        assert is_nan_iv(iv_cos(iv(float('nan'), 1.0)))


# --- iv_log ---

class TestLog:
    def test_containment(self):
        assert_containment(iv_log, math.log, iv(1.0, math.e ** 2))

    def test_log_of_one_is_zero(self):
        result = iv_log(iv(1.0, 1.0))
        assert lo(result) == pytest.approx(0.0, abs=1e-10)
        assert hi(result) == pytest.approx(0.0, abs=1e-10)

    def test_entirely_negative_returns_nan(self):
        assert is_nan_iv(iv_log(iv(-2.0, -1.0)))

    def test_zero_upper_bound_returns_neg_inf(self):
        result = iv_log(iv(-1.0, 0.0))
        assert math.isinf(lo(result)) and lo(result) < 0
        assert math.isinf(hi(result)) and hi(result) < 0

    def test_mixed_input_returns_nan(self):
        # lo < 0 means domain violation → (nan, nan) under Option A
        assert is_nan_iv(iv_log(iv(-1.0, 2.0)))

    def test_zero_crossing_unbounded(self):
        result = iv_log(iv(0.0, 1.0))
        assert math.isinf(lo(result)) and lo(result) < 0

    def test_nan_input(self):
        assert is_nan_iv(iv_log(iv(float('nan'), 1.0)))


# --- iv_sqrt ---

class TestSqrt:
    def test_containment(self):
        assert_containment(iv_sqrt, math.sqrt, iv(1.0, 4.0))

    def test_sqrt_zero_to_four(self):
        result = iv_sqrt(iv(0.0, 4.0))
        assert lo(result) == pytest.approx(0.0, abs=1e-10)
        assert hi(result) == pytest.approx(2.0, abs=1e-10)

    def test_entirely_negative_returns_nan(self):
        assert is_nan_iv(iv_sqrt(iv(-1.0, -0.5)))

    def test_partially_negative_returns_nan(self):
        # conservative: any negative values in range → signal singularity
        assert is_nan_iv(iv_sqrt(iv(-1.0, 4.0)))

    def test_nan_input(self):
        assert is_nan_iv(iv_sqrt(iv(float('nan'), 1.0)))


# --- iv_exp ---

class TestExp:
    def test_containment(self):
        assert_containment(iv_exp, math.exp, iv(-1.0, 2.0))

    def test_exp_at_zero_is_one(self):
        result = iv_exp(iv(0.0, 0.0))
        assert lo(result) == pytest.approx(1.0, abs=1e-10)
        assert hi(result) == pytest.approx(1.0, abs=1e-10)

    def test_very_large_input_overflows(self):
        result = iv_exp(iv(1000.0, 2000.0))
        assert math.isinf(hi(result))

    def test_nan_input(self):
        assert is_nan_iv(iv_exp(iv(float('nan'), 1.0)))


# --- interval arithmetic ---

class TestArithmetic:
    def test_addition_containment(self):
        result = iv(1.0, 3.0) + iv(2.0, 4.0)
        rng = np.random.default_rng(42)
        for px, py in zip(rng.uniform(1, 3, 100), rng.uniform(2, 4, 100)):
            assert_contains(result, px + py)

    def test_subtraction_containment(self):
        result = iv(3.0, 5.0) - iv(1.0, 2.0)
        rng = np.random.default_rng(42)
        for px, py in zip(rng.uniform(3, 5, 100), rng.uniform(1, 2, 100)):
            assert_contains(result, px - py)

    def test_multiplication_containment(self):
        result = iv(-2.0, 3.0) * iv(-1.0, 4.0)
        rng = np.random.default_rng(42)
        for px, py in zip(rng.uniform(-2, 3, 100), rng.uniform(-1, 4, 100)):
            assert_contains(result, px * py)

    def test_division_by_zero_crossing_unbounded(self):
        result = iv(1.0, 2.0) / iv(-1.0, 1.0)
        l, h = _bounds(result)
        assert math.isinf(l) and l < 0
        assert math.isinf(h) and h > 0

    def test_nan_propagates_through_add(self):
        assert is_nan_iv(iv(float('nan'), 1.0) + iv(1.0, 2.0))

    def test_nan_propagates_through_mul(self):
        assert is_nan_iv(iv(float('nan'), 1.0) * iv(1.0, 2.0))

    def test_scalar_operands(self):
        # scalars should be treated as degenerate intervals
        result = iv(1.0, 3.0) + 2.0
        assert lo(result) == pytest.approx(3.0)
        assert hi(result) == pytest.approx(5.0)

    def test_even_power(self):
        # x^2 on [-2, 3] → [0, 9]
        result = iv(-2.0, 3.0) ** 2
        assert lo(result) == pytest.approx(0.0)
        assert hi(result) == pytest.approx(9.0)

    def test_odd_power(self):
        # x^3 on [-2, 3] → [-8, 27]
        result = iv(-2.0, 3.0) ** 3
        assert lo(result) == pytest.approx(-8.0)
        assert hi(result) == pytest.approx(27.0)


# --- generate_bounds ---

class TestGenerateBounds:
    def test_output_length_matches_n_vars(self):
        X = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        assert len(generate_bounds(X)) == 3

    def test_contains_column_extremes(self):
        X = np.array([[1.0, 10.0], [3.0, 20.0], [2.0, 15.0]])
        bounds = generate_bounds(X)
        assert_contains(bounds[0], 1.0)
        assert_contains(bounds[0], 3.0)
        assert_contains(bounds[1], 10.0)
        assert_contains(bounds[1], 20.0)

    def test_exact_bounds_no_moe(self):
        X = np.array([[1.0, 10.0], [3.0, 20.0]])
        bounds = generate_bounds(X, moe=0.0)
        assert lo(bounds[0]) == pytest.approx(1.0)
        assert hi(bounds[0]) == pytest.approx(3.0)
        assert lo(bounds[1]) == pytest.approx(10.0)
        assert hi(bounds[1]) == pytest.approx(20.0)

    def test_moe_expands_bounds(self):
        X = np.array([[0.0, 0.0], [10.0, 10.0]])
        bounds = generate_bounds(X, moe=0.1)
        assert lo(bounds[0]) == pytest.approx(-1.0)
        assert hi(bounds[0]) == pytest.approx(11.0)

    def test_single_variable(self):
        X = np.array([[2.0], [5.0], [3.0]])
        bounds = generate_bounds(X)
        assert len(bounds) == 1
        assert_contains(bounds[0], 2.0)
        assert_contains(bounds[0], 5.0)

    def test_returns_list_not_array(self):
        X = np.array([[0.0], [1.0]])
        assert isinstance(generate_bounds(X), list)

    def test_returns_iv_objects(self):
        X = np.array([[0.0], [1.0]])
        bounds = generate_bounds(X)
        assert isinstance(bounds[0], Iv)
