"""Exact, opt-in sign test for one monomial projection of raw 12D sketches.

Only the Python standard library is needed. See MONOMIAL_SIGN_TEST.md for
the statistical model: the declarations below are assumptions, not checks
of ellipticity, independent selection, measurement error, or sampling.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Iterable, Sequence


REQUIRED_ASSUMPTIONS = frozenset({
    "clean_iid_elliptical_finite_covariance",
    "known_zero_centre",
    "pair_theta_fixed_without_clean_holdout",
    "fixed_n_single_test",
    "certified_intervals_for_raw_sketches",
})


def rational(value: int | str | Fraction) -> Fraction:
    """Reject binary floats and booleans rather than claiming exact measures."""
    if type(value) not in (int, str, Fraction):
        raise ValueError("Use integer or rational-string inputs, not floats/bools")
    return Fraction(value)


def binomial_upper_tail(n: int, k: int) -> Fraction:
    """P(Bin(n, 1/2) >= k), with exact integers and rational output."""
    if type(n) is not int or n < 0 or type(k) is not int:
        raise ValueError("Invalid binomial parameters")
    if k <= 0:
        return Fraction(1)
    if k > n:
        return Fraction(0)
    return Fraction(sum(comb(n, j) for j in range(k, n + 1)), 2**n)


def projection_pair(degree: int, direction: Sequence) -> tuple[list, list]:
    if type(degree) is not int or degree not in (0, 1, 2):
        raise ValueError("Choose one monomial degree: 0, 1, or 2")
    c = tuple(map(rational, direction))
    if len(c) != 2 or not any(c):
        raise ValueError("A nonzero two-channel direction is required")
    a, b = [Fraction(0)] * 12, [Fraction(0)] * 12
    for j in range(2):
        a[4 * degree + j] = c[j]
        b[4 * degree + 2 + j] = c[j]
    return a, b


def project_interval(row: Sequence, direction: Sequence) -> tuple[Fraction, Fraction]:
    if len(row) != 12:
        raise ValueError("Expected 12 raw sketch intervals, not quadratic features")
    lower = upper = Fraction(0)
    for endpoints, c in zip(row, direction):
        x, y = map(rational, endpoints)
        if x > y:
            raise ValueError("Reversed interval")
        lower += min(c * x, c * y)
        upper += max(c * x, c * y)
    return lower, upper


def square_bounds(lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    minimum = Fraction(0) if lower <= 0 <= upper else min(lower**2, upper**2)
    return minimum, max(lower**2, upper**2)


def evaluate(
    rows: Iterable, *, degree: int, direction: Sequence, theta,
    n: int, alpha="1/20", max_replaced: int = 0,
    enabled: bool = False, assumptions=None,
) -> dict:
    """Test one preselected monomial; keep every tie/ambiguous row in n.

    Each row supplies 12 [lower, upper] rational intervals guaranteed to
    contain the true raw sketch coordinates. No interval or model inference
    happens here. At most max_replaced arbitrary replacements are permitted
    only under a genuinely known bound and selection independent of the clean
    holdout. A declaration cannot establish those facts.
    """
    if enabled is not True or assumptions is None or set(assumptions) != REQUIRED_ASSUMPTIONS:
        raise ValueError("Explicit opt-in and all model declarations are required")
    if type(n) is not int or n <= 0:
        raise ValueError("Fix a positive sample size before observing the test sample")
    if type(max_replaced) is not int or not 0 <= max_replaced <= n:
        raise ValueError("Invalid known replacement bound")
    theta, alpha = rational(theta), rational(alpha)
    if not 0 < theta <= 1 or not 0 < alpha < 1:
        raise ValueError("Require 0 < theta <= 1 and 0 < alpha < 1")
    a, b = projection_pair(degree, direction)
    total = negatives = ambiguous = 0
    for row in rows:
        xlow, xhigh = square_bounds(*project_interval(row, a))
        ylow, yhigh = square_bounds(*project_interval(row, b))
        lower, upper = theta * xlow - yhigh, theta * xhigh - ylow
        negatives += upper < 0
        ambiguous += lower < 0 <= upper
        total += 1
    if total != n:
        raise ValueError("Observed rows differ from the fixed n; retain ties and ambiguous rows")
    corrected = max(0, negatives - max_replaced)
    p = binomial_upper_tail(n, corrected)
    return {
        "n": n, "degree": degree, "direction": [str(rational(c)) for c in direction],
        "theta": str(theta), "certified_negative": negatives, "ambiguous": ambiguous,
        "max_replaced": max_replaced, "corrected_negative": corrected,
        "p_exact": str(p), "alpha": str(alpha), "reject": p <= alpha,
        "model_verified_by_code": False,
        "nonrejection": "inconclusive; no support or positive-gap certificate",
        "scope": "one fixed monomial channel projection of raw 12D sketches",
    }


def _reject_json_float(value):
    raise ValueError("JSON decimals must be quoted rational strings with certified meaning")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON contract and enclosing rational intervals")
    parser.add_argument("--enable", action="store_true", help="explicitly enable this separate test")
    args = parser.parse_args(argv)
    try:
        raw = args.input.read_bytes()
        payload = json.loads(raw, parse_float=_reject_json_float, parse_constant=_reject_json_float)
        result = evaluate(enabled=args.enable, **payload)
    except (ValueError, TypeError, KeyError, ZeroDivisionError, OSError) as exc:
        parser.exit(2, f"Invalid input: {exc}\n")
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
