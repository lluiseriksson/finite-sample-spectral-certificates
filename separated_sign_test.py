"""Opt-in exact sign test for ONE certified separated pair of raw sketch blocks.

See SEPARATED_SIGN_TEST.md. Model declarations are assumptions, not checks.
The original monomial degree API remains unchanged.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from monomial_sign_test import (
    REQUIRED_ASSUMPTIONS, evaluate as evaluate_monomial,
    project_interval, rational, _reject_json_float,
)

CERTIFIED_PAIRS = frozenset({(0, 3), (2, 5), (0, 5)})


def projection_pair(time_pair, direction):
    """Return the two fixed raw linear projections, never a data-selected pair."""
    if (type(time_pair) not in (list, tuple) or len(time_pair) != 2
            or any(type(k) is not int for k in time_pair)
            or tuple(time_pair) not in CERTIFIED_PAIRS):
        raise ValueError("Choose exactly one certified pair: [0,3], [2,5], or [0,5]")
    if type(direction) not in (list, tuple) or len(direction) != 2:
        raise ValueError("A nonzero two-channel direction is required")
    v = tuple(map(rational, direction))
    if not any(v):
        raise ValueError("A nonzero two-channel direction is required")
    a, b = [Fraction(0)] * 12, [Fraction(0)] * 12
    for j in range(2):
        a[2 * time_pair[0] + j] = b[2 * time_pair[1] + j] = v[j]
    return a, b


def evaluate(rows, *, time_pair, direction, theta, n, alpha="1/20",
             max_replaced=0, enabled=False, assumptions=None):
    """Evaluate theta**(b-a) X^2-Y^2 with conservative rational intervals.

Exactly one pair, theta and channel vector must be fixed without the clean
holdout. Multiple calls do not acquire simultaneous coverage. All original
raw intervals are validated, including coordinates outside the projections.
    """
    a, b = projection_pair(time_pair, direction)
    theta = rational(theta)
    if not 0 < theta <= 1:
        raise ValueError("Require 0 < theta <= 1")
    weight = theta ** (time_pair[1] - time_pair[0])

    def projected_rows():
        for row in rows:
            x, y = project_interval(row, a), project_interval(row, b)
            # Reuse the existing conservative counting kernel, not its null
            # interpretation. The SOS proof supplies this pair's null mean.
            reduced = [[0, 0] for _ in range(12)]
            reduced[0], reduced[2] = x, y
            yield reduced

    result = evaluate_monomial(projected_rows(), degree=0, direction=[1, 0],
        theta=weight, n=n, alpha=alpha, max_replaced=max_replaced,
        enabled=enabled, assumptions=assumptions)
    result.pop("degree")
    result.update(time_pair=list(time_pair), theta=str(theta),
        direction=[str(rational(c)) for c in direction],
        quadratic_weight=str(weight),
        scope="one fixed certified separated-time pair of raw 12D sketches")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON contract and certified raw intervals")
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
