"""Pilot for a strict noncommuting affine matrix-filter advantage.

The rational construction is the proof artifact.  The CVXPY solve is only an
optimality diagnostic; it is not used to certify the strict separation.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

import cvxpy as cp
import numpy as np


def rational_certificate() -> dict[str, object]:
    """Return exact interpolation, norm and commutator data."""

    a = (
        (Fraction(4, 5), Fraction(1, 5)),
        (Fraction(1, 5), Fraction(3, 10)),
    )
    b = (
        (Fraction(2, 5), Fraction(-2, 5)),
        (Fraction(-2, 5), Fraction(9, 10)),
    )

    p_half_e1 = (a[0][0] + b[0][0] / 2, a[1][0] + b[1][0] / 2)
    p_one_v2 = (
        a[0][0] + a[0][1] + b[0][0] + b[0][1],
        a[1][0] + a[1][1] + b[1][0] + b[1][1],
    )
    commutator_01 = (
        a[0][0] * b[0][1]
        + a[0][1] * b[1][1]
        - b[0][0] * a[0][1]
        - b[0][1] * a[1][1]
    )

    assert p_half_e1 == (Fraction(1), Fraction(0))
    assert p_one_v2 == (Fraction(1), Fraction(1))
    assert commutator_01 == Fraction(-1, 10)

    norm_at_zero = (11 + np.sqrt(41.0)) / 20
    norm_at_minus_one = (1 + np.sqrt(61.0)) / 10
    assert norm_at_minus_one < 1

    return {
        "A": [[str(value) for value in row] for row in a],
        "B": [[str(value) for value in row] for row in b],
        "P_half_e1": [str(value) for value in p_half_e1],
        "P_one_v2": [str(value) for value in p_one_v2],
        "commutator_01": str(commutator_01),
        "norm_at_zero": norm_at_zero,
        "norm_at_minus_one": norm_at_minus_one,
        "certified_upper_bound": norm_at_minus_one,
    }


def solve_grid(angle_degrees: float = 45.0, grid_size: int = 161) -> dict[str, object]:
    """Solve the affine minimax SDP on a verification grid."""

    phi = np.deg2rad(angle_degrees)
    v1 = np.array([1.0, 0.0])
    v2 = np.array([np.cos(phi), np.sin(phi)])
    a = cp.Variable((2, 2), symmetric=True)
    b = cp.Variable((2, 2), symmetric=True)
    bound = cp.Variable(nonneg=True)
    identity = np.eye(2)
    constraints = [(a + 0.5 * b) @ v1 == v1, (a + b) @ v2 == v2]
    for x in np.linspace(-1.0, 0.0, grid_size):
        value = a + x * b
        constraints.extend((value << bound * identity, -value << bound * identity))

    problem = cp.Problem(cp.Minimize(bound), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_gap_rel=1e-10,
        tol_feas=1e-10,
        max_iter=1000,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"unexpected solver status: {problem.status}")

    a_value = np.asarray(a.value)
    b_value = np.asarray(b.value)
    commutator = a_value @ b_value - b_value @ a_value
    return {
        "status": problem.status,
        "angle_degrees": angle_degrees,
        "grid_size": grid_size,
        "optimum": float(bound.value),
        "commutator_norm": float(np.linalg.norm(commutator, 2)),
        "A": a_value.tolist(),
        "B": b_value.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--angle", type=float, default=45.0)
    parser.add_argument("--grid-size", type=int, default=161)
    parser.add_argument("--skip-sdp", action="store_true")
    args = parser.parse_args()

    payload: dict[str, object] = {"rational_certificate": rational_certificate()}
    if not args.skip_sdp:
        payload["sdp_diagnostic"] = solve_grid(args.angle, args.grid_size)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
