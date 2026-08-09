"""Numerically check the robust-localizer primal/dual minimax identity."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np
from scipy.stats import chi2

from pilot_covariance_sdp import (
    covariance_model,
    localizer_design,
    localizer_matrix,
    margins,
    moments,
)


def primal_value(
    observed: np.ndarray, covariance: np.ndarray, theta: float, degree: int, q: float
) -> tuple[float, str]:
    root = np.linalg.cholesky(covariance)
    design = localizer_design(theta, degree)[1:]
    z = cp.Variable(covariance.shape[0])
    level = cp.Variable()
    perturbation = root @ z
    matrix = localizer_matrix(observed, theta, degree)
    for k in range(design.shape[0]):
        matrix = matrix + perturbation[k] * design[k]
    problem = cp.Problem(
        cp.Maximize(level),
        [cp.norm(z, 2) <= q, matrix - level * np.eye(degree + 1) >> 0],
    )
    problem.solve(solver="CLARABEL", tol_gap_abs=1e-9, tol_feas=1e-9, max_iter=500)
    return float(problem.value), problem.status


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--degree", type=int, default=5)
    parser.add_argument("--sigma", type=float, default=1e-6)
    parser.add_argument("--correlation", type=float, default=0.92)
    parser.add_argument("--seed", type=int, default=26082026)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("pilot_primal_dual_equivalence.json"))
    args = parser.parse_args()
    theta, target, gamma, coverage = 0.8, 0.9, 0.01, 0.95
    low_atoms = np.array([0.15, 0.45, 0.78])
    low_weights = np.array([0.2, 0.3, 0.5])
    cases = {
        "boundary_null": moments(np.array([theta]), np.array([1.0]), 2 * args.degree + 2),
        "visible_alternative": moments(
            np.append(low_atoms, target),
            np.append((1.0 - gamma) * low_weights, gamma),
            2 * args.degree + 2,
        ),
    }
    covariance = covariance_model(2 * args.degree + 1, args.sigma, args.correlation)
    q = float(np.sqrt(chi2.ppf(coverage, covariance.shape[0])))
    rng = np.random.default_rng(args.seed)
    output: dict[str, object] = {"setup": vars(args) | {"output": str(args.output)}}
    for name, truth in cases.items():
        rows = []
        mismatched_signs = 0
        max_abs_difference = 0.0
        for _ in range(args.trials):
            observed = truth.copy()
            observed[1:] += rng.multivariate_normal(np.zeros(covariance.shape[0]), covariance)
            dual = margins(observed, covariance, theta, target, args.degree, coverage)
            primal, status = primal_value(observed, covariance, theta, args.degree, q)
            dual_value = float(dual["optimized_ellipsoid_upper"])
            mismatched_signs += (primal < -1e-7) != (dual_value < -1e-7)
            max_abs_difference = max(max_abs_difference, abs(primal - dual_value))
            if len(rows) < 5:
                rows.append({"primal": primal, "dual": dual_value, "status": status})
        output[name] = {
            "trials": args.trials,
            "mismatched_signs": mismatched_signs,
            "max_absolute_primal_dual_difference": max_abs_difference,
            "first_five": rows,
        }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
