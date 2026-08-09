"""Finite-sample unknown-covariance localizer test from a Wishart Loewner band."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np

from pilot_annni_block_hotelling import (
    block_design,
    exact_half_moments,
    joint_sketch_covariance,
    parameter_index,
)


def localizer_from_joint(variable: cp.Expression, degree: int, channels: int, theta: float) -> cp.Expression:
    """Build theta H_N-G_N from diagonal time blocks of a joint Gram matrix."""
    blocks = [
        variable[k * channels : (k + 1) * channels, k * channels : (k + 1) * channels]
        for k in range(2 * degree + 2)
    ]
    rows = []
    for i in range(degree + 1):
        rows.append(
            [theta * blocks[i + j] - blocks[i + j + 1] for j in range(degree + 1)]
        )
    return cp.bmat(rows)


def wishart_band(sample_gram: np.ndarray, alpha: float, sample_count: int) -> tuple[float, float, float]:
    dimension = sample_gram.shape[0]
    eta = (np.sqrt(dimension) + np.sqrt(2.0 * np.log(2.0 / alpha))) / np.sqrt(sample_count)
    if eta >= 1.0:
        raise ValueError(f"Wishart band requires eta<1, got {eta:.3f}")
    return float(eta), float((1.0 + eta) ** -2), float((1.0 - eta) ** -2)


def compatible(sample_gram: np.ndarray, degree: int, channels: int, theta: float, alpha: float, sample_count: int) -> tuple[bool, str]:
    _, lower, upper = wishart_band(sample_gram, alpha, sample_count)
    dimension = sample_gram.shape[0]
    # Work in the empirical correlation scaling.  The time-Gram matrices are
    # extremely ill-conditioned; congruence scaling preserves every Loewner
    # inequality and avoids treating a solver ridge as statistical uncertainty.
    scales = np.sqrt(np.maximum(np.diag(sample_gram), np.finfo(float).tiny))
    inverse = np.diag(1.0 / scales)
    empirical = inverse @ sample_gram @ inverse
    empirical = (empirical + empirical.T) / 2
    q = cp.Variable((dimension, dimension), symmetric=True)
    k = np.diag(scales) @ q @ np.diag(scales)
    localizer = localizer_from_joint(k, degree, channels, theta)
    problem = cp.Problem(
        cp.Minimize(0.0),
        [q - lower * empirical >> 0, upper * empirical - q >> 0, localizer >> 0],
    )
    try:
        problem.solve(solver="CLARABEL", tol_gap_abs=1e-8, tol_feas=1e-8, max_iter=500)
    except cp.error.SolverError:
        problem.solve(solver="SCS", eps=1e-6, max_iters=100_000)
    if problem.status in {"infeasible", "infeasible_inaccurate"}:
        return False, problem.status
    if problem.status not in {"optimal", "optimal_inaccurate"}:
        raise RuntimeError(problem.status)
    return True, problem.status


def exact_localizer_min(exact_moments: np.ndarray, degree: int, channels: int, theta: float) -> float:
    maximum = 2 * degree + 1
    params = parameter_index(maximum, channels)
    flat = np.array([exact_moments[k, a, b] for k, a, b in params])
    design = block_design(theta, degree, channels)
    return float(np.linalg.eigvalsh(np.einsum("k,kij->ij", flat, design))[0])


def run(args: argparse.Namespace) -> dict[str, object]:
    channels = 2
    maximum = 2 * args.degree + 1
    half, metadata = exact_half_moments(args.length, args.degree, args.tau, args.j1, args.j2, args.hx)
    joint = joint_sketch_covariance(half, maximum)
    exact_moments = half[::2][: maximum + 1]
    true_edge = float(np.exp(-args.tau * metadata["true_gap"]))
    rng = np.random.default_rng(args.seed)
    rows = []
    for sample_count in args.sample_counts:
        eta, lower, upper = wishart_band(joint, args.alpha, sample_count)
        for factor in args.gap_factors:
            theta = float(np.exp(-args.tau * factor * metadata["true_gap"]))
            detections = 0
            statuses: dict[str, int] = {}
            for _ in range(args.ensembles):
                z = rng.multivariate_normal(np.zeros(joint.shape[0]), joint, size=sample_count)
                sample_gram = z.T @ z / sample_count
                is_compatible, status = compatible(
                    sample_gram, args.degree, channels, theta, args.alpha, sample_count
                )
                detections += not is_compatible
                statuses[status] = statuses.get(status, 0) + 1
            rows.append(
                {
                    "sample_count": sample_count,
                    "gap_factor": factor,
                    "eta": eta,
                    "loewner_lower": lower,
                    "loewner_upper": upper,
                    "exact_localizer_min": exact_localizer_min(exact_moments, args.degree, channels, theta),
                    "detection_rate": detections / args.ensembles,
                    "statuses": statuses,
                }
            )
    return {
        "schema_version": 1,
        "setup": vars(args) | {"output": str(args.output)},
        "model": metadata,
        "joint_dimension": int(joint.shape[0]),
        "true_edge": true_edge,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=12)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--sample-counts", nargs="+", type=int, default=[250, 500, 1000, 2000, 5000])
    parser.add_argument("--gap-factors", nargs="+", type=float, default=[1.0, 1.15, 1.35])
    parser.add_argument("--ensembles", type=int, default=100)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260813)
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--j1", type=float, default=1.0)
    parser.add_argument("--j2", type=float, default=0.37)
    parser.add_argument("--hx", type=float, default=2.2)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("pilot_wishart_loewner.json"))
    args = parser.parse_args()
    payload = run(args)
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(payload["rows"], indent=2))


if __name__ == "__main__":
    main()
