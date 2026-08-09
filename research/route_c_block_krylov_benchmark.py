"""Adversarial block-filter benchmark for the exact quadratic Route-C witness.

The experiment distinguishes two objectives that are often conflated:

1. suppress a stopband in a block subspace iteration; and
2. suppress it while preserving several prescribed target spectral rows.

The noncommuting quadratic certificate satisfies objective 2 exactly.  A
degree-two scalar Chebyshev filter is included as a deliberately strong
objective-1 baseline; it is normalized at only one target eigenvalue and is
therefore not feasible for the four tangential constraints.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


PASS_EIGENVALUES = np.array([0.5, 1.0, 1.5, 2.0])
PASS_ROWS = np.array([[1.0, -2.0], [1.0, -1.0], [1.0, 1.0], [1.0, 2.0]])

C0 = np.array([[2 / 5, -3 / 16], [-3 / 16, 29 / 32]])
C1 = np.array([[15 / 16, 3 / 20], [3 / 20, 3 / 32]])
C2 = np.array([[-3 / 8, 0.0], [0.0, -3 / 80]])


def matrix_polynomial(x: float) -> np.ndarray:
    return C0 + x * C1 + x * x * C2


def scalar_chebyshev(x: np.ndarray) -> np.ndarray:
    """Degree-two stopband Chebyshev filter normalized at lambda=1/2."""
    return (8 * x * x + 8 * x + 1) / 7


def orthonormal_columns(a: np.ndarray) -> np.ndarray:
    q, r = np.linalg.qr(a, mode="reduced")
    if np.min(np.abs(np.diag(r))) < 1e-12:
        raise RuntimeError("benchmark block lost numerical rank")
    return q


def largest_principal_sine(reference: np.ndarray, candidate: np.ndarray) -> float:
    u = orthonormal_columns(reference)
    v = orthonormal_columns(candidate)
    sigma_min = np.linalg.svd(u.T @ v, compute_uv=False)[-1]
    return float(np.sqrt(max(0.0, 1.0 - min(1.0, sigma_min) ** 2)))


def metrics(
    spectral_rows: np.ndarray,
    initial_target_rows: np.ndarray,
    clean_reference: np.ndarray,
    eigenvectors: np.ndarray,
) -> dict[str, float]:
    target = spectral_rows[:4]
    stop = spectral_rows[4:]
    block = eigenvectors @ spectral_rows
    best_scale = float(
        np.vdot(target, initial_target_rows).real / np.vdot(target, target).real
    )
    return {
        "pass_geometry_relative_error": float(
            np.linalg.norm(target - initial_target_rows)
            / np.linalg.norm(initial_target_rows)
        ),
        "best_global_rescaled_pass_error": float(
            np.linalg.norm(best_scale * target - initial_target_rows)
            / np.linalg.norm(initial_target_rows)
        ),
        "best_global_rescaling": best_scale,
        "stop_to_pass_frobenius_ratio": float(
            np.linalg.norm(stop) / np.linalg.norm(target)
        ),
        "desired_subspace_largest_sine": largest_principal_sine(
            clean_reference, block
        ),
    }


def run(n: int, seed: int, iterations: list[int]) -> dict[str, object]:
    if n < 8:
        raise ValueError("n must be at least 8")
    rng = np.random.default_rng(seed)

    stop_eigenvalues = np.linspace(-1.0, 0.0, n - 4)
    eigenvalues = np.concatenate([PASS_EIGENVALUES, stop_eigenvalues])

    raw = rng.standard_normal((n, n))
    eigenvectors, _ = np.linalg.qr(raw)

    stop_rows = rng.standard_normal((n - 4, 2))
    stop_rows *= 5 * np.linalg.norm(PASS_ROWS) / np.linalg.norm(stop_rows)
    initial = np.vstack([PASS_ROWS, stop_rows])
    clean_reference = eigenvectors[:, :4] @ PASS_ROWS

    matrix_rows = initial.copy()
    scalar_rows = initial.copy()
    commuting_rows = initial.copy()

    matrix_values = np.stack([matrix_polynomial(x) for x in eigenvalues])
    scalar_values = scalar_chebyshev(eigenvalues)

    requested = sorted(set(iterations))
    snapshots: dict[str, dict[str, dict[str, float]]] = {}
    for step in range(max(requested) + 1):
        if step in requested:
            snapshots[str(step)] = {
                "noncommuting_exact_tangential": metrics(
                    matrix_rows, PASS_ROWS, clean_reference, eigenvectors
                ),
                "scalar_chebyshev_one_point_normalized": metrics(
                    scalar_rows, PASS_ROWS, clean_reference, eigenvectors
                ),
                "commuting_exact_tangential": metrics(
                    commuting_rows, PASS_ROWS, clean_reference, eigenvectors
                ),
            }
        if step == max(requested):
            break
        matrix_rows = np.einsum("ni,nij->nj", matrix_rows, matrix_values)
        scalar_rows *= scalar_values[:, None]
        # The exact-feasible commuting quadratic is identically the identity.

    target_matrix_residual = max(
        np.linalg.norm(PASS_ROWS[j] @ matrix_values[j] - PASS_ROWS[j])
        for j in range(4)
    )
    stop_grid = np.linspace(-1.0, 0.0, 10001)
    empirical_matrix_stop_norm = max(
        np.linalg.norm(matrix_polynomial(x), 2) for x in stop_grid
    )
    empirical_scalar_stop_max = float(np.max(np.abs(scalar_chebyshev(stop_grid))))

    return {
        "experiment": "route_c_block_krylov_adversarial_benchmark",
        "seed": seed,
        "dimension": n,
        "block_width": 2,
        "degree_per_iteration": 2,
        "initial_stop_to_pass_ratio": 5.0,
        "exact_certificate_bound": 301 / 304,
        "target_matrix_residual": float(target_matrix_residual),
        "empirical_matrix_stop_norm_10001_grid": float(empirical_matrix_stop_norm),
        "empirical_scalar_stop_max_10001_grid": empirical_scalar_stop_max,
        "scalar_pass_multipliers": scalar_chebyshev(PASS_EIGENVALUES).tolist(),
        "interpretation": {
            "fair_exact_constraint_comparator": (
                "The commuting degree-two comparator is the identity because four "
                "full-spark tangential constraints force it exactly."
            ),
            "adversarial_unconstrained_comparator": (
                "The scalar Chebyshev filter suppresses the stopband much faster, "
                "but violates three of the four pass constraints and changes the "
                "desired clean block subspace."
            ),
        },
        "snapshots": snapshots,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--iterations", type=int, nargs="+", default=[0, 1, 5, 20, 100])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/route_c/block_krylov_benchmark.json"),
    )
    args = parser.parse_args()
    if min(args.iterations) < 0:
        raise ValueError("iterations must be nonnegative")
    result = run(args.dimension, args.seed, args.iterations)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
