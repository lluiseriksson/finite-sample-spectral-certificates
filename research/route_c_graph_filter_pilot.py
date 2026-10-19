"""Natural-application pilot for tangential MIMO polynomial graph filters.

The graph filter has the standard multi-feature form

    Y = sum_{k=0}^N S^k X C_k,

where S is a symmetric normalized adjacency, X has two node features, and the
2x2 taps C_k mix the features.  Four low-frequency spectral rows of two smooth
fields are preserved exactly while the negative graph-frequency band is
attenuated.  The experiment compares symmetric noncommuting taps with a fully
general MIMO filter, a scalar Chebyshev filter, and the exact-feasible commuting
degree-two class (which is the identity by the full-spark root count).

This is an exploratory SDP plus dense-grid validation, not an exact certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np


def rectangular_grid_adjacency(rows: int, cols: int) -> np.ndarray:
    n = rows * cols
    adjacency = np.zeros((n, n))

    def index(r: int, c: int) -> int:
        return r * cols + c

    for r in range(rows):
        for c in range(cols):
            i = index(r, c)
            if r + 1 < rows:
                j = index(r + 1, c)
                adjacency[i, j] = adjacency[j, i] = 1.0
            if c + 1 < cols:
                j = index(r, c + 1)
                adjacency[i, j] = adjacency[j, i] = 1.0
    degree = adjacency.sum(axis=1)
    inv_sqrt = np.diag(1.0 / np.sqrt(degree))
    return inv_sqrt @ adjacency @ inv_sqrt


def smooth_node_features(rows: int, cols: int) -> np.ndarray:
    yy, xx = np.meshgrid(
        np.linspace(-1.0, 1.0, rows),
        np.linspace(-1.0, 1.0, cols),
        indexing="ij",
    )
    bump_1 = np.exp(-6.0 * ((xx - 0.35) ** 2 + (yy + 0.20) ** 2))
    bump_2 = np.exp(-7.0 * ((xx + 0.25) ** 2 + (yy - 0.30) ** 2))
    feature_1 = 1.0 + 0.35 * xx - 0.15 * yy + 0.40 * bump_1
    feature_2 = 0.8 - 0.20 * xx + 0.40 * yy + 0.35 * bump_2
    return np.column_stack([feature_1.ravel(), feature_2.ravel()])


def polynomial_value(coefficients: list[np.ndarray], x: float) -> np.ndarray:
    return sum((x**k) * coefficient for k, coefficient in enumerate(coefficients))


def solve_mimo_filter(
    pass_eigenvalues: np.ndarray,
    pass_directions: np.ndarray,
    stop_grid: np.ndarray,
    symmetric: bool,
) -> tuple[list[np.ndarray], float, str]:
    if symmetric:
        coefficients = [cp.Variable((2, 2), symmetric=True) for _ in range(3)]
    else:
        coefficients = [cp.Variable((2, 2)) for _ in range(3)]
    leakage = cp.Variable(nonneg=True)
    constraints: list[cp.Constraint] = []
    for eigenvalue, direction in zip(pass_eigenvalues, pass_directions):
        value = sum(
            (float(eigenvalue) ** k) * coefficient
            for k, coefficient in enumerate(coefficients)
        )
        constraints.append(value.T @ direction == direction)
    identity = np.eye(2)
    for x in stop_grid:
        value = sum((float(x) ** k) * coefficient for k, coefficient in enumerate(coefficients))
        constraints.append(
            cp.bmat([[leakage * identity, value], [value.T, leakage * identity]])
            >> 0
        )
    problem = cp.Problem(cp.Minimize(leakage), constraints)
    problem.solve(
        solver=cp.CLARABEL,
        tol_gap_abs=1e-9,
        tol_feas=1e-9,
        tol_gap_rel=1e-9,
        max_iter=500,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"MIMO graph-filter SDP failed: {problem.status}")
    numeric = [np.asarray(coefficient.value) for coefficient in coefficients]
    return numeric, float(problem.value), problem.status


def scalar_chebyshev(x: np.ndarray, normalization_node: float) -> np.ndarray:
    numerator = 2.0 * (2.0 * x + 1.0) ** 2 - 1.0
    denominator = 2.0 * (2.0 * normalization_node + 1.0) ** 2 - 1.0
    return numerator / denominator


def orthonormal_columns(matrix: np.ndarray) -> np.ndarray:
    q, r = np.linalg.qr(matrix, mode="reduced")
    if np.min(np.abs(np.diag(r))) < 1e-12:
        raise RuntimeError("graph benchmark block lost numerical rank")
    return q


def largest_principal_sine(reference: np.ndarray, candidate: np.ndarray) -> float:
    u = orthonormal_columns(reference)
    v = orthonormal_columns(candidate)
    sigma_min = np.linalg.svd(u.T @ v, compute_uv=False)[-1]
    return float(np.sqrt(max(0.0, 1.0 - min(1.0, sigma_min) ** 2)))


def block_metrics(
    spectral_rows: np.ndarray,
    pass_indices: np.ndarray,
    stop_indices: np.ndarray,
    initial_pass_rows: np.ndarray,
    eigenvectors: np.ndarray,
    clean_reference: np.ndarray,
) -> dict[str, float]:
    target = spectral_rows[pass_indices]
    stop = spectral_rows[stop_indices]
    best_scale = float(
        np.vdot(target, initial_pass_rows).real / np.vdot(target, target).real
    )
    return {
        "pass_geometry_relative_error": float(
            np.linalg.norm(target - initial_pass_rows) / np.linalg.norm(initial_pass_rows)
        ),
        "best_global_rescaled_pass_error": float(
            np.linalg.norm(best_scale * target - initial_pass_rows)
            / np.linalg.norm(initial_pass_rows)
        ),
        "stop_to_pass_frobenius_ratio": float(
            np.linalg.norm(stop) / np.linalg.norm(target)
        ),
        "desired_subspace_largest_sine": largest_principal_sine(
            clean_reference, eigenvectors @ spectral_rows
        ),
    }


def run(
    rows: int,
    cols: int,
    seed: int,
    optimization_points: int,
    validation_points: int,
    iterations: list[int],
) -> dict[str, object]:
    shift = rectangular_grid_adjacency(rows, cols)
    eigenvalues, eigenvectors = np.linalg.eigh(shift)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    pass_indices = np.arange(4)
    pass_eigenvalues = eigenvalues[pass_indices]
    if np.min(np.diff(pass_eigenvalues[::-1])) < 1e-8:
        raise RuntimeError("top four graph eigenvalues are not distinct")
    stop_indices = np.flatnonzero(eigenvalues <= 0.0)

    features = smooth_node_features(rows, cols)
    natural_rows = eigenvectors.T @ features
    initial_pass_rows = natural_rows[pass_indices]
    pass_directions = initial_pass_rows / np.linalg.norm(
        initial_pass_rows, axis=1, keepdims=True
    )
    minors = [
        abs(np.linalg.det(pass_directions[[i, j]]))
        for i in range(4)
        for j in range(i + 1, 4)
    ]
    if min(minors) < 1e-6:
        raise RuntimeError("natural pass signatures are not numerically full spark")

    optimization_grid = np.linspace(-1.0, 0.0, optimization_points)
    validation_grid = np.linspace(-1.0, 0.0, validation_points)
    symmetric_coefficients, symmetric_objective, symmetric_status = solve_mimo_filter(
        pass_eigenvalues, pass_directions, optimization_grid, symmetric=True
    )
    general_coefficients, general_objective, general_status = solve_mimo_filter(
        pass_eigenvalues, pass_directions, optimization_grid, symmetric=False
    )

    symmetric_validation = max(
        np.linalg.norm(polynomial_value(symmetric_coefficients, x), 2)
        for x in validation_grid
    )
    general_validation = max(
        np.linalg.norm(polynomial_value(general_coefficients, x), 2)
        for x in validation_grid
    )
    symmetric_residual = max(
        np.linalg.norm(
            polynomial_value(symmetric_coefficients, eigenvalue).T @ direction
            - direction
        )
        for eigenvalue, direction in zip(pass_eigenvalues, pass_directions)
    )
    general_residual = max(
        np.linalg.norm(
            polynomial_value(general_coefficients, eigenvalue).T @ direction
            - direction
        )
        for eigenvalue, direction in zip(pass_eigenvalues, pass_directions)
    )
    commutators = [
        np.linalg.norm(
            symmetric_coefficients[i] @ symmetric_coefficients[j]
            - symmetric_coefficients[j] @ symmetric_coefficients[i],
            2,
        )
        for i in range(3)
        for j in range(i + 1, 3)
    ]

    rng = np.random.default_rng(seed)
    initial_rows = np.zeros_like(natural_rows)
    initial_rows[pass_indices] = initial_pass_rows
    stop_noise = rng.standard_normal((len(stop_indices), 2))
    stop_noise *= 5.0 * np.linalg.norm(initial_pass_rows) / np.linalg.norm(stop_noise)
    initial_rows[stop_indices] = stop_noise
    clean_reference = eigenvectors[:, pass_indices] @ initial_pass_rows

    symmetric_rows = initial_rows.copy()
    general_rows = initial_rows.copy()
    scalar_rows = initial_rows.copy()
    commuting_rows = initial_rows.copy()
    symmetric_values = np.stack(
        [polynomial_value(symmetric_coefficients, x) for x in eigenvalues]
    )
    general_values = np.stack(
        [polynomial_value(general_coefficients, x) for x in eigenvalues]
    )
    scalar_values = scalar_chebyshev(eigenvalues, pass_eigenvalues[0])

    requested = sorted(set(iterations))
    snapshots: dict[str, dict[str, dict[str, float]]] = {}
    for step in range(max(requested) + 1):
        if step in requested:
            snapshots[str(step)] = {
                "symmetric_noncommuting_mimo": block_metrics(
                    symmetric_rows,
                    pass_indices,
                    stop_indices,
                    initial_pass_rows,
                    eigenvectors,
                    clean_reference,
                ),
                "general_mimo": block_metrics(
                    general_rows,
                    pass_indices,
                    stop_indices,
                    initial_pass_rows,
                    eigenvectors,
                    clean_reference,
                ),
                "scalar_chebyshev_one_point_normalized": block_metrics(
                    scalar_rows,
                    pass_indices,
                    stop_indices,
                    initial_pass_rows,
                    eigenvectors,
                    clean_reference,
                ),
                "commuting_exact_tangential": block_metrics(
                    commuting_rows,
                    pass_indices,
                    stop_indices,
                    initial_pass_rows,
                    eigenvectors,
                    clean_reference,
                ),
            }
        if step == max(requested):
            break
        symmetric_rows = np.einsum("ni,nij->nj", symmetric_rows, symmetric_values)
        general_rows = np.einsum("ni,nij->nj", general_rows, general_values)
        scalar_rows *= scalar_values[:, None]

    return {
        "experiment": "route_c_natural_mimo_graph_filter_pilot",
        "graph": {
            "type": "rectangular_4_neighbor_grid",
            "rows": rows,
            "cols": cols,
            "nodes": rows * cols,
            "shift": "symmetric_normalized_adjacency",
        },
        "filter_degree": 2,
        "feature_channels": 2,
        "pass_eigenvalues": pass_eigenvalues.tolist(),
        "minimum_full_spark_minor": float(min(minors)),
        "stop_eigenvalue_count": int(len(stop_indices)),
        "optimization_points": optimization_points,
        "validation_points": validation_points,
        "symmetric_mimo": {
            "status": symmetric_status,
            "objective": symmetric_objective,
            "validation_maximum": float(symmetric_validation),
            "maximum_pass_residual": float(symmetric_residual),
            "maximum_commutator_norm": float(max(commutators)),
            "coefficients": [coefficient.tolist() for coefficient in symmetric_coefficients],
        },
        "general_mimo": {
            "status": general_status,
            "objective": general_objective,
            "validation_maximum": float(general_validation),
            "maximum_pass_residual": float(general_residual),
            "coefficients": [coefficient.tolist() for coefficient in general_coefficients],
        },
        "scalar_chebyshev": {
            "validation_maximum": float(
                np.max(np.abs(scalar_chebyshev(validation_grid, pass_eigenvalues[0])))
            ),
            "pass_multipliers": scalar_chebyshev(
                pass_eigenvalues, pass_eigenvalues[0]
            ).tolist(),
        },
        "commuting_exact_tangential": {
            "leakage": 1.0,
            "reason": (
                "For d=N=2 and four distinct full-spark pass signatures, every "
                "pairwise-commuting symmetric feasible polynomial is I_2."
            ),
        },
        "seed": seed,
        "initial_stop_to_pass_ratio": 5.0,
        "snapshots": snapshots,
        "limitations": [
            "The MIMO leakage bounds use a finite optimization grid and dense validation, not exact arithmetic.",
            "The grid graph and smooth fields are natural signal-processing objects, but the four exact pass constraints are task-designed.",
            "No supervised learning or downstream classification claim is made.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=18)
    parser.add_argument("--cols", type=int, default=23)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--optimization-points", type=int, default=121)
    parser.add_argument("--validation-points", type=int, default=20001)
    parser.add_argument("--iterations", type=int, nargs="+", default=[0, 1, 2, 5, 20])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/route_c/graph_filter_pilot.json"),
    )
    args = parser.parse_args()
    result = run(
        args.rows,
        args.cols,
        args.seed,
        args.optimization_points,
        args.validation_points,
        args.iterations,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    assert result["minimum_full_spark_minor"] > 0.02
    assert 0.99 < result["symmetric_mimo"]["validation_maximum"] < 1.01
    assert result["general_mimo"]["validation_maximum"] < 0.61
    assert result["scalar_chebyshev"]["validation_maximum"] < 0.06
    scalar_two = result["snapshots"]["2"]["scalar_chebyshev_one_point_normalized"]
    assert scalar_two["best_global_rescaled_pass_error"] < 0.006
    assert scalar_two["desired_subspace_largest_sine"] < 0.04
    print(
        "PASS: natural graph pilot remains a negative result "
        f"(symmetric={result['symmetric_mimo']['validation_maximum']:.8f}, "
        f"scalar={result['scalar_chebyshev']['validation_maximum']:.8f})"
    )


if __name__ == "__main__":
    main()
