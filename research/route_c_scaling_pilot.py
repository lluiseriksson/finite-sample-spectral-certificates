"""Numerical stress test for the full-spark tangential-filter construction.

This is exploratory evidence only.  Full-spark checks and SDP output are
floating point; the exact theorem must use an explicit perturbation and a
certified interpolation bound.
"""

from __future__ import annotations

import argparse
import itertools
import json

import cvxpy as cp
import numpy as np


def target_data(
    degree: int,
    dimension: int,
    perturbation: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    if dimension < degree:
        raise ValueError("dimension must be at least the polynomial degree")
    number_targets = dimension + degree
    groups = list(range(dimension)) + list(range(degree))
    pass_points = np.linspace(0.4, 1.4, number_targets)
    rng = np.random.default_rng(seed)
    vectors = []
    for group in groups:
        vector = np.zeros(dimension)
        vector[group] = 1.0
        vector += perturbation * rng.normal(size=dimension)
        vector /= np.linalg.norm(vector)
        vectors.append(vector)
    return pass_points, np.asarray(vectors), groups


def minimum_full_spark_singular_value(vectors: np.ndarray) -> float:
    number_targets, dimension = vectors.shape
    minimum = np.inf
    for subset in itertools.combinations(range(number_targets), dimension):
        matrix = vectors[list(subset)].T
        minimum = min(minimum, np.linalg.svd(matrix, compute_uv=False)[-1])
    return float(minimum)


def solve(
    degree: int = 2,
    dimension: int = 4,
    perturbation: float = 1e-2,
    seed: int = 7,
    grid_size: int = 81,
    validation_grid_size: int = 4001,
) -> dict[str, object]:
    pass_points, vectors, groups = target_data(
        degree, dimension, perturbation, seed
    )
    full_spark_margin = minimum_full_spark_singular_value(vectors)

    # Match the theorem class: coefficients are real symmetric, but are not
    # constrained to commute with one another.
    coefficients = [
        cp.Variable((dimension, dimension), symmetric=True)
        for _ in range(degree + 1)
    ]
    bound = cp.Variable(nonneg=True)
    constraints = []
    for point, vector in zip(pass_points, vectors, strict=True):
        value = sum((point**power) * coefficients[power] for power in range(degree + 1))
        constraints.append(value @ vector == vector)

    identity = np.eye(dimension)
    for point in np.linspace(-1.0, 0.0, grid_size):
        value = sum((point**power) * coefficients[power] for power in range(degree + 1))
        constraints.append(
            cp.bmat([[bound * identity, value], [value.T, bound * identity]]) >> 0
        )

    problem = cp.Problem(cp.Minimize(bound), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-9,
        tol_gap_rel=1e-9,
        tol_feas=1e-9,
        max_iter=1000,
    )
    if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        raise RuntimeError(f"unexpected solver status: {problem.status}")

    values = [np.asarray(coefficient.value) for coefficient in coefficients]
    commutator_norm = max(
        np.linalg.norm(values[left] @ values[right] - values[right] @ values[left], 2)
        for left in range(degree + 1)
        for right in range(left)
    )
    interpolation_residual = max(
        np.linalg.norm(
            sum((point**power) * values[power] for power in range(degree + 1))
            @ vector
            - vector
        )
        for point, vector in zip(pass_points, vectors, strict=True)
    )
    validation_stopband_max = max(
        np.linalg.norm(
            sum((point**power) * values[power] for power in range(degree + 1)),
            2,
        )
        for point in np.linspace(-1.0, 0.0, validation_grid_size)
    )
    return {
        "degree": degree,
        "dimension": dimension,
        "number_targets": dimension + degree,
        "groups": groups,
        "perturbation": perturbation,
        "seed": seed,
        "grid_size": grid_size,
        "validation_grid_size": validation_grid_size,
        "coefficient_class": "real symmetric (not constrained to commute)",
        "status": problem.status,
        "full_spark_margin": full_spark_margin,
        "interpolation_residual": float(interpolation_residual),
        "stopband_optimum": float(bound.value),
        "validation_stopband_max": float(validation_stopband_max),
        "validation_grid_overshoot": float(validation_stopband_max - bound.value),
        "max_pairwise_commutator_norm": float(commutator_norm),
        "commuting_lower_bound_if_full_spark": 1.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--dimension", type=int, default=4)
    parser.add_argument("--perturbation", type=float, default=1e-2)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--grid-size", type=int, default=81)
    parser.add_argument("--validation-grid-size", type=int, default=4001)
    args = parser.parse_args()
    print(
        json.dumps(
            solve(
                degree=args.degree,
                dimension=args.dimension,
                perturbation=args.perturbation,
                seed=args.seed,
                grid_size=args.grid_size,
                validation_grid_size=args.validation_grid_size,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
