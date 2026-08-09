#!/usr/bin/env python3
"""Offline replay of the VBL-VA001 measured-calibration application gate."""

import json
import sys
from itertools import combinations
from pathlib import Path

import cvxpy as cp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verification.verify_vbl_va001_witness import build_coefficients  # noqa: E402

INPUT = ROOT / "data" / "route_c" / "vbl_va001_calibration.json"
OUTPUT = ROOT / "results" / "route_c" / "vbl_va001_calibration.json"
STOP_NODES = np.linspace(-1.0, -0.929, 151)
VALIDATION_NODES = np.linspace(-1.0, -0.929, 20_001)


def matrix_value(coefficients, node):
    return sum(node**degree * coefficients[degree] for degree in range(3))


def pass_residuals(coefficients, nodes, directions):
    return [
        float(np.linalg.norm(matrix_value(coefficients, node) @ direction - direction))
        for node, direction in zip(nodes, directions, strict=True)
    ]


def leading_direction(matrix):
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    direction = eigenvectors[:, -1]
    if direction[np.argmax(np.abs(direction))] < 0:
        direction = -direction
    return eigenvalues, direction


def acute_angle_degrees(left, right):
    return float(np.degrees(np.arccos(np.clip(abs(left @ right), 0, 1))))


def fdd_invariance_metrics(coefficients, nodes, training, held_out_cospectra):
    rows = []
    for node, target, cospectrum in zip(nodes, training, held_out_cospectra, strict=True):
        eigenvalues, original = leading_direction(cospectrum)
        response = matrix_value(coefficients, node)
        filtered_cospectrum = response @ cospectrum @ response.T
        filtered_eigenvalues, filtered = leading_direction(filtered_cospectrum)
        rows.append(
            {
                "unfiltered_to_filtered_angle_degrees": acute_angle_degrees(original, filtered),
                "training_to_filtered_angle_degrees": acute_angle_degrees(target, filtered),
                "leading_eigenvalue_ratio": float(filtered_eigenvalues[-1] / eigenvalues[-1]),
                "rayleigh_ratio_along_unfiltered_direction": float(
                    original @ filtered_cospectrum @ original / (original @ cospectrum @ original)
                ),
            }
        )
    return {
        "per_frequency": rows,
        "maximum_unfiltered_to_filtered_angle_degrees": max(
            row["unfiltered_to_filtered_angle_degrees"] for row in rows
        ),
        "maximum_absolute_leading_eigenvalue_ratio_error": max(
            abs(row["leading_eigenvalue_ratio"] - 1) for row in rows
        ),
    }


def solve_matrix(nodes, directions, symmetric):
    dimension = directions.shape[1]
    coefficients = [cp.Variable((dimension, dimension), symmetric=symmetric) for _ in range(3)]
    objective = cp.Variable(nonneg=True)
    constraints = []
    for node, direction in zip(nodes, directions, strict=True):
        value = sum(node**degree * coefficients[degree] for degree in range(3))
        constraints.append(value @ direction == direction)
    for node in STOP_NODES:
        constraints.append(cp.norm(sum(node**degree * coefficients[degree] for degree in range(3)), 2) <= objective)
    problem = cp.Problem(cp.Minimize(objective), constraints)
    problem.solve(solver="CLARABEL")
    values = np.array([coefficient.value for coefficient in coefficients])
    validation = max(np.linalg.norm(matrix_value(values, node), 2) for node in VALIDATION_NODES)
    commutator = max(
        np.linalg.norm(values[left] @ values[right] - values[right] @ values[left], 2)
        for left, right in combinations(range(3), 2)
    )
    return {
        "status": problem.status,
        "objective": float(problem.value),
        "validation_leakage": float(validation),
        "coefficients": values,
        "coefficient_frobenius_norms": [float(np.linalg.norm(value)) for value in values],
        "maximum_commutator_norm": float(commutator),
    }


def solve_scalar(nodes, tolerance):
    coefficients = cp.Variable(3)
    objective = cp.Variable(nonneg=True)
    constraints = [
        cp.abs(sum(node**degree * coefficients[degree] for degree in range(3)) - 1) <= tolerance
        for node in nodes
    ]
    constraints.extend(
        cp.abs(sum(node**degree * coefficients[degree] for degree in range(3))) <= objective
        for node in STOP_NODES
    )
    problem = cp.Problem(cp.Minimize(objective), constraints)
    problem.solve(solver="CLARABEL")
    values = coefficients.value
    validation = max(abs(sum(node**degree * values[degree] for degree in range(3))) for node in VALIDATION_NODES)
    residual = max(abs(sum(node**degree * values[degree] for degree in range(3)) - 1) for node in nodes)
    return {
        "status": problem.status,
        "objective": float(problem.value),
        "validation_leakage": float(validation),
        "maximum_training_residual": float(residual),
        "coefficients": values.tolist(),
    }


def solve_fixed_basis(nodes, directions, tolerance, basis):
    dimension = directions.shape[1]
    coefficients = cp.Variable((dimension, 3))
    objective = cp.Variable(nonneg=True)
    constraints = []
    for node, direction in zip(nodes, directions, strict=True):
        multipliers = coefficients @ np.array([1.0, node, node * node])
        coordinates = basis.T @ direction
        constraints.append(cp.norm(cp.multiply(multipliers - 1, coordinates), 2) <= tolerance)
    for node in STOP_NODES:
        constraints.append(cp.abs(coefficients @ np.array([1.0, node, node * node])) <= objective)
    problem = cp.Problem(cp.Minimize(objective), constraints)
    problem.solve(solver="CLARABEL")
    values = coefficients.value
    validation = max(np.max(np.abs(values @ np.array([1.0, node, node * node]))) for node in VALIDATION_NODES)
    return float(problem.value), float(validation), values


def fixed_basis_stress(nodes, directions, held_out, tolerance=0.01):
    rng = np.random.default_rng(20260809)
    bases = [np.eye(3)]
    for _ in range(63):
        basis, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        bases.append(basis)
    solutions = [solve_fixed_basis(nodes, directions, tolerance, basis) for basis in bases]
    index = int(np.argmin([solution[1] for solution in solutions]))
    objective, validation, coefficients = solutions[index]
    basis = bases[index]

    def residual(node, direction):
        multipliers = coefficients @ np.array([1.0, node, node * node])
        return np.linalg.norm((multipliers - 1) * (basis.T @ direction))

    return {
        "label": "randomized stress test, not a global certificate over orthogonal bases",
        "seed": 20260809,
        "basis_count": len(bases),
        "tolerance": tolerance,
        "best_index": index,
        "best_objective": objective,
        "best_validation_leakage": validation,
        "median_validation_leakage": float(np.median([solution[1] for solution in solutions])),
        "maximum_training_residual": float(max(residual(node, direction) for node, direction in zip(nodes, directions, strict=True))),
        "maximum_held_out_residual": float(max(residual(node, direction) for node, direction in zip(nodes, held_out, strict=True))),
        "best_basis": basis.tolist(),
    }


def main():
    manifest = json.loads(INPUT.read_text(encoding="utf-8"))
    calibration = manifest["calibration"]
    nodes = np.array(calibration["measured_nodes"])
    rational_nodes = np.array([float(value) for value in calibration["rational_nodes"]])
    training = np.array(calibration["training_directions"])
    held_out = np.array(calibration["held_out_directions"])
    held_out_cospectra = np.array(calibration["held_out_cospectral_matrices"])
    rational_training = np.array(calibration["rational_training_directions"], dtype=float)

    exact = np.array(
        [[[float(value) for value in row] for row in coefficient] for coefficient in build_coefficients()]
    )
    exact_dense_leakage = max(np.linalg.norm(matrix_value(exact, node), 2) for node in VALIDATION_NODES)
    exact_rational_residuals = pass_residuals(exact, rational_nodes, rational_training)
    measured_training_residuals = pass_residuals(exact, nodes, training)
    held_out_residuals = pass_residuals(exact, nodes, held_out)
    exact_commutator = max(
        np.linalg.norm(exact[left] @ exact[right] - exact[right] @ exact[left], 2)
        for left, right in combinations(range(3), 2)
    )
    normalized = rational_training / np.linalg.norm(rational_training, axis=1)[:, None]
    normalized_minors = [
        abs(np.linalg.det(normalized[list(subset)].T)) for subset in combinations(range(5), 3)
    ]
    angles = np.degrees(np.arccos(np.clip(np.abs(np.sum(training * held_out, axis=1)), 0, 1)))

    symmetric = solve_matrix(nodes, training, symmetric=True)
    general = solve_matrix(nodes, training, symmetric=False)
    for result in (symmetric, general):
        result["training_residuals"] = pass_residuals(result["coefficients"], nodes, training)
        result["held_out_residuals"] = pass_residuals(result["coefficients"], nodes, held_out)
        result["coefficients"] = result["coefficients"].tolist()

    result = {
        "schema_version": 1,
        "input": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "dataset_doi": manifest["dataset"]["doi"],
        "problem": {
            "degree": 2,
            "channels": 3,
            "training_calibrations": 5,
            "stop_node_interval": [-1.0, -0.929],
            "stop_frequency_interval_hz": [1099.1629585051116, 1250.0],
        },
        "data_stability": {
            "train_to_held_out_angles_degrees": angles.tolist(),
            "maximum_angle_degrees": float(angles.max()),
        },
        "exact_rational_witness": {
            "minimum_normalized_full_spark_minor": float(min(normalized_minors)),
            "continuum_certificate": "verification/verify_vbl_va001_witness.py proves leakage < 24/25 by exact rational arithmetic",
            "certified_leakage_upper_bound": 24 / 25,
            "dense_validation_leakage": float(exact_dense_leakage),
            "maximum_exact_rational_residual": float(max(exact_rational_residuals)),
            "maximum_measured_training_residual": float(max(measured_training_residuals)),
            "maximum_held_out_residual": float(max(held_out_residuals)),
            "held_out_residuals": held_out_residuals,
            "held_out_fdd_invariance": fdd_invariance_metrics(
                exact, nodes, training, held_out_cospectra
            ),
            "coefficient_frobenius_norms": [float(np.linalg.norm(value)) for value in exact],
            "maximum_commutator_norm": float(exact_commutator),
            "commuting_exact_calibration_lower_bound": 1.0,
        },
        "numerical_symmetric_optimum": symmetric,
        "numerical_general_nonsymmetric_optimum": general,
        "scalar_one_percent": solve_scalar(nodes, 0.01),
        "fixed_basis_one_percent_stress": fixed_basis_stress(nodes, training, held_out, 0.01),
        "interpretation": {
            "supported": "The rationalized measured calibration has an exact reciprocal symmetric noncommuting degree-two witness below 0.96, whereas every exactly calibrated commuting symmetric filter is the identity. The six-record held-out split is not used in fitting.",
            "not_supported": "The 64-basis search is not a global approximate-commuting optimum; no fault-classification, hardware, runtime, or deployment advantage is inferred.",
        },
        "versions": {"cvxpy": cp.__version__, "numpy": np.__version__},
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
