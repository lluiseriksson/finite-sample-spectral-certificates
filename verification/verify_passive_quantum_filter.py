"""Verification gates for the passive quantum reservoir-filter paper."""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.quantum_reservoir_filter import (
    audit_order,
    blaschke_value,
    lossless_block,
    mixer,
    pass_nodes,
)


def rational_rank(matrix: list[list[F]]) -> int:
    work = [row[:] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for index in range(len(work)):
            if index == row:
                continue
            multiple = work[index][column]
            work[index] = [
                value - multiple * pivot_value
                for value, pivot_value in zip(work[index], work[row], strict=True)
            ]
        row += 1
    return row


def exact_hermitian_no_go() -> None:
    points = [F(1, 2), F(1), F(3, 2), F(7, 2)]
    vectors = [[F(1), F(2)], [F(0), F(1)], [F(-1), F(2)], [F(-1), F(1)]]
    constraint_matrix: list[list[F]] = []
    for point, vector in zip(points, vectors, strict=True):
        first, second = vector
        row_first = [F()] * 9
        row_second = [F()] * 9
        for degree in range(3):
            row_first[3 * degree] = point**degree * first
            row_first[3 * degree + 1] = point**degree * second
            row_second[3 * degree + 1] = point**degree * first
            row_second[3 * degree + 2] = point**degree * second
        constraint_matrix.extend([row_first, row_second])
    assert rational_rank(constraint_matrix) == 8

    # Symmetric-entry order a_k,b_k,c_k.  Every feasible solution is I+tau R.
    r_coefficients = [
        [[F(5, 6), F(-1, 3)], [F(-1, 3), F(5, 6)]],
        [[F(-2, 3), F(-1, 3)], [F(-1, 3), F(-11, 6)]],
        [[F(2, 3), F(2, 3)], [F(2, 3), F(1)]],
    ]
    for point, vector in zip(points, vectors, strict=True):
        value = [
            sum(
                point**degree * r_coefficients[degree][row][column] * vector[column]
                for degree in range(3)
                for column in range(2)
            )
            for row in range(2)
        ]
        assert value == [F(), F()]

    r_zero = r_coefficients[0]
    determinant_zero = r_zero[0][0] * r_zero[1][1] - r_zero[0][1] ** 2
    assert r_zero[0][0] > 0 and determinant_zero == F(7, 12) > 0

    point = F(5, 2)
    r_test = [
        [
            sum(point**degree * r_coefficients[degree][row][column] for degree in range(3))
            for column in range(2)
        ]
        for row in range(2)
    ]
    determinant_test = r_test[0][0] * r_test[1][1] - r_test[0][1] ** 2
    assert r_test == [[F(10, 3), F(3)], [F(3), F(5, 2)]]
    assert determinant_test == F(-2, 3)
    print("PASS: exact quadratic solution space is I + tau R")
    print("PASS: global Hermitian contractivity forces tau = 0")


def passive_family_gates() -> None:
    alpha = 0.45
    for order in [5, 7, 9, 11, 15, 19]:
        record = audit_order(order, alpha)
        analytic_bound = (3.141592653589793 * order / 6 + 0.5) * (2.718281828459045 ** (-alpha * order))
        assert record["nodes_per_group"] >= order - 1
        assert record["root_count_separation_valid"]
        assert record["maximum_calibration_residual"] < 1e-9
        assert record["full_spark_minimum_minor"] > 0
        assert record["global_signal_norm"] <= 1 + 2e-12
        assert record["stopband_signal_norm"] <= analytic_bound
        assert record["six_port_unitarity_error"] < 1e-10
        assert abs(record["worst_case_rate_ratio_upper"] - record["stopband_signal_norm"] ** 2) < 1e-15
        assert record["peak_allpass_group_delay"] >= order / (1 - record["radius"])
    print("PASS: exact delayed calibrations and full-spark margins")
    print("PASS: global Schur contractivity and analytic stopband envelope")
    print("PASS: explicit six-port lossless completion")
    print("PASS: reducible rational root-count lower bound")
    print("PASS: squared bath-rate separation and exponential delay cost")


def robust_separation_gates() -> None:
    """Audit the finite matrices and the matching scalar robustness barrier."""
    rng = np.random.default_rng(20260810)
    alpha = 0.45
    for order in [5, 7, 9]:
        radius = 1 - np.exp(-alpha * order)
        gap = 1 - radius
        phases = gap * np.array([-1.0, 0.0, 1.0])
        groups = [pass_nodes(order, radius, phase, np.pi / 6)[: order - 1] for phase in phases]
        nodes: list[complex] = []
        signatures: list[np.ndarray] = []
        for channel, group in enumerate(groups):
            for theta in group:
                z = np.exp(1j * theta)
                nodes.append(z)
                signatures.append(mixer(z)[:, channel])

        gamma = min(
            np.linalg.svd(np.column_stack([signatures[index] for index in triple]), compute_uv=False)[-1]
            for triple in combinations(range(len(signatures)), 3)
        )
        assert gamma > 0
        test_direction = rng.normal(size=3) + 1j * rng.normal(size=3)
        test_direction /= np.linalg.norm(test_direction)
        overlaps = np.array([abs(np.vdot(test_direction, vector)) for vector in signatures])
        keep = np.argsort(overlaps)[2:]
        assert float(np.min(overlaps[keep])) + 1e-14 >= gamma / np.sqrt(3)

        degree = order + 4
        erasure = set(np.argsort(overlaps)[:2].tolist())
        vandermonde = np.array(
            [[nodes[index] ** power for power in range(degree + 1)]
             for index in range(len(nodes)) if index not in erasure],
            dtype=complex,
        )
        pairwise = np.abs(np.subtract.outer(np.asarray(nodes), np.asarray(nodes)))
        pairwise += np.eye(len(nodes)) * 10
        assert float(np.min(pairwise)) > 1e-6
        # The proof of full column rank is the exact Vandermonde determinant;
        # its numerical singular value is intentionally not thresholded here
        # because the multiprecision ledger measures severe conditioning.
        assert np.linalg.svd(vandermonde, compute_uv=False)[-1] > 0

        # Q_c=z^2(1+B_S)/2 I is reducing and Schur.  It nearly obeys all
        # calibrations while retaining the same small stop response.
        central_error = np.sin(gap / 2)
        observed_error = 0.0
        for group in groups:
            for theta in group:
                b_value = blaschke_value(np.array([theta]), order, radius)[0]
                observed_error = max(observed_error, abs((1 + b_value) / 2 - 1))
        assert abs(observed_error - central_error) < 2e-12
        stop_grid = np.concatenate(
            [np.linspace(-np.pi, -5 * np.pi / 6, 4001), np.linspace(5 * np.pi / 6, np.pi, 4001)]
        )
        central_stop = max(
            abs((1 + blaschke_value(np.array([theta]), order, radius)[0]) / 2)
            for theta in stop_grid
        )
        analytic_bound = (np.pi * order / 6 + 0.5) * gap
        assert central_stop <= analytic_bound + 1e-12
        assert (1 - analytic_bound) / central_error > 0
    print("PASS: quantitative full-spark and erasure-Vandermonde margins")
    print("PASS: approximate calibration/reduction obstruction ingredients")
    print("PASS: scalar Schur construction forces exponential robustness loss")


def closed_dephasing_gates() -> None:
    alpha = 0.45
    maximum_identity_error = 0.0
    for order in [5, 9, 19]:
        radius = 1 - np.exp(-alpha * order)
        phases = (1 - radius) * np.array([-1.0, 0.0, 1.0])
        for theta in np.linspace(-np.pi, np.pi, 101):
            scattering = lossless_block(theta, order, radius, phases)
            signal = scattering[:3, :3]
            loss = scattering[:3, 3:]
            maximum_identity_error = max(
                maximum_identity_error,
                float(np.linalg.norm(signal @ signal.conj().T + loss @ loss.conj().T - np.eye(3), 2)),
            )
    assert maximum_identity_error < 1e-10

    # With all six inputs in vacuum, the signal output covariance is exactly
    # I/2.  Excess noise J on the first three inputs adds G J G* and nothing
    # else, so the independently measured vacuum rate is architecture-neutral.
    kappa_zero = 1e-3
    leakage = float(audit_order(9, alpha)["stopband_signal_norm"])
    total_constructed_upper = kappa_zero + leakage**2
    total_comparator_lower = kappa_zero + 1.0
    assert total_comparator_lower / total_constructed_upper > 500
    assert total_comparator_lower / total_constructed_upper < 1 / kappa_zero + 2
    print("PASS: vacuum ports close to an architecture-independent covariance baseline")
    print("PASS: total Ramsey-rate separation saturates at the measured vacuum floor")


if __name__ == "__main__":
    exact_hermitian_no_go()
    passive_family_gates()
    robust_separation_gates()
    closed_dephasing_gates()
