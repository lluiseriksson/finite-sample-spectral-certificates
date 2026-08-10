"""Verification gates for the passive quantum reservoir-filter paper."""

from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.quantum_reservoir_filter import audit_order


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


if __name__ == "__main__":
    exact_hermitian_no_go()
    passive_family_gates()
