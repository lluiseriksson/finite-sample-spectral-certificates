"""Independently replay the exact-rational certificate stored as JSON."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/")
    return Fraction(int(numerator), int(denominator))


def gram_from_hex(rows: list[list[str]]) -> list[list[Fraction]]:
    z = [[Fraction.from_float(float.fromhex(value)) for value in row] for row in rows]
    dimension = len(z[0])
    return [
        [sum((row[i] * row[j] for row in z), Fraction()) / len(z)
         for j in range(dimension)]
        for i in range(dimension)
    ]


def ldl_positive(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    lower = [[Fraction() for _ in range(size)] for _ in range(size)]
    pivots: list[Fraction] = []
    for i in range(size):
        lower[i][i] = Fraction(1)
        pivot = matrix[i][i] - sum(
            (lower[i][k] ** 2 * pivots[k] for k in range(i)), Fraction()
        )
        assert pivot > 0, (i, pivot)
        pivots.append(pivot)
        for j in range(i + 1, size):
            lower[j][i] = (
                matrix[j][i]
                - sum((lower[j][k] * lower[i][k] * pivots[k] for k in range(i)), Fraction())
            ) / pivot
    return min(pivots)


def multiply_trace(left: list[list[Fraction]], right: list[list[Fraction]]) -> Fraction:
    return sum(
        (left[i][j] * right[j][i] for i in range(len(left)) for j in range(len(left))),
        Fraction(),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", nargs="?", type=Path, default=Path("results/verification/rational_witness.json"))
    args = parser.parse_args()
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    parameters = artifact["parameters"]
    degree, channels = parameters["degree"], parameters["channels"]
    theta = parse_fraction(parameters["theta"])
    band_lower = parse_fraction(parameters["band_lower"])
    band_upper = parse_fraction(parameters["band_upper"])

    gram = gram_from_hex(artifact["data"]["z_float_hex"])
    scales = [Fraction.from_float(float.fromhex(value)) for value in artifact["data"]["scales_float_hex"]]
    empirical = [
        [gram[i][j] / (scales[i] * scales[j]) for j in range(len(scales))]
        for i in range(len(scales))
    ]
    interval_lower = [[band_lower * value for value in row] for row in empirical]
    interval_upper = [[band_upper * value for value in row] for row in empirical]

    rational = artifact["rational_multipliers"]
    denominator = rational["common_denominator"]
    matrices = []
    for name in ("lower_numerators", "upper_numerators", "localizer_numerators"):
        matrices.append([[Fraction(value, denominator) for value in row] for row in rational[name]])
    y_lower, y_upper, y_local = matrices

    minimum_pivots = [
        ldl_positive(y_lower),
        ldl_positive(y_upper),
        ldl_positive(y_local),
        ldl_positive(interval_lower),
        ldl_positive(interval_upper),
    ]
    dimension = len(scales)
    local_adjoint = [[Fraction() for _ in range(dimension)] for _ in range(dimension)]
    for k in range(2 * degree + 2):
        for a in range(channels):
            for b in range(channels):
                block_value = Fraction()
                for i in range(degree + 1):
                    for j in range(degree + 1):
                        coefficient = theta * int(i + j == k) - int(i + j + 1 == k)
                        block_value += coefficient * y_local[i * channels + a][j * channels + b]
                row, column = k * channels + a, k * channels + b
                local_adjoint[row][column] = scales[row] * block_value * scales[column]
    residual = [
        [y_lower[i][j] - y_upper[i][j] + local_adjoint[i][j] for j in range(dimension)]
        for i in range(dimension)
    ]
    beta = -multiply_trace(y_lower, interval_lower) + multiply_trace(y_upper, interval_upper)
    maximum_residual = max(abs(value) for row in residual for value in row)
    trace_upper = sum(interval_upper[i][i] for i in range(dimension))
    residual_bound = dimension * maximum_residual * trace_upper
    certificate_upper = beta + residual_bound
    assert certificate_upper < 0
    assert parse_fraction(artifact["exact_checks"]["dual_margin_beta"]) == beta
    assert parse_fraction(artifact["exact_checks"]["residual_bound"]) == residual_bound
    assert parse_fraction(artifact["exact_checks"]["certificate_upper"]) == certificate_upper
    print(f"PASS: all five matrices have positive exact LDL pivots (minimum {float(min(minimum_pivots)):.3e})")
    print(f"PASS: exact conservative certificate upper = {float(certificate_upper):.12e} < 0")


if __name__ == "__main__":
    main()
