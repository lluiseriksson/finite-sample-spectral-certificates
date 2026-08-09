"""Export one replayable rational SDP contradiction certificate.

The statistical theorem uses eta=(sqrt(d)+sqrt(2 log(2/alpha)))/sqrt(n).
For d=12, n=500, alpha=0.05 we replace it by the rigorously larger rational
eta_bar=3/10.  The resulting Loewner interval is wider, so infeasibility of
that interval is also a valid rejection for the analytic confidence band.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from pilot_annni_block_hotelling import exact_half_moments, joint_sketch_covariance
from pilot_wishart_loewner import dual_certificate


RATIONAL_DENOMINATOR = 10**14
DIAGONAL_SHIFT_UNITS = 1000


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def rational_matrix(values: np.ndarray) -> list[list[Fraction]]:
    symmetric = (values + values.T) / 2.0
    answer = [
        [Fraction(round(float(symmetric[i, j]) * RATIONAL_DENOMINATOR), RATIONAL_DENOMINATOR)
         for j in range(symmetric.shape[1])]
        for i in range(symmetric.shape[0])
    ]
    shift = Fraction(DIAGONAL_SHIFT_UNITS, RATIONAL_DENOMINATOR)
    for i in range(len(answer)):
        answer[i][i] += shift
    return answer


def exact_gram(z: np.ndarray) -> list[list[Fraction]]:
    rows = [[Fraction.from_float(float(value)) for value in row] for row in z]
    dimension = z.shape[1]
    return [
        [sum((row[i] * row[j] for row in rows), Fraction()) / z.shape[0]
         for j in range(dimension)]
        for i in range(dimension)
    ]


def scaled_gram(
    gram: list[list[Fraction]], scales: list[Fraction]
) -> list[list[Fraction]]:
    return [
        [gram[i][j] / (scales[i] * scales[j]) for j in range(len(gram))]
        for i in range(len(gram))
    ]


def scale_matrix(matrix: list[list[Fraction]], scalar: Fraction) -> list[list[Fraction]]:
    return [[scalar * value for value in row] for row in matrix]


def trace_product(left: list[list[Fraction]], right: list[list[Fraction]]) -> Fraction:
    return sum(
        (left[i][j] * right[j][i] for i in range(len(left)) for j in range(len(left))),
        Fraction(),
    )


def ldl_pivots(matrix: list[list[Fraction]]) -> list[Fraction]:
    size = len(matrix)
    lower = [[Fraction() for _ in range(size)] for _ in range(size)]
    pivots = [Fraction() for _ in range(size)]
    for i in range(size):
        lower[i][i] = Fraction(1)
        pivots[i] = matrix[i][i] - sum(
            (lower[i][k] * lower[i][k] * pivots[k] for k in range(i)), Fraction()
        )
        if pivots[i] <= 0:
            raise RuntimeError(f"matrix is not positive definite at exact LDL pivot {i}")
        for j in range(i + 1, size):
            numerator = matrix[j][i] - sum(
                (lower[j][k] * lower[i][k] * pivots[k] for k in range(i)), Fraction()
            )
            lower[j][i] = numerator / pivots[i]
    return pivots


def adjoint(
    local: list[list[Fraction]],
    scales: list[Fraction],
    degree: int,
    channels: int,
    theta: Fraction,
) -> list[list[Fraction]]:
    dimension = len(scales)
    answer = [[Fraction() for _ in range(dimension)] for _ in range(dimension)]
    for k in range(2 * degree + 2):
        for a in range(channels):
            for b in range(channels):
                value = Fraction()
                for i in range(degree + 1):
                    for j in range(degree + 1):
                        coefficient = theta * int(i + j == k) - int(i + j + 1 == k)
                        if coefficient:
                            value += coefficient * local[i * channels + a][j * channels + b]
                row, column = k * channels + a, k * channels + b
                answer[row][column] = scales[row] * value * scales[column]
    return answer


def add_residual(
    lower: list[list[Fraction]],
    upper: list[list[Fraction]],
    local_adjoint: list[list[Fraction]],
) -> list[list[Fraction]]:
    return [
        [lower[i][j] - upper[i][j] + local_adjoint[i][j] for j in range(len(lower))]
        for i in range(len(lower))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/verification/rational_witness.json"))
    parser.add_argument("--seed", type=int, default=20260831)
    args = parser.parse_args()

    length, degree, channels, sample_count = 8, 2, 2, 500
    alpha, theta = 0.05, Fraction(3, 5)
    eta_bar = Fraction(3, 10)
    band_lower = Fraction(1, 1) / (1 + eta_bar) ** 2
    band_upper = Fraction(1, 1) / (1 - eta_bar) ** 2

    half, metadata = exact_half_moments(length, degree, 0.2, 1.0, 0.37, 2.2)
    population = joint_sketch_covariance(half, 2 * degree + 1)
    rng = np.random.default_rng(args.seed)
    z = rng.multivariate_normal(np.zeros(population.shape[0]), population, size=sample_count)
    sample_gram_float = z.T @ z / sample_count
    solved = dual_certificate(
        sample_gram_float,
        degree,
        channels,
        float(theta),
        alpha,
        sample_count,
        return_matrices=True,
        band_eta=float(eta_bar),
    )
    payload = solved.pop("matrix_payload")

    gram = exact_gram(z)
    scales = [Fraction.from_float(float(value)) for value in payload["scales"]]
    empirical = scaled_gram(gram, scales)
    interval_lower = scale_matrix(empirical, band_lower)
    interval_upper = scale_matrix(empirical, band_upper)
    y_lower = rational_matrix(np.asarray(payload["lower_multiplier"], dtype=float))
    y_upper = rational_matrix(np.asarray(payload["upper_multiplier"], dtype=float))
    y_local = rational_matrix(np.asarray(payload["localizer_multiplier"], dtype=float))

    pivots_lower = ldl_pivots(y_lower)
    pivots_upper = ldl_pivots(y_upper)
    pivots_local = ldl_pivots(y_local)
    pivots_interval_lower = ldl_pivots(interval_lower)
    pivots_interval_upper = ldl_pivots(interval_upper)

    residual = add_residual(y_lower, y_upper, adjoint(y_local, scales, degree, channels, theta))
    beta = -trace_product(y_lower, interval_lower) + trace_product(y_upper, interval_upper)
    max_residual = max(abs(value) for row in residual for value in row)
    trace_upper = sum(interval_upper[i][i] for i in range(len(interval_upper)))
    # For 0 <= Q <= B, sum_ij |Q_ij| <= d tr(Q) <= d tr(B).
    residual_bound = len(interval_upper) * max_residual * trace_upper
    certificate_upper = beta + residual_bound
    if certificate_upper >= 0:
        raise RuntimeError(f"rational certificate failed: {float(certificate_upper):.3e}")

    # Exact elementary proof that eta(d=12,n=500,alpha=.05) < 3/10.
    # exp(15/4)>40 follows from a finite positive Taylor partial sum.
    exp_lower = sum((Fraction(15, 4) ** k / math.factorial(k) for k in range(13)), Fraction())
    if not exp_lower > 40:
        raise RuntimeError("Taylor bound did not prove log(40)<15/4")
    if not Fraction(49, 4) > 12 or not Fraction(121, 16) > Fraction(15, 2):
        raise RuntimeError("radical upper bounds failed")
    if not Fraction(45) > Fraction(625, 16):
        raise RuntimeError("eta comparison failed")

    output = {
        "schema_version": 1,
        "claim": "Exact-rational infeasibility witness for the eta_bar=3/10 Loewner interval.",
        "parameters": {
            "length": length,
            "degree": degree,
            "channels": channels,
            "sample_count": sample_count,
            "alpha": alpha,
            "theta": fraction_text(theta),
            "eta_upper": fraction_text(eta_bar),
            "band_lower": fraction_text(band_lower),
            "band_upper": fraction_text(band_upper),
            "seed": args.seed,
        },
        "eta_bound_proof": {
            "exp_15_over_4_taylor_0_to_12": fraction_text(exp_lower),
            "logic": [
                "exp(15/4)>40, hence log(40)<15/4",
                "sqrt(12)<7/2 and sqrt(2 log 40)<11/4",
                "their sum is <25/4<3 sqrt(5)=(3/10)sqrt(500)",
            ],
        },
        "data": {
            "z_float_hex": [[float(value).hex() for value in row] for row in z],
            "scales_float_hex": [float(value).hex() for value in payload["scales"]],
        },
        "rational_multipliers": {
            "common_denominator": RATIONAL_DENOMINATOR,
            "diagonal_shift_units": DIAGONAL_SHIFT_UNITS,
            "lower_numerators": [[int(value * RATIONAL_DENOMINATOR) for value in row] for row in y_lower],
            "upper_numerators": [[int(value * RATIONAL_DENOMINATOR) for value in row] for row in y_upper],
            "localizer_numerators": [[int(value * RATIONAL_DENOMINATOR) for value in row] for row in y_local],
        },
        "exact_checks": {
            "minimum_ldl_pivot_lower": fraction_text(min(pivots_lower)),
            "minimum_ldl_pivot_upper": fraction_text(min(pivots_upper)),
            "minimum_ldl_pivot_localizer": fraction_text(min(pivots_local)),
            "minimum_ldl_pivot_interval_lower": fraction_text(min(pivots_interval_lower)),
            "minimum_ldl_pivot_interval_upper": fraction_text(min(pivots_interval_upper)),
            "dual_margin_beta": fraction_text(beta),
            "maximum_absolute_stationarity_residual": fraction_text(max_residual),
            "residual_bound": fraction_text(residual_bound),
            "certificate_upper": fraction_text(certificate_upper),
            "certificate_upper_float": float(certificate_upper),
        },
        "floating_solver_diagnostics": solved,
        "model": metadata,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(f"exact certificate upper: {float(certificate_upper):.12e}")


if __name__ == "__main__":
    main()
