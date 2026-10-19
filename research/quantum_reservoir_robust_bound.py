#!/usr/bin/env python3
"""Multiprecision certificate for the approximate separation theorem.

The positive constant is deliberately computed from the finite matrices used
in the proof, rather than inferred from double-precision ranks.  The script
also evaluates the scalar central-branch construction that forces any such
constant to grow at least exponentially.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "quantum_reservoir" / "robust_separation.json"


def smallest_singular_value(matrix: mp.matrix) -> mp.mpf:
    values = mp.svd(matrix, compute_uv=False)
    return mp.mpf(values[values.rows - 1])


def pass_nodes(order: int, radius: mp.mpf, phase: mp.mpf) -> list[mp.mpf]:
    normalization_phase = mp.mpf(0) if order % 2 else mp.pi
    scale = (1 + radius) / (1 - radius)
    lower = 2 * mp.atan(scale * mp.tan(-mp.pi / 12))
    upper = -lower
    nodes: list[mp.mpf] = []
    for winding in range(-2 * order, 2 * order + 1):
        target = (2 * mp.pi * winding - normalization_phase - phase) / order
        if lower < target < upper:
            nodes.append(
                2 * mp.atan(((1 - radius) / (1 + radius)) * mp.tan(target / 2))
            )
    return sorted(nodes)[: order - 1]


def signatures_and_nodes(order: int, alpha: mp.mpf) -> tuple[list[mp.matrix], list[mp.mpc]]:
    radius = 1 - mp.exp(-alpha * order)
    gap = 1 - radius
    groups = [pass_nodes(order, radius, phase) for phase in (-gap, mp.mpf(0), gap)]
    omega = mp.exp(2j * mp.pi / 3)
    fourier = mp.matrix(
        [[omega ** (row * column) / mp.sqrt(3) for column in range(3)] for row in range(3)]
    )
    signatures: list[mp.matrix] = []
    nodes: list[mp.mpc] = []
    for channel, group in enumerate(groups):
        for theta in group:
            z = mp.exp(1j * theta)
            mixer = fourier * mp.diag([1, z, z * z]) * fourier.transpose_conj()
            signatures.append(mixer[:, channel])
            nodes.append(z)
    return signatures, nodes


def audit_order(order: int, alpha: mp.mpf) -> dict[str, str | int]:
    signatures, nodes = signatures_and_nodes(order, alpha)
    count = len(nodes)
    numerator_degree = order + 4

    gamma = mp.inf
    for triple in combinations(range(count), 3):
        matrix = mp.matrix(3, 3)
        for column, index in enumerate(triple):
            for row in range(3):
                matrix[row, column] = signatures[index][row]
        gamma = min(gamma, smallest_singular_value(matrix))

    # Adding rows can only increase the least singular value, so the minimum
    # over erasures of size at most two occurs for an erasure pair.
    chi = mp.inf
    for erased in combinations(range(count), 2):
        keep = [index for index in range(count) if index not in erased]
        vandermonde = mp.matrix(
            [[nodes[index] ** degree for degree in range(numerator_degree + 1)] for index in keep]
        )
        chi = min(chi, smallest_singular_value(vandermonde))

    sample_count = 8 * (order + 1)
    stop_nodes: list[mp.mpc] = []
    for lower, upper in ((-mp.pi, -5 * mp.pi / 6), (5 * mp.pi / 6, mp.pi)):
        for index in range(sample_count // 2):
            theta = lower + (index + mp.mpf("0.5")) * (upper - lower) / (sample_count // 2)
            stop_nodes.append(mp.exp(1j * theta))
    stop_vandermonde = mp.matrix(
        [[z ** degree for degree in range(order + 1)] for z in stop_nodes]
    )
    nu = smallest_singular_value(stop_vandermonde) / (
        mp.sqrt(sample_count) * mp.sqrt(order + 1)
    )

    robust_constant = mp.sqrt(3 * (count - 2) * (numerator_degree + 1)) / (
        gamma * chi * nu
    )
    pole_gap = mp.exp(-alpha * order)
    analytic_leakage = (mp.pi * order / 6 + mp.mpf("0.5")) * pole_gap
    scalar_calibration_error = mp.sin(pole_gap / 2)
    necessary_constant = max(
        mp.mpf(0), (1 - analytic_leakage) / scalar_calibration_error
    )

    def decimal(value: mp.mpf, digits: int = 18) -> str:
        return mp.nstr(value, digits, strip_zeros=False)

    return {
        "order": order,
        "calibration_count": count,
        "numerator_degree": numerator_degree,
        "full_spark_singular_margin_gamma": decimal(gamma),
        "pass_vandermonde_margin_chi": decimal(chi),
        "stop_denominator_margin_nu": decimal(nu),
        "proved_robust_constant_upper": decimal(robust_constant),
        "proved_robust_constant_log10_upper": decimal(mp.log10(robust_constant), 12),
        "central_scalar_calibration_error": decimal(scalar_calibration_error),
        "central_scalar_stop_bound": decimal(analytic_leakage),
        "necessary_robust_constant_lower": decimal(necessary_constant),
        "necessary_robust_constant_log10_lower": decimal(mp.log10(necessary_constant), 12),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", type=int, default=[5, 7, 9])
    parser.add_argument("--alpha", type=str, default="0.45")
    parser.add_argument("--precision", type=int, default=60)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    mp.mp.dps = args.precision
    alpha = mp.mpf(args.alpha)
    records = [audit_order(order, alpha) for order in args.orders]
    result = {
        "schema_version": 1,
        "precision_decimal_digits": args.precision,
        "alpha": args.alpha,
        "stop_denominator_sample_count": "8(S+1) midpoint samples",
        "records": records,
        "interpretation": {
            "positive_result": "For a comparator with calibration defect delta and sampled left-reduction defect beta, stop norm is at least max(0,1-C_S(delta+beta)).",
            "conditioning_warning": "The displayed upper constants are proof certificates, not claims of sharpness; the pass Vandermonde interpolation is severely ill-conditioned.",
            "impossibility_result": "The scalar central branch proves that every valid C_S is at least the recorded necessary lower bound, hence no S-uniform robustness is possible for these calibrations.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
