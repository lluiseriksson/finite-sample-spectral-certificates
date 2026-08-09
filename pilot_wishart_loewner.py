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


def dual_certificate(
    sample_gram: np.ndarray,
    degree: int,
    channels: int,
    theta: float,
    alpha: float,
    sample_count: int,
) -> dict[str, object]:
    """Solve the normalized semidefinite alternative and report its residual.

    For A<=Q<=B and L(DQD)>=0, PSD multipliers Y1,Y2,Y3 obey
    Y1-Y2+L_D^*(Y3)=0.  A negative value of
    -<Y1,A>+<Y2,B> is an explicit incompatibility certificate.
    """
    _, lower, upper = wishart_band(sample_gram, alpha, sample_count)
    dimension = sample_gram.shape[0]
    scales = np.sqrt(np.maximum(np.diag(sample_gram), np.finfo(float).tiny))
    inverse = np.diag(1.0 / scales)
    empirical = inverse @ sample_gram @ inverse
    empirical = (empirical + empirical.T) / 2
    localizer_size = (degree + 1) * channels

    y_lower = cp.Variable((dimension, dimension), PSD=True)
    y_upper = cp.Variable((dimension, dimension), PSD=True)
    y_local = cp.Variable((localizer_size, localizer_size), PSD=True)

    diagonal_blocks: list[cp.Expression] = []
    for k in range(2 * degree + 2):
        block: cp.Expression = cp.Constant(np.zeros((channels, channels)))
        for i in range(degree + 1):
            for j in range(degree + 1):
                coefficient = theta * (i + j == k) - (i + j + 1 == k)
                if coefficient:
                    block = block + coefficient * y_local[
                        i * channels : (i + 1) * channels,
                        j * channels : (j + 1) * channels,
                    ]
        block_scales = np.diag(scales[k * channels : (k + 1) * channels])
        diagonal_blocks.append(block_scales @ block @ block_scales)

    zero = np.zeros((channels, channels))
    adjoint = cp.bmat(
        [
            [diagonal_blocks[i] if i == j else zero for j in range(2 * degree + 2)]
            for i in range(2 * degree + 2)
        ]
    )
    stationarity = y_lower - y_upper + adjoint
    objective = -cp.trace(y_lower @ (lower * empirical)) + cp.trace(
        y_upper @ (upper * empirical)
    )
    normalization = cp.trace(y_lower) + cp.trace(y_upper) + cp.trace(y_local)
    problem = cp.Problem(cp.Minimize(objective), [stationarity == 0, normalization == 1])
    problem.solve(solver="CLARABEL", tol_gap_abs=1e-9, tol_feas=1e-9, max_iter=1000)
    if problem.status not in {"optimal", "optimal_inaccurate"}:
        raise RuntimeError(f"dual certificate: {problem.status}")

    yl = (y_lower.value + y_lower.value.T) / 2
    yu = (y_upper.value + y_upper.value.T) / 2
    yz = (y_local.value + y_local.value.T) / 2

    def numerical_adjoint(multiplier: np.ndarray) -> np.ndarray:
        answer = np.zeros((dimension, dimension))
        for k in range(2 * degree + 2):
            block_value = np.zeros((channels, channels))
            for i in range(degree + 1):
                for j in range(degree + 1):
                    coefficient = theta * (i + j == k) - (i + j + 1 == k)
                    if coefficient:
                        block_value += coefficient * multiplier[
                            i * channels : (i + 1) * channels,
                            j * channels : (j + 1) * channels,
                        ]
            indices = slice(k * channels, (k + 1) * channels)
            block_scales = np.diag(scales[indices])
            answer[indices, indices] = block_scales @ block_value @ block_scales
        return answer

    adjoint_value = numerical_adjoint(yz)
    residual = yl - yu + adjoint_value
    value = float(-np.trace(yl @ (lower * empirical)) + np.trace(yu @ (upper * empirical)))

    # Construct a feasible dual witness rather than merely measuring residuals.
    # First project each multiplier onto the PSD cone.  If the resulting
    # stationarity residual is R=R_+-R_-, adding R_- to Y_lower and R_+ to
    # Y_upper cancels R exactly while preserving positive semidefiniteness.
    def positive_part(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2)
        positive = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T
        negative = (eigenvectors * np.maximum(-eigenvalues, 0.0)) @ eigenvectors.T
        return positive, negative

    yl_psd, _ = positive_part(yl)
    yu_psd, _ = positive_part(yu)
    yz_psd, _ = positive_part(yz)
    repair_residual = yl_psd - yu_psd + numerical_adjoint(yz_psd)
    residual_positive, residual_negative = positive_part(repair_residual)
    yl_repaired = yl_psd + residual_negative
    yu_repaired = yu_psd + residual_positive
    repaired_scale = float(np.trace(yl_repaired) + np.trace(yu_repaired) + np.trace(yz_psd))
    yl_repaired /= repaired_scale
    yu_repaired /= repaired_scale
    yz_repaired = yz_psd / repaired_scale
    repaired_residual = yl_repaired - yu_repaired + numerical_adjoint(yz_repaired)
    repaired_value = float(
        -np.trace(yl_repaired @ (lower * empirical))
        + np.trace(yu_repaired @ (upper * empirical))
    )
    return {
        "dual_status": problem.status,
        "dual_value": value,
        "dual_stationarity_fro": float(np.linalg.norm(residual)),
        "dual_normalization_error": float(abs(np.trace(yl) + np.trace(yu) + np.trace(yz) - 1.0)),
        "dual_min_psd_eigenvalue": float(
            min(np.linalg.eigvalsh(yl)[0], np.linalg.eigvalsh(yu)[0], np.linalg.eigvalsh(yz)[0])
        ),
        "repaired_dual_value": repaired_value,
        "repaired_stationarity_fro": float(np.linalg.norm(repaired_residual)),
        "repaired_normalization_error": float(
            abs(np.trace(yl_repaired) + np.trace(yu_repaired) + np.trace(yz_repaired) - 1.0)
        ),
        "repaired_min_psd_eigenvalue": float(
            min(
                np.linalg.eigvalsh(yl_repaired)[0],
                np.linalg.eigvalsh(yu_repaired)[0],
                np.linalg.eigvalsh(yz_repaired)[0],
            )
        ),
    }


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
            certificates = []
            for _ in range(args.ensembles):
                z = rng.multivariate_normal(np.zeros(joint.shape[0]), joint, size=sample_count)
                sample_gram = z.T @ z / sample_count
                certificate = dual_certificate(
                    sample_gram, args.degree, channels, theta, args.alpha, sample_count
                )
                detected = certificate["repaired_dual_value"] < -1e-9
                detections += detected
                status = str(certificate["dual_status"])
                statuses[status] = statuses.get(status, 0) + 1
                certificates.append(certificate)
            rows.append(
                {
                    "sample_count": sample_count,
                    "gap_factor": factor,
                    "eta": eta,
                    "loewner_lower": lower,
                    "loewner_upper": upper,
                    "exact_localizer_min": exact_localizer_min(exact_moments, args.degree, channels, theta),
                    "detection_rate": detections / args.ensembles,
                    "dual_detection_threshold": -1e-9,
                    "dual_statuses": statuses,
                    "dual_value_range": [
                        min(item["dual_value"] for item in certificates),
                        max(item["dual_value"] for item in certificates),
                    ],
                    "repaired_dual_value_range": [
                        min(item["repaired_dual_value"] for item in certificates),
                        max(item["repaired_dual_value"] for item in certificates),
                    ],
                    "maximum_stationarity_fro": max(
                        item["dual_stationarity_fro"] for item in certificates
                    ),
                    "maximum_normalization_error": max(
                        item["dual_normalization_error"] for item in certificates
                    ),
                    "minimum_reported_psd_eigenvalue": min(
                        item["dual_min_psd_eigenvalue"] for item in certificates
                    ),
                    "maximum_repaired_stationarity_fro": max(
                        item["repaired_stationarity_fro"] for item in certificates
                    ),
                    "maximum_repaired_normalization_error": max(
                        item["repaired_normalization_error"] for item in certificates
                    ),
                    "minimum_repaired_psd_eigenvalue": min(
                        item["repaired_min_psd_eigenvalue"] for item in certificates
                    ),
                    "first_certificate": certificates[0],
                }
            )
    return {
        "schema_version": 3,
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
