"""Production-grid driver for distribution-free ANNNI block certificates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Chebyshev, Polynomial

from pilot_annni_block_hotelling import (
    block_design,
    exact_half_moments,
    flatten_moments,
    joint_sketch_covariance,
    parameter_index,
    robust_tests,
    sketch_feature_covariance,
)


def feature_samples(
    rng: np.random.Generator, joint: np.ndarray, count: int, maximum: int, channels: int
) -> np.ndarray:
    z = rng.multivariate_normal(np.zeros(joint.shape[0]), joint, size=count)
    z = z.reshape(count, maximum + 1, channels)
    params = parameter_index(maximum, channels)
    return np.column_stack([z[:, k, a] * z[:, k, b] for k, a, b in params])


def unflatten(values: np.ndarray, maximum: int, channels: int) -> np.ndarray:
    answer = np.zeros((maximum + 1, channels, channels))
    for value, (k, a, b) in zip(values, parameter_index(maximum, channels), strict=True):
        answer[k, a, b] = answer[k, b, a] = value
    return answer


def ritz_edge(moments: np.ndarray, degree: int) -> float:
    channels = moments.shape[1]
    size = (degree + 1) * channels
    h = np.empty((size, size))
    g = np.empty_like(h)
    for i in range(degree + 1):
        for j in range(degree + 1):
            h[i * channels : (i + 1) * channels, j * channels : (j + 1) * channels] = moments[i + j]
            g[i * channels : (i + 1) * channels, j * channels : (j + 1) * channels] = moments[i + j + 1]
    h = (h + h.T) / 2
    values, vectors = np.linalg.eigh(h)
    keep = values > max(values.max(), 1.0) * 1e-9
    if not np.any(keep):
        return float("nan")
    whitening = vectors[:, keep] / np.sqrt(values[keep])
    reduced = (whitening.T @ g @ whitening + whitening.T @ g.T @ whitening) / 2
    return float(np.linalg.eigvalsh(reduced)[-1])


def fixed_chebyshev_upper(
    observed: np.ndarray,
    covariance_mean: np.ndarray,
    design: np.ndarray,
    theta: float,
    target: float,
    degree: int,
    channels: int,
    q: float,
) -> float:
    polynomial = Chebyshev.basis(degree, domain=[0.0, theta]).convert(kind=Polynomial)
    coefficients = np.pad(
        polynomial.coef / polynomial(target), (0, degree + 1 - len(polynomial.coef))
    )
    vector = np.zeros((degree + 1) * channels)
    vector[1::channels] = coefficients  # parity-odd Z channel
    x = np.outer(vector, vector)
    x /= np.trace(x)
    matrix = np.einsum("k,kij->ij", observed, design)
    sensitivity = np.einsum("kij,ij->k", design, x)
    root = np.linalg.cholesky(covariance_mean + 1e-16 * np.eye(len(observed)))
    return float(np.trace(matrix @ x) + q * np.linalg.norm(root.T @ sensitivity))


def run_grid(args: argparse.Namespace) -> dict[str, object]:
    channels = 2
    rows = []
    for length in args.lengths:
        half_all, metadata = exact_half_moments(
            length, max(args.degrees), args.tau, args.j1, args.j2, args.hx
        )
        true_edge = float(np.exp(-args.tau * metadata["true_gap"]))
        for degree in args.degrees:
            maximum = 2 * degree + 1
            half = half_all[: 2 * maximum + 1]
            exact_matrices = half[::2][: maximum + 1]
            exact_flat = flatten_moments(exact_matrices)
            joint = joint_sketch_covariance(half, maximum)
            feature_covariance = sketch_feature_covariance(joint, maximum, channels)
            dimension = len(exact_flat)
            q = float(np.sqrt(dimension / (1.0 - args.coverage)))
            for sample_count in args.sample_counts:
                rng = np.random.default_rng(
                    np.random.SeedSequence([args.seed, length, degree, sample_count])
                )
                designs = {
                    factor: block_design(
                        float(np.exp(-args.tau * factor * metadata["true_gap"])), degree, channels
                    )
                    for factor in args.gap_factors
                }
                counters = {
                    factor: {"psd": 0, "rank_one": 0, "fixed_chebyshev": 0, "ritz": 0}
                    for factor in args.gap_factors
                }
                first = []
                for ensemble in range(args.ensembles):
                    observed = feature_samples(rng, joint, sample_count, maximum, channels).mean(axis=0)
                    matrices = unflatten(observed, maximum, channels)
                    edge_hat = ritz_edge(matrices, degree)
                    for factor, design in designs.items():
                        theta = float(np.exp(-args.tau * factor * metadata["true_gap"]))
                        robust = robust_tests(
                            observed, feature_covariance / sample_count, design, q
                        )
                        fixed = fixed_chebyshev_upper(
                            observed,
                            feature_covariance / sample_count,
                            design,
                            theta,
                            true_edge,
                            degree,
                            channels,
                            q,
                        )
                        counters[factor]["psd"] += robust["optimized_psd_upper"] < 0
                        counters[factor]["rank_one"] += robust["adaptive_rank_one_upper"] < 0
                        counters[factor]["fixed_chebyshev"] += fixed < 0
                        counters[factor]["ritz"] += edge_hat > theta
                        if ensemble == 0:
                            exact_localizer = np.einsum("k,kij->ij", exact_flat, design)
                            first.append(
                                {
                                    "gap_factor": factor,
                                    "theta": theta,
                                    "exact_localizer_min": float(np.linalg.eigvalsh(exact_localizer)[0]),
                                    "ritz_edge": edge_hat,
                                    "fixed_chebyshev_upper": fixed,
                                    **robust,
                                }
                            )
                for factor in args.gap_factors:
                    rows.append(
                        {
                            "length": length,
                            "degree": degree,
                            "sample_count": sample_count,
                            "gap_factor": factor,
                            "true_gap": metadata["true_gap"],
                            "true_edge": true_edge,
                            "moment_dimension": dimension,
                            "ensembles": args.ensembles,
                            "detection_rates": {
                                key: value / args.ensembles for key, value in counters[factor].items()
                            },
                            "first_ensemble": next(item for item in first if item["gap_factor"] == factor),
                        }
                    )
    return {
        "schema_version": 1,
        "setup": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in vars(args).items()
        },
        "method_boundary": {
            "psd_rank_one_chebyshev": "distribution-free covariance-known certificates",
            "ritz": "point estimator; not a confidence certificate",
            "fixed_chebyshev_target": "oracle true edge, deliberately favorable baseline",
        },
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", nargs="+", type=int, default=[6, 8, 10, 12, 14, 16])
    parser.add_argument("--degrees", nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--sample-counts", nargs="+", type=int, default=[250, 500, 1000, 2000, 5000])
    parser.add_argument("--gap-factors", nargs="+", type=float, default=[1.0, 1.05, 1.15, 1.35])
    parser.add_argument("--ensembles", type=int, default=500)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--j1", type=float, default=1.0)
    parser.add_argument("--j2", type=float, default=0.37)
    parser.add_argument("--hx", type=float, default=2.2)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("production_annni_campaign.json"))
    args = parser.parse_args()
    payload = run_grid(args)
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"wrote {len(payload['rows'])} rows to {args.output}")


if __name__ == "__main__":
    main()
