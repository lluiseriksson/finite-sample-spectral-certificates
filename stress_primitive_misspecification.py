"""Stress the Gaussian Wishart band under deliberate primitive misspecification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

from pilot_annni_block_hotelling import exact_half_moments, joint_sketch_covariance
from pilot_wishart_loewner import dual_certificate, wishart_band


def draw(rng: np.random.Generator, root: np.ndarray, count: int, law: str, dof: float) -> np.ndarray:
    dimension = root.shape[0]
    if law == "gaussian":
        primitive = rng.normal(size=(count, dimension))
    elif law == "student":
        primitive = rng.normal(size=(count, dimension))
        radial = np.sqrt(rng.chisquare(dof, size=count) / (dof - 2.0))
        primitive /= radial[:, None]
    elif law == "rademacher":
        primitive = rng.choice([-1.0, 1.0], size=(count, dimension))
    else:
        raise ValueError(law)
    return primitive @ root.T


def band_covers(sample_gram: np.ndarray, population: np.ndarray, alpha: float, count: int) -> bool:
    eta, _, _ = wishart_band(sample_gram, alpha, count)
    scale = np.sqrt(np.maximum(np.diag(population), np.finfo(float).tiny))
    inverse = np.diag(1.0 / scale)
    empirical = inverse @ sample_gram @ inverse
    truth = inverse @ population @ inverse
    generalized = eigvalsh(empirical, truth, check_finite=True)
    return bool(generalized[0] >= (1.0 - eta) ** 2 and generalized[-1] <= (1.0 + eta) ** 2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=8)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--sample-counts", nargs="+", type=int, default=[500, 1000])
    parser.add_argument("--gap-factors", nargs="+", type=float, default=[1.0, 1.15])
    parser.add_argument("--laws", nargs="+", choices=["gaussian", "student", "rademacher"], default=["gaussian", "student", "rademacher"])
    parser.add_argument("--student-dof", type=float, default=5.0)
    parser.add_argument("--ensembles", type=int, default=200)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--output", type=Path, default=Path("results/primitive_misspecification.json"))
    args = parser.parse_args()

    half, metadata = exact_half_moments(args.length, args.degree, 0.2, 1.0, 0.37, 2.2)
    joint = joint_sketch_covariance(half, 2 * args.degree + 1)
    values, vectors = np.linalg.eigh(joint)
    root = (vectors * np.sqrt(np.maximum(values, 0.0))) @ vectors.T
    rng = np.random.default_rng(args.seed)
    rows = []
    for law in args.laws:
        for count in args.sample_counts:
            for factor in args.gap_factors:
                theta = float(np.exp(-0.2 * factor * metadata["true_gap"]))
                detections = 0
                coverage = 0
                upper_values = []
                for _ in range(args.ensembles):
                    z = draw(rng, root, count, law, args.student_dof)
                    sample_gram = z.T @ z / count
                    coverage += band_covers(sample_gram, joint, args.alpha, count)
                    certificate = dual_certificate(sample_gram, args.degree, 2, theta, args.alpha, count)
                    upper_values.append(certificate["safe_certificate_upper"])
                    detections += certificate["safe_certificate_upper"] < -1e-9
                rows.append(
                    {
                        "law": law,
                        "sample_count": count,
                        "gap_factor": factor,
                        "coverage_rate": coverage / args.ensembles,
                        "detection_rate": detections / args.ensembles,
                        "safe_upper_range": [min(upper_values), max(upper_values)],
                    }
                )
                print(f"completed {law}, n={count}, factor={factor}", flush=True)
    payload = {
        "schema_version": 1,
        "setup": vars(args) | {"output": str(args.output)},
        "model": metadata,
        "rows": rows,
        "warning": "Only the Gaussian rows are covered by the theorem; the other laws are deliberate misspecification tests.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

