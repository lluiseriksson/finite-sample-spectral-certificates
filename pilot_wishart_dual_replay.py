"""Replay primal classifications with explicit Wishart--Loewner dual witnesses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pilot_annni_block_hotelling import exact_half_moments, joint_sketch_covariance
from pilot_wishart_loewner import compatible, dual_certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=8)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--ensembles", type=int, default=100)
    parser.add_argument("--gap-factors", nargs="+", type=float, default=[1.0, 1.35])
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--output", type=Path, default=Path("results/pilots/wishart_dual_replay.json"))
    args = parser.parse_args()

    tau, j1, j2, hx = 0.2, 1.0, 0.37, 2.2
    half, metadata = exact_half_moments(args.length, args.degree, tau, j1, j2, hx)
    joint = joint_sketch_covariance(half, 2 * args.degree + 1)
    rng = np.random.default_rng(args.seed)
    rows = []
    for factor in args.gap_factors:
        theta = float(np.exp(-tau * factor * metadata["true_gap"]))
        for ensemble in range(args.ensembles):
            z = rng.multivariate_normal(np.zeros(joint.shape[0]), joint, size=args.samples)
            sample_gram = z.T @ z / args.samples
            primal_compatible, primal_status = compatible(
                sample_gram, args.degree, 2, theta, args.alpha, args.samples
            )
            dual = dual_certificate(sample_gram, args.degree, 2, theta, args.alpha, args.samples)
            rows.append(
                {
                    "gap_factor": factor,
                    "ensemble": ensemble,
                    "primal_compatible": primal_compatible,
                    "primal_status": primal_status,
                    **dual,
                }
            )

    for row in rows:
        row["dual_detected_at_1e-7"] = row["dual_value"] < -1e-7
        row["sign_agrees"] = (not row["primal_compatible"]) == row["dual_detected_at_1e-7"]
    payload = {
        "schema_version": 1,
        "setup": vars(args) | {"output": str(args.output)},
        "model": metadata,
        "summary": {
            "trials": len(rows),
            "sign_disagreements": sum(not row["sign_agrees"] for row in rows),
            "maximum_stationarity_fro": max(row["dual_stationarity_fro"] for row in rows),
            "maximum_normalization_error": max(row["dual_normalization_error"] for row in rows),
            "minimum_reported_psd_eigenvalue": min(row["dual_min_psd_eigenvalue"] for row in rows),
            "minimum_null_dual_value": min(
                row["dual_value"] for row in rows if row["gap_factor"] == 1.0
            ),
            "maximum_alternative_dual_value": max(
                row["dual_value"] for row in rows if row["gap_factor"] > 1.0
            ),
        },
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
