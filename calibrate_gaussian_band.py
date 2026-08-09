"""Monte Carlo audit of conservatism in the analytic Gaussian Loewner band."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def scores(rng: np.random.Generator, ensembles: int, count: int, dimension: int, batch: int) -> np.ndarray:
    result = np.empty(ensembles)
    offset = 0
    while offset < ensembles:
        take = min(batch, ensembles - offset)
        g = rng.normal(size=(take, count, dimension))
        gram = np.einsum("bni,bnj->bij", g, g) / count
        eigenvalues = np.linalg.eigvalsh(gram)
        singular_min = np.sqrt(eigenvalues[:, 0])
        singular_max = np.sqrt(eigenvalues[:, -1])
        result[offset : offset + take] = np.maximum(1.0 - singular_min, singular_max - 1.0)
        offset += take
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, default=12)
    parser.add_argument("--sample-counts", nargs="+", type=int, default=[500, 1000, 2000, 5000])
    parser.add_argument("--ensembles", type=int, default=20000)
    parser.add_argument("--batch", type=int, default=500)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260819)
    parser.add_argument("--output", type=Path, default=Path("results/calibration/gaussian_band.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = []
    for count in args.sample_counts:
        values = scores(rng, args.ensembles, count, args.dimension, args.batch)
        analytic = (np.sqrt(args.dimension) + np.sqrt(2.0 * np.log(2.0 / args.alpha))) / np.sqrt(count)
        rows.append(
            {
                "sample_count": count,
                "analytic_eta": float(analytic),
                "analytic_empirical_coverage": float(np.mean(values <= analytic)),
                "score_quantiles": {
                    "0.90": float(np.quantile(values, 0.90)),
                    "0.95": float(np.quantile(values, 0.95)),
                    "0.99": float(np.quantile(values, 0.99)),
                    "0.999": float(np.quantile(values, 0.999)),
                },
                "analytic_to_empirical_95_width_ratio": float(analytic / np.quantile(values, 0.95)),
            }
        )
        print(f"completed n={count}", flush=True)
    payload = {"schema_version": 1, "setup": vars(args) | {"output": str(args.output)}, "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
