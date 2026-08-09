"""Grid wrapper for the unknown-covariance Wishart-Loewner experiment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from pilot_wishart_loewner import run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", nargs="+", type=int, default=[6, 8, 10, 12, 14, 16])
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--sample-counts", nargs="+", type=int, default=[500, 1000, 2000, 5000])
    parser.add_argument("--gap-factors", nargs="+", type=float, default=[1.0, 1.15, 1.35])
    parser.add_argument("--ensembles", type=int, default=100)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--j1", type=float, default=1.0)
    parser.add_argument("--j2", type=float, default=0.37)
    parser.add_argument("--hx", type=float, default=2.2)
    parser.add_argument("--output", type=Path, default=Path("results/wishart_campaign.json"))
    args = parser.parse_args()

    campaigns = []
    for length in args.lengths:
        child = SimpleNamespace(
            length=length,
            degree=args.degree,
            sample_counts=args.sample_counts,
            gap_factors=args.gap_factors,
            ensembles=args.ensembles,
            alpha=args.alpha,
            seed=args.seed + length,
            tau=args.tau,
            j1=args.j1,
            j2=args.j2,
            hx=args.hx,
            output=args.output,
        )
        campaigns.append(run(child))
        print(f"completed L={length}", flush=True)

    payload = {
        "schema_version": 1,
        "setup": {key: (str(value) if isinstance(value, Path) else value) for key, value in vars(args).items()},
        "campaigns": campaigns,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

