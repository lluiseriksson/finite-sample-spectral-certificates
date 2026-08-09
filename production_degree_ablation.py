"""Moment-degree ablation for the safe Wishart--Loewner certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from pilot_wishart_loewner import run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", nargs="+", type=int, default=[8, 12, 16])
    parser.add_argument("--degrees", nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--sample-counts", nargs="+", type=int, default=[500, 1000, 2000, 5000])
    parser.add_argument("--gap-factors", nargs="+", type=float, default=[1.0, 1.15])
    parser.add_argument("--ensembles", type=int, default=100)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--output", type=Path, default=Path("results/degree_ablation.json"))
    args = parser.parse_args()
    campaigns = []
    for length in args.lengths:
        for degree in args.degrees:
            child = SimpleNamespace(
                length=length,
                degree=degree,
                sample_counts=args.sample_counts,
                gap_factors=args.gap_factors,
                ensembles=args.ensembles,
                alpha=args.alpha,
                seed=args.seed + 100 * length + degree,
                tau=0.2,
                j1=1.0,
                j2=0.37,
                hx=2.2,
                output=args.output,
            )
            campaigns.append(run(child))
            print(f"completed L={length}, N={degree}", flush=True)
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

