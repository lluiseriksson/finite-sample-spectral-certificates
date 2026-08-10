#!/usr/bin/env python3
"""Audit the closed Markov pure-dephasing model, including vacuum ports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from research.quantum_reservoir_filter import lossless_block


DEFAULT_SCALING = ROOT / "results" / "quantum_reservoir" / "passive_filter_scaling.json"
DEFAULT_OUTPUT = ROOT / "results" / "quantum_reservoir" / "closed_dephasing.json"
DEFAULT_FIGURE = ROOT / "paper_quantum_reservoir" / "figures" / "closed_dephasing.pdf"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scaling", type=Path, default=DEFAULT_SCALING)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    args = parser.parse_args()

    scaling = json.loads(args.scaling.read_text(encoding="utf-8"))
    records = scaling["records"]
    baselines = [0.0, 1e-4, 1e-3, 1e-2]
    rate_records: list[dict[str, float | int | list[float]]] = []
    covariance_residual = 0.0
    alpha = float(scaling["alpha"])
    for record in records:
        order = int(record["order"])
        radius = float(record["radius"])
        gap = 1 - radius
        phases = gap * np.array([-1.0, 0.0, 1.0])
        for theta in np.linspace(-np.pi, np.pi, 257):
            scattering = lossless_block(theta, order, radius, phases)
            signal = scattering[:3, :3]
            loss = scattering[:3, 3:]
            covariance_residual = max(
                covariance_residual,
                float(np.linalg.norm(signal @ signal.conj().T + loss @ loss.conj().T - np.eye(3), 2)),
            )
        leakage = float(record["stopband_signal_norm"])
        advantages = [(baseline + 1.0) / (baseline + leakage**2) for baseline in baselines]
        rate_records.append(
            {
                "order": order,
                "stopband_signal_norm": leakage,
                "vacuum_baselines_over_j": baselines,
                "certified_total_rate_advantages": advantages,
            }
        )

    figure, axis = plt.subplots(figsize=(6.4, 3.8))
    orders = np.array([item["order"] for item in rate_records])
    for index, baseline in enumerate(baselines):
        values = np.array([item["certified_total_rate_advantages"][index] for item in rate_records])
        label = r"$\kappa_0/j=0$" if baseline == 0 else rf"$\kappa_0/j={baseline:g}$"
        axis.semilogy(orders, values, "o-", label=label)
    axis.set_xlabel("Blaschke order $S$")
    axis.set_ylabel("certified comparator / constructed total rate")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(ncol=2, fontsize=8)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, metadata={"CreationDate": None, "ModDate": None})
    plt.close(figure)

    result = {
        "schema_version": 1,
        "model": "six-input Markov Gaussian pure dephasing with three excess-noise inputs and three vacuum auxiliary inputs",
        "alpha": alpha,
        "vacuum_identity": "G G* + L L* = I, so V_out = I/2 + G J G*",
        "maximum_covariance_identity_residual": covariance_residual,
        "normalization": "j_minus = j_plus = 1; kappa_0 is the independently calibrated vacuum rate",
        "records": rate_records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
