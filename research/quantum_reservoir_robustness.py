#!/usr/bin/env python3
"""Passive-preserving component-tolerance audit for the six-port construction.

This is a sensitivity experiment, not part of the exact separation proof.  It
perturbs pole distances, arm phases, and balanced-interferometer angles while
keeping every sampled component lossless.  Performance is evaluated at the
nominal calibration nodes and across the frozen stop grid.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from research.quantum_reservoir_filter import blaschke_value, pass_nodes


DEFAULT_OUTPUT = ROOT / "results" / "quantum_reservoir" / "passive_tolerance.json"
DEFAULT_FIGURE = ROOT / "paper_quantum_reservoir" / "figures" / "passive_tolerance.pdf"


def perturbed_scalar(
    theta: np.ndarray,
    order: int,
    radius: float,
    phase: float,
    angle: float,
) -> np.ndarray:
    allpass = np.exp(1j * phase) * blaschke_value(theta, order, radius)
    return np.cos(angle) ** 2 + np.sin(angle) ** 2 * allpass


def audit(
    order: int,
    alpha: float,
    trials: int,
    pole_log_sigma: float,
    phase_sigma_units: float,
    beamsplitter_sigma: float,
    seed: int,
) -> dict[str, object]:
    nominal_radius = 1 - np.exp(-alpha * order)
    gap = 1 - nominal_radius
    nominal_phases = gap * np.array([-1.0, 0.0, 1.0])
    node_groups = [
        pass_nodes(order, nominal_radius, phase, np.pi / 6)
        for phase in nominal_phases
    ]
    stop_grid = np.concatenate(
        [np.linspace(-np.pi, -5 * np.pi / 6, 2001), np.linspace(5 * np.pi / 6, np.pi, 2001)]
    )
    rng = np.random.default_rng(seed)
    stop_norms = np.empty(trials)
    calibration_residuals = np.empty(trials)
    for trial in range(trials):
        radius_gaps = gap * np.exp(rng.normal(0.0, pole_log_sigma, size=3))
        radii = 1 - radius_gaps
        phases = nominal_phases + gap * rng.normal(0.0, phase_sigma_units, size=3)
        angles = np.pi / 4 + rng.normal(0.0, beamsplitter_sigma, size=3)
        stop_norms[trial] = max(
            float(np.max(np.abs(perturbed_scalar(stop_grid, order, radii[g], phases[g], angles[g]))))
            for g in range(3)
        )
        residual = 0.0
        for g, nodes in enumerate(node_groups):
            values = perturbed_scalar(nodes, order, radii[g], phases[g], angles[g])
            residual = max(residual, float(np.max(np.abs(values - 1))))
        calibration_residuals[trial] = residual

    def summary(values: np.ndarray) -> dict[str, float]:
        return {
            "median": float(np.quantile(values, 0.50)),
            "q90": float(np.quantile(values, 0.90)),
            "q95": float(np.quantile(values, 0.95)),
            "q99": float(np.quantile(values, 0.99)),
            "maximum": float(np.max(values)),
        }

    return {
        "schema_version": 1,
        "order": order,
        "alpha": alpha,
        "trials": trials,
        "seed": seed,
        "uncertainty_model": {
            "pole_gap_log_sigma": pole_log_sigma,
            "phase_sigma_in_nominal_gap_units": phase_sigma_units,
            "beamsplitter_angle_sigma_radians": beamsplitter_sigma,
            "passivity_preserved": True,
        },
        "stopband_signal_norm": summary(stop_norms),
        "nominal_node_calibration_residual": summary(calibration_residuals),
        "raw": {
            "stopband_signal_norm": stop_norms.tolist(),
            "nominal_node_calibration_residual": calibration_residuals.tolist(),
        },
    }


def make_figure(result: dict[str, object], output: Path) -> None:
    raw = result["raw"]
    stop = np.asarray(raw["stopband_signal_norm"], dtype=float)
    residual = np.asarray(raw["nominal_node_calibration_residual"], dtype=float)
    figure, axes = plt.subplots(1, 2, figsize=(6.6, 3.0))
    axes[0].hist(stop, bins=35, color="#1769aa", alpha=0.85)
    axes[0].axvline(np.quantile(stop, 0.95), color="#c44e52", linestyle="--", label="95th percentile")
    axes[0].set_xlabel("stopband signal norm")
    axes[0].set_ylabel("passive trials")
    axes[0].legend(fontsize=8)
    axes[1].scatter(residual, stop, s=7, alpha=0.35, color="#00a087", edgecolors="none")
    axes[1].set_xlabel("residual at nominal pass nodes")
    axes[1].set_ylabel("stopband signal norm")
    for axis in axes:
        axis.grid(True, alpha=0.22)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=11)
    parser.add_argument("--alpha", type=float, default=0.45)
    parser.add_argument("--trials", type=int, default=2000)
    parser.add_argument("--pole-log-sigma", type=float, default=0.01)
    parser.add_argument("--phase-sigma-units", type=float, default=0.01)
    parser.add_argument("--beamsplitter-sigma", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    args = parser.parse_args()
    result = audit(
        args.order,
        args.alpha,
        args.trials,
        args.pole_log_sigma,
        args.phase_sigma_units,
        args.beamsplitter_sigma,
        args.seed,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    make_figure(result, args.figure)
    print(json.dumps({key: value for key, value in result.items() if key != "raw"}, indent=2))


if __name__ == "__main__":
    main()
