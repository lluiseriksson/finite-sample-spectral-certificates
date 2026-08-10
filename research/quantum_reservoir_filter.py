#!/usr/bin/env python3
"""Construct and audit the passive multichannel reservoir-filter family.

The signal block is a rational Schur function embedded explicitly in a
six-port lossless network.  The construction combines three Mach-Zehnder loss
ports with a fixed three-channel paraunitary mixer.  It prints a JSON ledger
and optionally writes the ledger and a scaling figure.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "quantum_reservoir" / "passive_filter_scaling.json"
DEFAULT_FIGURE = ROOT / "paper_quantum_reservoir" / "figures" / "passive_scaling.pdf"
DEFAULT_ARCHITECTURE = ROOT / "paper_quantum_reservoir" / "figures" / "passive_architecture.pdf"


def dft3() -> np.ndarray:
    omega = np.exp(2j * np.pi / 3)
    return np.array([[omega ** (row * col) for col in range(3)] for row in range(3)]) / np.sqrt(3)


def boundary_phase(theta: np.ndarray | float, radius: float) -> np.ndarray | float:
    """Unwrapped phase of (z-r)/(1-rz) on -pi < theta < pi."""
    scale = (1 + radius) / (1 - radius)
    return 2 * np.arctan(scale * np.tan(np.asarray(theta) / 2))


def blaschke_value(theta: np.ndarray, order: int, radius: float) -> np.ndarray:
    z = np.exp(1j * theta)
    base = (z - radius) / (1 - radius * z)
    normalization = 1 if order % 2 else -1
    return normalization * base**order


def pass_nodes(order: int, radius: float, phase: float, half_width: float) -> np.ndarray:
    normalization_phase = 0.0 if order % 2 else np.pi
    lower = float(boundary_phase(-half_width, radius))
    upper = float(boundary_phase(half_width, radius))
    candidates = []
    for winding in range(-2 * order, 2 * order + 1):
        target = (2 * np.pi * winding - normalization_phase - phase) / order
        if lower < target < upper:
            theta = 2 * np.arctan(
                ((1 - radius) / (1 + radius)) * np.tan(target / 2)
            )
            candidates.append(float(theta))
    return np.array(sorted(set(candidates)))


def mixer(z: complex) -> np.ndarray:
    fourier = dft3()
    return fourier @ np.diag([1, z, z * z]) @ fourier.conj().T


def reverse_mixer(z: complex) -> np.ndarray:
    fourier = dft3()
    return fourier @ np.diag([z * z, z, 1]) @ fourier.conj().T


def signal_block(theta: float, order: int, radius: float, phases: np.ndarray) -> np.ndarray:
    z = np.exp(1j * theta)
    b_value = blaschke_value(np.array([theta]), order, radius)[0]
    diagonal = np.diag((1 + np.exp(1j * phases) * b_value) / 2)
    return mixer(z) @ diagonal @ reverse_mixer(z)


def lossless_block(theta: float, order: int, radius: float, phases: np.ndarray) -> np.ndarray:
    """Return the full 6x6 unitary scattering matrix at one frequency."""
    z = np.exp(1j * theta)
    b_value = blaschke_value(np.array([theta]), order, radius)[0]
    interferometer = np.zeros((6, 6), dtype=complex)
    for channel, phase in enumerate(phases):
        value = np.exp(1j * phase) * b_value
        block = 0.5 * np.array([[1 + value, 1 - value], [1 - value, 1 + value]])
        indices = [channel, channel + 3]
        interferometer[np.ix_(indices, indices)] = block
    pre = np.block(
        [[reverse_mixer(z), np.zeros((3, 3))], [np.zeros((3, 3)), np.eye(3)]]
    )
    post = np.block([[mixer(z), np.zeros((3, 3))], [np.zeros((3, 3)), np.eye(3)]])
    return post @ interferometer @ pre


def audit_order(order: int, alpha: float) -> dict[str, float | int | bool]:
    if order < 5:
        raise ValueError("order must be at least five for the root-count separation")
    radius = 1 - np.exp(-alpha * order)
    phase_scale = 1 - radius
    phases = phase_scale * np.array([-1.0, 0.0, 1.0])
    pass_half_width = np.pi / 6
    groups = [pass_nodes(order, radius, phase, pass_half_width) for phase in phases]
    nodes_per_group = min(len(group) for group in groups)
    groups = [group[:nodes_per_group] for group in groups]

    signatures: list[np.ndarray] = []
    calibration_residual = 0.0
    for channel, group in enumerate(groups):
        for theta in group:
            z = np.exp(1j * theta)
            vector = mixer(z)[:, channel]
            signatures.append(vector)
            calibration_residual = max(
                calibration_residual,
                float(np.linalg.norm(signal_block(theta, order, radius, phases) @ vector - z**2 * vector)),
            )

    minimum_minor = min(
        abs(np.linalg.det(np.column_stack([signatures[index] for index in triple])))
        for triple in combinations(range(len(signatures)), 3)
    )

    stop_grid = np.concatenate(
        [np.linspace(-np.pi, -5 * np.pi / 6, 4001), np.linspace(5 * np.pi / 6, np.pi, 4001)]
    )
    global_grid = np.linspace(-np.pi, np.pi, 12001)
    stop_norm = max(
        np.linalg.svd(signal_block(theta, order, radius, phases), compute_uv=False)[0]
        for theta in stop_grid
    )
    global_norm = max(
        np.linalg.svd(signal_block(theta, order, radius, phases), compute_uv=False)[0]
        for theta in global_grid
    )
    unitary_error = max(
        np.linalg.norm(
            lossless_block(theta, order, radius, phases).conj().T
            @ lossless_block(theta, order, radius, phases)
            - np.eye(6),
            ord=2,
        )
        for theta in np.linspace(-np.pi, np.pi, 1001)
    )
    peak_group_delay = order * (1 + radius) / (1 - radius)
    numerator_degree = order + 4
    root_count_margin = len(signatures) - 2 - numerator_degree
    return {
        "order": order,
        "radius": radius,
        "nodes_per_group": nodes_per_group,
        "calibration_count": len(signatures),
        "numerator_degree": numerator_degree,
        "root_count_margin": root_count_margin,
        "full_spark_minimum_minor": float(minimum_minor),
        "maximum_calibration_residual": calibration_residual,
        "global_signal_norm": float(global_norm),
        "stopband_signal_norm": float(stop_norm),
        "worst_case_rate_ratio_upper": float(stop_norm**2),
        "reducible_stopband_lower_bound": 1.0,
        "rate_separation_factor_lower": float(1 / stop_norm**2),
        "peak_allpass_group_delay": float(peak_group_delay),
        "six_port_unitarity_error": float(unitary_error),
        "root_count_separation_valid": root_count_margin > 0,
    }


def make_figure(records: list[dict[str, float | int | bool]], output: Path) -> None:
    orders = np.array([record["order"] for record in records], dtype=float)
    leakage = np.array([record["stopband_signal_norm"] for record in records], dtype=float)
    delay = np.array([record["peak_allpass_group_delay"] for record in records], dtype=float)
    figure, left = plt.subplots(figsize=(6.4, 3.8))
    right = left.twinx()
    left.semilogy(orders, leakage, "o-", color="#1769aa", label="signal leakage")
    left.semilogy(orders, leakage**2, "s--", color="#00a087", label="rate upper bound")
    right.semilogy(orders, delay, "^-", color="#c44e52", label="peak group delay")
    left.set_xlabel("Blaschke order S")
    left.set_ylabel("stopband leakage / rate ratio")
    right.set_ylabel("peak all-pass group delay")
    left.grid(True, which="both", alpha=0.25)
    handles_left, labels_left = left.get_legend_handles_labels()
    handles_right, labels_right = right.get_legend_handles_labels()
    left.legend(handles_left + handles_right, labels_left + labels_right, loc="center left")
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)


def make_architecture_figure(output: Path) -> None:
    figure, axis = plt.subplots(figsize=(7.0, 2.8))
    axis.set_xlim(0, 10)
    axis.set_ylim(-0.2, 4.3)
    axis.axis("off")
    box_style = dict(boxstyle="round,pad=0.25", facecolor="#eef4fb", edgecolor="#1769aa", linewidth=1.2)
    loss_style = dict(boxstyle="round,pad=0.22", facecolor="#f8eeee", edgecolor="#c44e52", linewidth=1.0)
    axis.text(0.35, 2.15, "3 bath\nsignals", ha="center", va="center")
    axis.text(1.8, 2.15, "$U_{\\rm rev}(z)$\nDFT + delays", ha="center", va="center", bbox=box_style)
    axis.text(8.2, 2.15, "$U(z)$\nDFT + delays", ha="center", va="center", bbox=box_style)
    axis.text(9.65, 2.15, "3 filtered\nsignals", ha="center", va="center")
    for y, channel in zip([3.35, 2.15, 0.95], [0, 1, 2], strict=True):
        axis.text(5.0, y, f"$W_{channel}(z)=H\\,\\mathrm{{diag}}(1,B_{channel})\\,H$", ha="center", va="center", bbox=box_style)
        axis.text(5.0, y - 0.55, "vacuum/loss port", ha="center", va="center", fontsize=8, color="#8b2e2e", bbox=loss_style)
        axis.annotate("", xy=(6.8, y), xytext=(3.2, y), arrowprops=dict(arrowstyle="->", color="#555555"))
        axis.annotate("", xy=(5.0, y - 0.28), xytext=(5.0, y - 0.47), arrowprops=dict(arrowstyle="->", color="#c44e52"))
    axis.annotate("", xy=(1.15, 2.15), xytext=(0.75, 2.15), arrowprops=dict(arrowstyle="->", linewidth=1.4))
    axis.annotate("", xy=(3.05, 2.15), xytext=(2.45, 2.15), arrowprops=dict(arrowstyle="->", linewidth=1.4))
    axis.annotate("", xy=(7.55, 2.15), xytext=(6.95, 2.15), arrowprops=dict(arrowstyle="->", linewidth=1.4))
    axis.annotate("", xy=(9.25, 2.15), xytext=(8.85, 2.15), arrowprops=dict(arrowstyle="->", linewidth=1.4))
    axis.text(5.0, 4.05, "explicit causal inner six-port completion", ha="center", va="center", weight="bold")
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, default=0.45)
    parser.add_argument("--orders", type=int, nargs="+", default=[5, 7, 9, 11, 15, 19])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--architecture-figure", type=Path, default=DEFAULT_ARCHITECTURE)
    args = parser.parse_args()

    records = [audit_order(order, args.alpha) for order in args.orders]
    result = {
        "schema_version": 1,
        "construction": "three-channel Blaschke--Mach-Zehnder passive reservoir filter",
        "alpha": args.alpha,
        "pass_arc": [-np.pi / 6, np.pi / 6],
        "stop_arcs": [[-np.pi, -5 * np.pi / 6], [5 * np.pi / 6, np.pi]],
        "records": records,
        "interpretation": {
            "proved_by_construction": "The signal block is globally contractive and is embedded in an explicitly unitary six-port scattering matrix. Exact delayed tangential calibrations are full spark, and the rational root count excludes every constant-channel reducing architecture of the stated bidegree.",
            "resource_caveat": "Exponential stopband suppression uses poles exponentially close to the unit circle and therefore incurs exponential peak group delay; it moves the resource boundary rather than removing it.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    make_figure(records, args.figure)
    make_architecture_figure(args.architecture_figure)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
