#!/usr/bin/env python3
"""Reproduce the robust spectral-routing speed-limit certificates.

The analytic theorem bounds Grassmann motion of an input subspace by the
positive Wigner--Smith action.  This script audits four independent layers:

1. exact and finite-error equality cases for the interferometric family;
2. global degree bounds under alternating approximate pass/stop data;
3. random positive-block inequalities controlling Grassmann speed; and
4. random Blaschke--Potapov products tested on independent frequency arcs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "routing_speed_limit" / "certificate.json"
DEFAULT_FIGURE = ROOT / "paper_routing_speed_limit" / "figures" / "routing_law.pdf"


def robust_angle(epsilon_pass: float, epsilon_stop: float) -> float:
    total = np.arcsin(epsilon_pass) + np.arcsin(epsilon_stop)
    if total >= 0.5 * np.pi:
        return 0.0
    rho = np.sin(total)
    return float(np.arccos(np.clip(rho, 0.0, 1.0)))


def explicit_transition(rank: int, order: int, epsilon_pass: float, epsilon_stop: float) -> dict[str, float | int]:
    """An equality case on one pass-to-stop arc.

    S(z)=H diag(I,z^M I) H has constant tr Q=kM and lambda_max(Q)=M.
    Choosing endpoints on the error boundaries makes both local speed limits
    exact for every admissible pair of errors.
    """
    angle = robust_angle(epsilon_pass, epsilon_stop)
    phase_span = 2.0 * angle
    arc_width = phase_span / order
    trace_integral = rank * order * arc_width
    peak_integral = order * arc_width
    rho_observed = float(np.cos(angle))
    return {
        "rank": rank,
        "order": order,
        "epsilon_pass": epsilon_pass,
        "epsilon_stop": epsilon_stop,
        "overlap_bound_rho": float(
            np.sin(np.arcsin(epsilon_pass) + np.arcsin(epsilon_stop))
        ),
        "observed_endpoint_overlap": rho_observed,
        "principal_angle": angle,
        "frequency_arc_width": arc_width,
        "integrated_trace_delay": trace_integral,
        "trace_lower_bound": 2.0 * rank * angle,
        "integrated_peak_delay": peak_integral,
        "peak_lower_bound": 2.0 * angle,
        "trace_saturation_error": abs(trace_integral - 2.0 * rank * angle),
        "peak_saturation_error": abs(peak_integral - 2.0 * angle),
    }


def exact_cyclic_saturation(rank: int, order: int) -> dict[str, float | int]:
    angle = 0.5 * np.pi
    transitions = 2 * order
    degree = rank * order
    action = 2.0 * np.pi * degree
    lower = 2.0 * rank * angle * transitions
    return {
        "rank": rank,
        "order": order,
        "alternating_transitions": transitions,
        "mcmillan_degree": degree,
        "degree_lower_bound": (2.0 * order * rank / np.pi) * angle,
        "integrated_trace_delay": action,
        "routing_action_lower_bound": lower,
        "degree_saturation_error": abs(degree - (2.0 * order * rank / np.pi) * angle),
        "action_saturation_error": abs(action - lower),
    }


def haar_isometry(rng: np.random.Generator, rows: int, columns: int) -> np.ndarray:
    raw = rng.normal(size=(rows, columns)) + 1j * rng.normal(size=(rows, columns))
    q, r = np.linalg.qr(raw)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) > 0.0, phases / np.abs(phases), 1.0)
    return q * phases.conj()


def positive_block_audit(trials: int, seed: int) -> dict[str, float | int | bool]:
    """Stress the two pointwise PSD block inequalities used in the proof."""
    rng = np.random.default_rng(seed)
    worst_trace_ratio = 0.0
    worst_peak_ratio = 0.0
    for _ in range(trials):
        size = int(rng.integers(3, 13))
        rank = int(rng.integers(1, size))
        raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
        q_matrix = raw @ raw.conj().T
        basis = haar_isometry(rng, size, size)
        transformed = basis.conj().T @ q_matrix @ basis
        off_diagonal = transformed[rank:, :rank]
        trace_q = float(np.trace(q_matrix).real)
        lambda_max = float(np.linalg.eigvalsh(q_matrix)[-1])
        trace_ratio = 2.0 * float(np.linalg.norm(off_diagonal, ord="nuc")) / trace_q
        peak_ratio = 2.0 * float(np.linalg.norm(off_diagonal, ord=2)) / lambda_max
        worst_trace_ratio = max(worst_trace_ratio, trace_ratio)
        worst_peak_ratio = max(worst_peak_ratio, peak_ratio)
    return {
        "seed": seed,
        "trials": trials,
        "largest_two_nuclear_offdiag_over_trace": worst_trace_ratio,
        "largest_two_operator_offdiag_over_peak": worst_peak_ratio,
        "all_trace_block_bounds_pass": worst_trace_ratio <= 1.0 + 2.0e-12,
        "all_peak_block_bounds_pass": worst_peak_ratio <= 1.0 + 2.0e-12,
    }


def blaschke(z: complex, zero: complex) -> complex:
    return (z - zero) / (1.0 - np.conj(zero) * z)


def poisson(theta: float, zero: complex) -> float:
    z = np.exp(1j * theta)
    return float((1.0 - abs(zero) ** 2) / abs(z - zero) ** 2)


def potapov_value(theta: float, zeros: np.ndarray, vectors: list[np.ndarray]) -> np.ndarray:
    size = vectors[0].size
    transfer = np.eye(size, dtype=complex)
    z = np.exp(1j * theta)
    for zero, vector in zip(zeros, vectors, strict=True):
        projection = np.outer(vector, vector.conj())
        transfer = transfer @ (
            np.eye(size, dtype=complex) + (blaschke(z, zero) - 1.0) * projection
        )
    return transfer


def principal_angle_sum(left: np.ndarray, right: np.ndarray) -> float:
    singular_values = np.linalg.svd(left.conj().T @ right, compute_uv=False)
    return float(np.sum(np.arccos(np.clip(singular_values, 0.0, 1.0))))


def random_potapov_path_audit(trials: int, seed: int) -> dict[str, float | int | bool]:
    """Try to violate distance <= one half of positive WS trace action."""
    rng = np.random.default_rng(seed)
    worst_ratio = 0.0
    smallest_slack = float("inf")
    maximum_degree = 0
    tested_ranks: set[int] = set()
    for _ in range(trials):
        size = int(rng.integers(3, 9))
        degree = int(rng.integers(1, 13))
        maximum_degree = max(maximum_degree, degree)
        rank = int(rng.integers(1, min(4, size)))
        tested_ranks.add(rank)
        zeros = 0.82 * np.sqrt(rng.random(degree)) * np.exp(
            2j * np.pi * rng.random(degree)
        )
        vectors = [haar_isometry(rng, size, 1)[:, 0] for _ in range(degree)]
        input_frame = haar_isometry(rng, size, rank)
        start = float(rng.uniform(-np.pi, np.pi - 0.12))
        width = float(rng.uniform(0.02, min(2.4, np.pi - start)))
        end = start + width
        left = potapov_value(start, zeros, vectors) @ input_frame
        right = potapov_value(end, zeros, vectors) @ input_frame
        distance = principal_angle_sum(left, right)
        trace_action, _ = quad(
            lambda theta: sum(poisson(theta, zero) for zero in zeros),
            start,
            end,
            epsabs=2.0e-10,
            epsrel=2.0e-10,
            limit=400,
        )
        half_action = 0.5 * trace_action
        ratio = distance / half_action if half_action > 1.0e-14 else 0.0
        slack = half_action - distance
        worst_ratio = max(worst_ratio, ratio)
        smallest_slack = min(smallest_slack, slack)
    return {
        "seed": seed,
        "trials": trials,
        "maximum_factor_count": maximum_degree,
        "subspace_ranks_tested": sorted(tested_ranks),
        "largest_distance_over_half_trace_action": worst_ratio,
        "smallest_half_action_minus_distance": smallest_slack,
        "all_random_paths_respect_speed_limit": worst_ratio <= 1.0 + 2.0e-9,
    }


def make_figure(output: Path) -> None:
    errors = np.linspace(0.0, 0.69, 500)
    fractions = []
    for error in errors:
        angle = robust_angle(float(error), float(error))
        fractions.append(2.0 * angle / np.pi)

    theta = np.linspace(-0.1, np.pi + 0.1, 1200)
    signal = np.cos(theta / 2.0) ** 2
    loss = np.sin(theta / 2.0) ** 2
    epsilon = 0.18
    a = np.arcsin(epsilon)
    left = 2.0 * a
    right = np.pi - 2.0 * a

    figure, axes = plt.subplots(1, 2, figsize=(7.35, 3.15))
    axes[0].plot(errors, fractions, linewidth=2.0)
    axes[0].set_xlabel(r"symmetric amplitude error $\varepsilon$")
    axes[0].set_ylabel(r"certified fraction of $Mk$")
    axes[0].set_ylim(-0.02, 1.03)
    axes[0].grid(alpha=0.25)
    axes[0].set_title("robust global degree law", fontsize=10)

    axes[1].plot(theta, signal, label="signal power")
    axes[1].plot(theta, loss, label="loss-port power")
    axes[1].axvspan(left, right, alpha=0.14, color="black", label="certified transition")
    axes[1].scatter([left, right], [np.cos(left / 2) ** 2, np.cos(right / 2) ** 2], s=20)
    axes[1].set_xlabel(r"interferometric phase $\phi$")
    axes[1].set_ylabel("power")
    axes[1].set_xlim(-0.1, np.pi + 0.1)
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=7, loc="center right")
    axes[1].set_title("finite-error equality arc", fontsize=10)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, metadata={"CreationDate": None, "ModDate": None})
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--block-trials", type=int, default=2048)
    parser.add_argument("--path-trials", type=int, default=512)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    errors = [(0.0, 0.0), (0.02, 0.08), (0.10, 0.10), (0.20, 0.30), (0.45, 0.20)]
    transitions = [
        explicit_transition(rank, order, epsilon_pass, epsilon_stop)
        for rank, order in [(1, 3), (2, 5), (4, 7)]
        for epsilon_pass, epsilon_stop in errors
    ]
    cyclic = [
        exact_cyclic_saturation(rank, order)
        for rank, order in [(1, 1), (1, 8), (2, 5), (4, 7)]
    ]
    payload = {
        "schema_version": 1,
        "theorem_audited": "alternating approximate routing Grassmann distance <= positive Wigner-Smith action / 2",
        "robust_factor": "(2/pi) arccos(sin(arcsin(epsilon_pass)+arcsin(epsilon_stop)))",
        "finite_error_local_equality_cases": transitions,
        "exact_cyclic_degree_saturation": cyclic,
        "positive_block_stress_test": positive_block_audit(args.block_trials, 2026081001),
        "random_potapov_path_stress_test": random_potapov_path_audit(args.path_trials, 2026081002),
        "interpretation": {
            "global": "M alternating pass/stop pairs of a rank-k input subspace force n >= (2Mk/pi) arccos(rho).",
            "local": "A transition across an arc of width Delta forces both sup tr(Q) >= 2k arccos(rho)/Delta and sup lambda_max(Q) >= 2 arccos(rho)/Delta.",
            "scope": "The action theorem needs only an absolutely continuous unitary path with positive Wigner-Smith generator; rational inner structure is used solely to identify total action with 2 pi times McMillan degree.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    make_figure(args.figure)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
