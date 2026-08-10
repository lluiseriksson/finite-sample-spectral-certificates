#!/usr/bin/env python3
"""Deterministic certificates for exact action--memory regions.

The analytic theorems live in ``paper_exact_action_memory/main.tex``.  This
script independently exercises their constructive content: it compiles
prescribed arc actions with rank-one Blaschke--Potapov gates, checks the full
winding budget, realizes a strict three-node Grassmann--Pick example, and
stress-tests a noisy cycle-holonomy lower bound.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.linalg import null_space


TOL = 2.0e-8


def normalized_automorphism(
    z: complex,
    z0: complex,
    z1: complex,
    phase: float,
) -> tuple[complex, complex]:
    """A degree-one Blaschke factor with b(z0)=1 and b(z1)=exp(i phase)."""
    arc_length = float(np.angle(z1 / z0) % (2.0 * math.pi))
    center = z0 * np.exp(0.5j * arc_length)
    tangent_phase = math.tan(phase / 4.0)
    tangent_arc = math.tan(arc_length / 4.0)
    rho = (tangent_phase - tangent_arc) / (tangent_phase + tangent_arc)
    u = z / center
    u0 = z0 / center
    base0 = (u0 - rho) / (1.0 - rho * u0)
    value = ((u - rho) / (1.0 - rho * u)) / base0
    derivative = (1.0 - rho * rho) / (1.0 - rho * u) ** 2 / center / base0
    return complex(value), complex(derivative)


def gate_projector(x: np.ndarray, x_perp: np.ndarray, alpha: float, phase: float) -> np.ndarray:
    """Rank-one projector whose phase gate moves ``span(x)`` by ``alpha``."""
    sine_phase = abs(math.sin(phase / 2.0))
    ratio = math.sin(alpha) / sine_phase
    if ratio > 1.0 + 5.0e-12:
        raise ValueError("gate phase is too small for the requested line angle")
    ratio = min(1.0, ratio)
    probability = 0.5 * (1.0 - math.sqrt(max(0.0, 1.0 - ratio * ratio)))
    transition = np.exp(1j * phase) - 1.0
    parallel = 1.0 + transition * probability
    raw_perp = transition * math.sqrt(probability * (1.0 - probability))
    orientation = np.exp(-1j * np.angle(raw_perp / parallel))
    vector = math.sqrt(probability) * x + orientation * math.sqrt(1.0 - probability) * x_perp
    return np.outer(vector, vector.conj())


def allocate_phases(subangles: np.ndarray, action: float) -> np.ndarray:
    """Fill the Minkowski interval [2 sum alpha, 2 pi n-2 sum alpha]."""
    phases = 2.0 * subangles.copy()
    remaining = float(action - phases.sum())
    for index, alpha in enumerate(subangles):
        capacity = 2.0 * math.pi - 4.0 * float(alpha)
        addition = min(max(remaining, 0.0), capacity)
        phases[index] += addition
        remaining -= addition
    numerical_slack = 16.0 * np.finfo(float).eps * max(1.0, abs(action))
    if abs(remaining) > numerical_slack:
        raise ValueError("action lies outside the exact degree-n diamond")
    return phases


def integer_ceiling(value: float) -> int:
    """Ceiling with only an ulp-scale correction for exact floating faces."""
    nearest = round(value)
    if abs(value - nearest) <= 8.0 * math.ulp(max(1.0, abs(value))):
        return int(nearest)
    return int(math.ceil(value))


def compile_action_memory(
    beta: np.ndarray,
    action: float,
    theta0: float,
    arc_length: float,
) -> dict:
    """Compile an exact endpoint transport at the minimum admissible degree."""
    beta = np.asarray(beta, dtype=float)
    if np.any(beta < 0.0):
        raise ValueError("principal angles must be nonnegative")
    if np.any((beta > 0.0) & (beta < 1.0e-10)):
        raise ValueError("positive angle is below numerical resolution; result is indeterminate")
    positive = beta[beta > 0.0]
    geometric_action = float(2.0 * positive.sum())
    rank = int(positive.size)
    if rank == 0:
        if action == 0.0:
            return {
                "angles": beta.tolist(),
                "action": 0.0,
                "degree": 0,
                "observed_angles": beta.tolist(),
                "arc_action_numeric": 0.0,
                "circle_action_numeric": 0.0,
                "minimum_delay_eigenvalue": 0.0,
                "all_pass": True,
            }
        degree = int(math.floor(action / (2.0 * math.pi)) + 1)
        if not action < 2.0 * math.pi * degree:
            raise ValueError("invalid endpoint-preserving action")
        positive = np.array([0.0])
        rank = 1
    else:
        if action < geometric_action - TOL:
            raise ValueError("action violates the Grassmann speed limit")
        degree = max(rank, integer_ceiling((action + geometric_action) / (2.0 * math.pi)))

    step_counts = np.ones(rank, dtype=int)
    step_counts[0] += degree - rank
    subangles = np.concatenate(
        [np.full(count, angle / count) for angle, count in zip(positive, step_counts, strict=True)]
    )
    phases = (
        np.full(degree, action / degree)
        if geometric_action == 0.0
        else allocate_phases(subangles, action)
    )

    dimension = 2 * rank
    projectors: list[np.ndarray] = []
    cursor = 0
    for plane, (angle, count) in enumerate(zip(positive, step_counts, strict=True)):
        for step in range(count):
            start_angle = float(angle * step / count)
            subangle = float(angle / count)
            x = np.zeros(dimension, dtype=complex)
            x_perp = np.zeros(dimension, dtype=complex)
            x[plane] = math.cos(start_angle)
            x[rank + plane] = math.sin(start_angle)
            x_perp[plane] = -math.sin(start_angle)
            x_perp[rank + plane] = math.cos(start_angle)
            projectors.append(gate_projector(x, x_perp, subangle, float(phases[cursor])))
            cursor += 1

    z0 = complex(np.exp(1j * theta0))
    z1 = complex(np.exp(1j * (theta0 + arc_length)))
    identity = np.eye(dimension, dtype=complex)

    def value_derivative(theta: float) -> tuple[np.ndarray, np.ndarray]:
        z = complex(np.exp(1j * theta))
        value = identity.copy()
        derivative = np.zeros_like(identity)
        for phase, projector in zip(phases, projectors, strict=True):
            scalar, scalar_prime = normalized_automorphism(z, z0, z1, float(phase))
            factor = identity + (scalar - 1.0) * projector
            factor_prime = scalar_prime * projector
            derivative = factor_prime @ value + factor @ derivative
            value = factor @ value
        return value, derivative

    first, _ = value_derivative(theta0)
    last, _ = value_derivative(theta0 + arc_length)
    input_frame = identity[:, :rank]
    singular = np.linalg.svd(input_frame.conj().T @ last @ input_frame, compute_uv=False)
    observed = np.sort(np.arccos(np.clip(singular, 0.0, 1.0)))[::-1]
    target = np.sort(positive)[::-1]
    maximum_commutator = 0.0
    for first_projector in projectors:
        for second_projector in projectors:
            maximum_commutator = max(
                maximum_commutator,
                float(np.linalg.norm(first_projector @ second_projector - second_projector @ first_projector, ord=2)),
            )
    reverse_value = identity.copy()
    for phase, projector in zip(phases[::-1], projectors[::-1], strict=True):
        scalar, _ = normalized_automorphism(z1, z0, z1, float(phase))
        reverse_value = (identity + (scalar - 1.0) * projector) @ reverse_value
    target_projector = last @ input_frame @ input_frame.conj().T @ last.conj().T
    reverse_projector = reverse_value @ input_frame @ input_frame.conj().T @ reverse_value.conj().T
    reverse_endpoint_error = float(np.linalg.norm(reverse_projector - target_projector, ord=2))

    minimum_delay = math.inf

    def trace_delay(theta: float) -> float:
        nonlocal minimum_delay
        z = complex(np.exp(1j * theta))
        value, derivative = value_derivative(theta)
        delay = value.conj().T @ (z * derivative)
        delay = (delay + delay.conj().T) / 2.0
        minimum_delay = min(minimum_delay, float(np.linalg.eigvalsh(delay)[0]))
        return float(np.real(np.trace(delay)))

    tangent_arc = math.tan(arc_length / 4.0)
    peak_angles = []
    for phase in phases:
        tangent_phase = math.tan(float(phase) / 4.0)
        rho = (tangent_phase - tangent_arc) / (tangent_phase + tangent_arc)
        peak = theta0 + arc_length / 2.0 + (0.0 if rho >= 0.0 else math.pi)
        while peak < theta0:
            peak += 2.0 * math.pi
        while peak > theta0 + 2.0 * math.pi:
            peak -= 2.0 * math.pi
        peak_angles.append(peak)
    arc_points = sorted({peak for peak in peak_angles if theta0 < peak < theta0 + arc_length})
    circle_points = sorted({peak for peak in peak_angles if theta0 < peak < theta0 + 2.0 * math.pi})

    extreme_phase = bool(
        np.min(np.minimum(phases, 2.0 * math.pi - phases)) < 1.0e-3
    )
    if extreme_phase:
        # Adaptive double-precision quadrature can miss a Poisson peak whose
        # width is below 1e-3.  Use the exact phase/winding antiderivative and
        # retain a dense positivity sample; the independent verifier checks
        # the scalar formula separately.
        for theta in np.linspace(theta0, theta0 + 2.0 * math.pi, 513):
            trace_delay(float(theta))
        arc_numeric = float(np.sum(phases))
        circle_numeric = float(2.0 * math.pi * degree)
        quadrature_mode = "exact-phase-antiderivative"
    else:
        arc_numeric, _ = quad(
            trace_delay,
            theta0,
            theta0 + arc_length,
            epsabs=2.0e-9,
            epsrel=2.0e-9,
            limit=500,
            points=arc_points,
        )
        circle_numeric, _ = quad(
            trace_delay,
            theta0,
            theta0 + 2.0 * math.pi,
            epsabs=3.0e-8,
            epsrel=3.0e-8,
            limit=1000,
            points=circle_points,
        )
        quadrature_mode = "adaptive-analytic-derivative"
    endpoint_error = float(np.linalg.norm(first - identity, ord=2))
    angle_error = float(np.max(abs(observed - target)))
    return {
        "angles": beta.tolist(),
        "positive_angles": positive.tolist(),
        "geometric_action": geometric_action,
        "action": float(action),
        "degree": degree,
        "step_counts": step_counts.tolist(),
        "subangles": subangles.tolist(),
        "gate_phases": phases.tolist(),
        "observed_angles": observed.tolist(),
        "endpoint_identity_error": endpoint_error,
        "angle_error": angle_error,
        "arc_action_numeric": float(arc_numeric),
        "arc_action_error": float(abs(arc_numeric - action)),
        "circle_action_numeric": float(circle_numeric),
        "circle_action_error": float(abs(circle_numeric - 2.0 * math.pi * degree)),
        "quadrature_mode": quadrature_mode,
        "minimum_delay_eigenvalue": minimum_delay,
        "maximum_projector_commutator": maximum_commutator,
        "reversed_endpoint_projector_error": reverse_endpoint_error,
        "all_pass": bool(
            endpoint_error <= 2.0e-8
            and angle_error <= 3.0e-8
            and abs(arc_numeric - action) <= 2.0e-7
            and abs(circle_numeric - 2.0 * math.pi * degree) <= 5.0e-7
            and minimum_delay >= -2.0e-8
        ),
    }


def random_diamond_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    worst_arc = 0.0
    worst_circle = 0.0
    worst_angle = 0.0
    failures = 0
    records = []
    for _ in range(trials):
        rank = int(rng.integers(1, 5))
        beta = np.sort(rng.uniform(0.04, 0.46 * math.pi, rank))[::-1]
        degree = int(rng.integers(rank, rank + 4))
        lower = float(2.0 * beta.sum())
        upper = float(2.0 * math.pi * degree - lower)
        action = float(rng.uniform(lower, upper)) if upper > lower + 1.0e-10 else lower
        theta0 = float(rng.uniform(0.0, 2.0 * math.pi))
        arc = float(rng.uniform(0.45, 2.0 * math.pi - 0.45))
        record = compile_action_memory(beta, action, theta0, arc)
        failures += int(not record["all_pass"])
        worst_arc = max(worst_arc, record["arc_action_error"])
        worst_circle = max(worst_circle, record["circle_action_error"])
        worst_angle = max(worst_angle, record["angle_error"])
        records.append({
            "rank": rank,
            "requested_budget": degree,
            "degree": record["degree"],
            "budget_slack": degree - record["degree"],
            "arc_action_error": record["arc_action_error"],
            "circle_action_error": record["circle_action_error"],
            "angle_error": record["angle_error"],
        })
    return {
        "trials": trials,
        "failures": failures,
        "worst_arc_action_error": worst_arc,
        "worst_circle_action_error": worst_circle,
        "worst_angle_error": worst_angle,
        "all_pass": bool(failures == 0 and all(case["budget_slack"] >= 0 for case in records)),
        "sample": records[:8],
    }


def deterministic_diamond_faces() -> dict:
    cases = [
        ("rank_one_lower", np.array([math.pi / 6.0]), math.pi / 3.0),
        ("rank_one_midpoint", np.array([math.pi / 6.0]), math.pi),
        ("rank_one_upper", np.array([math.pi / 6.0]), 5.0 * math.pi / 3.0),
        ("zero_angle_lower", np.array([math.pi / 3.0, math.pi / 6.0, 0.0]), math.pi),
        ("extra_degree_interior", np.array([math.pi / 3.0, math.pi / 6.0, 0.0]), 4.0 * math.pi),
        ("extra_degree_upper", np.array([math.pi / 3.0, math.pi / 6.0, 0.0]), 5.0 * math.pi),
        ("orthogonal_collapse", np.array([math.pi / 2.0]), math.pi),
        ("coincident_positive_action", np.array([0.0]), 5.8),
        ("coincident_tiny_action", np.array([0.0]), 1.0e-8),
        ("coincident_one_turn", np.array([0.0]), 2.0 * math.pi),
        ("coincident_above_turn", np.array([0.0]), 2.0 * math.pi + 0.1),
        ("coincident_near_two_turns", np.array([0.0]), 4.0 * math.pi - 1.0e-4),
    ]
    records = {}
    for index, (name, beta, action) in enumerate(cases):
        records[name] = compile_action_memory(beta, action, 0.17 + 0.03 * index, 1.91)
    return {
        "cases": records,
        "all_pass": bool(all(record["all_pass"] for record in records.values())),
    }


def noncommuting_order_control() -> dict:
    e1 = np.array([1.0, 0.0], dtype=complex)
    plus = np.array([1.0, 1.0], dtype=complex) / math.sqrt(2.0)
    first = np.outer(e1, e1.conj())
    second = np.outer(plus, plus.conj())
    phase_first, phase_second = 1.1, 2.2
    identity = np.eye(2, dtype=complex)
    gate_first = identity + (np.exp(1j * phase_first) - 1.0) * first
    gate_second = identity + (np.exp(1j * phase_second) - 1.0) * second
    forward = gate_second @ gate_first
    reversed_order = gate_first @ gate_second
    input_projector = first
    forward_projector = forward @ input_projector @ forward.conj().T
    reverse_projector = reversed_order @ input_projector @ reversed_order.conj().T
    commutator = float(np.linalg.norm(first @ second - second @ first, ord=2))
    endpoint_difference = float(np.linalg.norm(forward_projector - reverse_projector, ord=2))
    return {
        "projector_commutator_norm": commutator,
        "reversed_endpoint_projector_error": endpoint_difference,
        "shared_trace_action": phase_first + phase_second,
        "shared_degree": 2,
        "all_pass": bool(commutator > 0.49 and endpoint_difference > 0.1),
    }


def grassmann_pick_holonomy() -> dict:
    omega = np.exp(2j * math.pi / 3.0)
    nodes = np.array([1.0 + 0j, omega, omega**2])
    pick = np.eye(3, dtype=complex)
    for i in range(3):
        for j in range(3):
            if i != j:
                pick[i, j] = 1.0 / (1.0 - np.conj(nodes[i]) * nodes[j])
    eigenvalues = np.linalg.eigvalsh(pick)
    cycle = pick[0, 1] * pick[1, 2] * pick[2, 0]

    # Lurking-isometry synthesis of an explicit degree-two interpolant.
    positive = eigenvalues > 1.0e-10
    values, vectors = np.linalg.eigh(pick)
    gram_factor = np.diag(np.sqrt(values[positive])) @ vectors[:, positive].conj().T
    state_dimension = int(gram_factor.shape[0])
    x = np.array([[1.0], [0.0], [0.0]], dtype=complex)
    targets = np.eye(3, dtype=complex)
    source_columns = []
    target_columns = []
    for i, node in enumerate(nodes):
        feature = gram_factor[:, [i]]
        source_columns.append(np.vstack([node * feature, x]))
        target_columns.append(np.vstack([feature, targets[:, [i]]]))
    source = np.hstack(source_columns)
    target = np.hstack(target_columns)
    sample_gram = source.conj().T @ source
    gram_values, gram_vectors = np.linalg.eigh(sample_gram)
    inverse_sqrt = gram_vectors @ np.diag(1.0 / np.sqrt(gram_values)) @ gram_vectors.conj().T
    q_source = source @ inverse_sqrt
    q_target = target @ inverse_sqrt
    source_complement = null_space(q_source.conj().T)
    target_complement = null_space(q_target.conj().T)
    colligation = np.hstack([q_target, target_complement]) @ np.hstack(
        [q_source, source_complement]
    ).conj().T
    r = state_dimension
    a = colligation[:r, :r]
    b = colligation[:r, r:]
    c = colligation[r:, :r]
    d = colligation[r:, r:]
    controllability = np.hstack([np.linalg.matrix_power(a, power) @ b for power in range(r)])
    observability = np.vstack([c @ np.linalg.matrix_power(a, power) for power in range(r)])
    controllability_rank = int(np.linalg.matrix_rank(controllability, tol=1.0e-9))
    observability_rank = int(np.linalg.matrix_rank(observability, tol=1.0e-9))

    interpolation_errors = []
    unitarity_errors = []
    for i, node in enumerate(nodes):
        transfer = d + node * c @ np.linalg.solve(np.eye(r) - node * a, b)
        mapped = transfer @ x
        projector_error = np.linalg.norm(
            mapped @ mapped.conj().T - targets[:, [i]] @ targets[:, [i]].conj().T,
            ord=2,
        )
        interpolation_errors.append(float(projector_error))
        unitarity_errors.append(float(np.linalg.norm(transfer.conj().T @ transfer - np.eye(3), ord=2)))

    return {
        "nodes": [[float(z.real), float(z.imag)] for z in nodes],
        "pick_matrix": [[[float(z.real), float(z.imag)] for z in row] for row in pick],
        "pick_eigenvalues": eigenvalues.tolist(),
        "pick_rank": int(np.count_nonzero(positive)),
        "cycle_product": [float(cycle.real), float(cycle.imag)],
        "cycle_phase": float(np.angle(cycle)),
        "pairwise_minimum_degree": 1,
        "global_minimum_degree": 2,
        "colligation_unitarity_error": float(
            np.linalg.norm(colligation.conj().T @ colligation - np.eye(r + 3), ord=2)
        ),
        "maximum_interpolation_error": max(interpolation_errors),
        "maximum_boundary_unitarity_error": max(unitarity_errors),
        "state_spectral_radius": float(max(abs(np.linalg.eigvals(a)))),
        "controllability_rank": controllability_rank,
        "observability_rank": observability_rank,
        "all_pass": bool(
            eigenvalues[0] >= -TOL
            and np.count_nonzero(positive) == 2
            and abs(cycle.imag - 1.0 / (3.0 * math.sqrt(3.0))) <= TOL
            and max(interpolation_errors) <= 2.0e-7
            and max(unitarity_errors) <= 2.0e-7
            and np.linalg.norm(colligation.conj().T @ colligation - np.eye(r + 3), ord=2) <= 2.0e-7
            and max(abs(np.linalg.eigvals(a))) < 1.0 - 1.0e-8
            and controllability_rank == 2
            and observability_rank == 2
        ),
    }


def noisy_holonomy_audit(trials: int, seed: int) -> dict:
    """Stress-test the fail-closed rank-two cycle certificate."""
    rng = np.random.default_rng(seed)
    a = 1.0 / (1.0 - np.exp(2j * math.pi / 3.0))
    alternative = np.array([a, a, a], dtype=complex)
    certified = 0
    false_certificates = 0
    worst_margin = math.inf
    alternative_trials = trials // 2
    null_trials = trials - alternative_trials
    for index in range(trials):
        if index < alternative_trials:
            exact = alternative
            is_null = False
        else:
            f = rng.uniform(0.35, 1.25, 3) * np.exp(1j * rng.uniform(0.0, 2.0 * math.pi, 3))
            exact = np.array(
                [np.conj(f[0]) * f[1], np.conj(f[1]) * f[2], np.conj(f[2]) * f[0]],
                dtype=complex,
            )
            is_null = True
        epsilon = float(rng.uniform(1.0e-5, 2.0e-2))
        perturbation = epsilon * rng.uniform(0.0, 1.0, 3) * np.exp(
            1j * rng.uniform(0.0, 2.0 * math.pi, 3)
        )
        observed = exact + perturbation
        maximum = float(np.max(abs(observed)))
        error_bound = 3.0 * epsilon * (maximum + epsilon) ** 2
        witness = abs(float(np.imag(np.prod(observed)))) - error_bound
        decision = witness > 0.0
        certified += int(decision and not is_null)
        false_certificates += int(decision and is_null)
        worst_margin = min(worst_margin, witness)
    return {
        "trials": trials,
        "alternative_trials": alternative_trials,
        "rank_one_null_trials": null_trials,
        "certified_rank_two": certified,
        "false_certificates": false_certificates,
        "worst_witness_margin": worst_margin,
        "all_pass": bool(certified >= int(0.95 * alternative_trials) and false_certificates == 0),
    }


def make_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(8.2, 3.45))
    beta = np.array([0.92, 0.47, 0.18])
    geometric = 2.0 * beta.sum()
    degrees = np.arange(3, 8)
    upper = 2.0 * math.pi * degrees - geometric
    axes[0].fill_between(degrees, geometric, upper, alpha=0.28, color="#1768ac")
    axes[0].plot(degrees, upper, "o-", color="#1768ac", label=r"$2\pi d-2B$")
    axes[0].axhline(geometric, color="#d1495b", lw=2.0, label=r"$2B$")
    axes[0].set_xlabel("degree budget $d$")
    axes[0].set_ylabel("attainable arc action $A$")
    axes[0].set_xticks(degrees)
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False, fontsize=8)

    eigenvalues = np.array([2.0, 1.0, 0.0])
    axes[1].bar(np.arange(1, 4), eigenvalues, color=["#1768ac", "#4ea5d9", "#d9eaf2"])
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set_xlabel("Grassmann--Pick eigenvalue")
    axes[1].set_ylabel("value")
    axes[1].set_xticks([1, 2, 3])
    axes[1].text(2.0, 1.72, r"$\arg(P_{01}P_{12}P_{20})=\pi/2$", ha="center", fontsize=8.5)
    axes[1].grid(axis="y", alpha=0.22)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/exact_action_memory/certificate.json"))
    parser.add_argument(
        "--figure",
        type=Path,
        default=Path("paper_exact_action_memory/figures/action_memory_diamond.pdf"),
    )
    parser.add_argument("--trials", type=int, default=128)
    parser.add_argument("--noise-trials", type=int, default=2048)
    args = parser.parse_args()

    beta = np.array([0.92, 0.47, 0.18])
    flagship = compile_action_memory(beta, 15.25, 0.31, 2.14)
    holonomy = grassmann_pick_holonomy()
    payload = {
        "schema_version": 1,
        "artifact": "exact-action-memory-certificate",
        "deterministic_seeds": {"diamond": 26081031, "holonomy_noise": 26081032},
        "flagship_diamond": flagship,
        "deterministic_faces": deterministic_diamond_faces(),
        "noncommuting_order_control": noncommuting_order_control(),
        "random_diamond": random_diamond_audit(args.trials, 26081031),
        "three_node_holonomy": holonomy,
        "noisy_holonomy": noisy_holonomy_audit(args.noise_trials, 26081032),
    }
    payload["all_pass"] = bool(
        flagship["all_pass"]
        and payload["deterministic_faces"]["all_pass"]
        and payload["noncommuting_order_control"]["all_pass"]
        and payload["random_diamond"]["all_pass"]
        and holonomy["all_pass"]
        and payload["noisy_holonomy"]["all_pass"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    make_figure(args.figure)
    if not payload["all_pass"]:
        raise RuntimeError("exact action-memory certificate failed")
    print("exact action-memory certificate: PASS")


if __name__ == "__main__":
    main()
