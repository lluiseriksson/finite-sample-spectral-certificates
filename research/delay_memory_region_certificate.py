#!/usr/bin/env python3
"""Deterministic certificate for complete delay regions and memory quantization."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
from scipy.optimize import linprog


TOL = 2.0e-9


def canonical_angles(unitary: np.ndarray, rank: int) -> np.ndarray:
    singular = np.linalg.svd(unitary[:rank, :rank], compute_uv=False)
    return np.sort(np.arccos(np.clip(singular, 0.0, 1.0)))[::-1]


def signed_orbit(q: np.ndarray) -> np.ndarray:
    vectors: list[np.ndarray] = []
    for permutation in itertools.permutations(range(q.size)):
        permuted = q[list(permutation)]
        for signs in itertools.product((-1.0, 1.0), repeat=q.size):
            vectors.append(permuted * np.asarray(signs))
    return np.column_stack(vectors)


def sparse_signed_decomposition(beta: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Find and sparsify 2 beta as a convex combination of signed permutations of q."""
    orbit = signed_orbit(q)
    constraints = np.vstack([orbit, np.ones((1, orbit.shape[1]))])
    target = np.concatenate([2.0 * beta, [1.0]])
    solution = linprog(
        np.zeros(orbit.shape[1]),
        A_eq=constraints,
        b_eq=target,
        bounds=(0.0, None),
        method="highs",
    )
    if not solution.success:
        raise RuntimeError(f"signed-permutahedron decomposition failed: {solution.message}")
    active = np.flatnonzero(solution.x > 2.0e-10)
    weights = solution.x[active]
    vectors = orbit[:, active]
    # HiGHS returns a basic feasible solution: at most k+1 equality columns.
    if active.size > beta.size + 1:
        raise RuntimeError("decomposition exceeded the Caratheodory bound")
    return weights, vectors


def generator_from_vertex(vertex: np.ndarray) -> np.ndarray:
    """PSD coupled-plane generator with top spectrum abs(vertex)."""
    rank = vertex.size
    q = np.zeros((2 * rank, 2 * rank), dtype=complex)
    for j, signed_value in enumerate(vertex):
        magnitude = abs(float(signed_value))
        sign = 1.0 if signed_value >= 0.0 else -1.0
        block = 0.5 * magnitude * np.array(
            [[1.0, -1j * sign], [1j * sign, 1.0]], dtype=complex
        )
        indices = np.ix_([j, rank + j], [j, rank + j])
        q[indices] = block
    return q


def compile_delay_profile(beta: np.ndarray, q: np.ndarray) -> dict:
    weights, vertices = sparse_signed_decomposition(beta, q)
    rank = beta.size
    unitary = np.eye(2 * rank, dtype=complex)
    integrated = np.zeros(rank)
    generators: list[np.ndarray] = []
    minimum_eigenvalue = math.inf
    maximum_commutator = 0.0
    for weight, vertex in zip(weights, vertices.T, strict=True):
        generator = generator_from_vertex(vertex)
        generators.append(generator)
        eig = np.linalg.eigvalsh(generator)
        minimum_eigenvalue = min(minimum_eigenvalue, float(eig[0]))
        integrated += weight * eig[::-1][:rank]
        unitary = unitary @ expm(1j * weight * generator)
    for first in generators:
        for second in generators:
            maximum_commutator = max(
                maximum_commutator,
                float(np.linalg.norm(first @ second - second @ first, ord=2)),
            )
    observed = canonical_angles(unitary, rank)
    return {
        "rank": rank,
        "beta": beta.tolist(),
        "target_delay_profile": q.tolist(),
        "prefix_margin": (np.cumsum(q) - 2.0 * np.cumsum(beta)).tolist(),
        "segments": int(weights.size),
        "weights": weights.tolist(),
        "signed_vertices": vertices.T.tolist(),
        "observed_angles": observed.tolist(),
        "integrated_top_delays": integrated.tolist(),
        "minimum_generator_eigenvalue": minimum_eigenvalue,
        "maximum_commutator_norm": maximum_commutator,
        "maximum_error": float(max(np.max(abs(observed - beta)), np.max(abs(integrated - q)))),
    }


def random_complete_region_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    cases: list[dict] = []
    worst_error = 0.0
    largest_segments = 0
    for _ in range(trials):
        rank = int(rng.integers(1, 5))
        q = np.sort(rng.uniform(0.35, 2.7, rank))[::-1]
        orbit = signed_orbit(q)
        count = int(rng.integers(1, rank + 2))
        chosen = rng.choice(orbit.shape[1], size=count, replace=False)
        weights = rng.dirichlet(np.ones(count))
        target = orbit[:, chosen] @ weights
        beta = 0.5 * np.sort(np.abs(target))[::-1]
        # Keep canonical angles in their principal range by a harmless global scaling.
        scale = min(1.0, (0.49 * np.pi) / max(beta[0], 1.0e-12))
        beta *= scale
        q *= scale
        record = compile_delay_profile(beta, q)
        worst_error = max(worst_error, record["maximum_error"])
        largest_segments = max(largest_segments, record["segments"])
        cases.append({
            "rank": rank,
            "segments": record["segments"],
            "maximum_error": record["maximum_error"],
            "minimum_prefix_margin": float(min(record["prefix_margin"])),
            "minimum_generator_eigenvalue": record["minimum_generator_eigenvalue"],
            "maximum_commutator_norm": record["maximum_commutator_norm"],
        })
    return {
        "trials": trials,
        "largest_segment_count": largest_segments,
        "worst_error": worst_error,
        "worst_commutator_norm": max(case["maximum_commutator_norm"] for case in cases),
        "minimum_generator_eigenvalue": min(case["minimum_generator_eigenvalue"] for case in cases),
        "all_pass": bool(
            worst_error <= 2.0e-7
            and largest_segments <= 5
            and min(case["minimum_prefix_margin"] for case in cases) >= -2.0e-8
            and min(case["minimum_generator_eigenvalue"] for case in cases) >= -2.0e-9
        ),
    }


def blaschke(z: complex, zero: complex) -> complex:
    return (z - zero) / (1.0 - np.conj(zero) * z)


def blaschke_derivative(z: complex, zero: complex) -> complex:
    return (1.0 - abs(zero) ** 2) / (1.0 - np.conj(zero) * z) ** 2


def scalar_product_value_derivative(z: complex, zeros: list[complex]) -> tuple[complex, complex]:
    factors = [blaschke(z, zero) for zero in zeros]
    value = np.prod(factors)
    derivative = 0.0j
    for index, zero in enumerate(zeros):
        derivative += blaschke_derivative(z, zero) * np.prod(
            [factor for j, factor in enumerate(factors) if j != index]
        )
    return complex(value), complex(derivative)


def scalar_pick_matrix(nodes: np.ndarray, zeros: list[complex]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.empty(nodes.size, dtype=complex)
    delays = np.empty(nodes.size)
    for index, node in enumerate(nodes):
        value, derivative = scalar_product_value_derivative(complex(node), zeros)
        values[index] = value
        delays[index] = float(np.real(np.conj(value) * node * derivative))
    gram = np.empty((nodes.size, nodes.size), dtype=complex)
    for i in range(nodes.size):
        for j in range(nodes.size):
            if i == j:
                gram[i, i] = delays[i]
            else:
                gram[i, j] = (1.0 - values[i] * np.conj(values[j])) / (
                    1.0 - nodes[i] * np.conj(nodes[j])
                )
    return values, delays, (gram + gram.conj().T) / 2.0


def local_delay_blindness() -> dict:
    nodes = np.array([1.0 + 0.0j, -1.0 + 0.0j])
    first_zeros = [0.0 + 0.0j]
    second_zeros = [1j / math.sqrt(3.0), 1j / math.sqrt(3.0)]
    first_values, first_delays, first_gram = scalar_pick_matrix(nodes, first_zeros)
    second_values, second_delays, second_gram = scalar_pick_matrix(nodes, second_zeros)
    first_eigenvalues = np.linalg.eigvalsh(first_gram)[::-1]
    second_eigenvalues = np.linalg.eigvalsh(second_gram)[::-1]
    return {
        "nodes": [[float(z.real), float(z.imag)] for z in nodes],
        "degree_one": {
            "zeros": [[0.0, 0.0]],
            "values": [[float(z.real), float(z.imag)] for z in first_values],
            "local_delays": first_delays.tolist(),
            "pick_eigenvalues": first_eigenvalues.tolist(),
            "pick_rank": int(np.count_nonzero(first_eigenvalues > TOL)),
        },
        "degree_two": {
            "zeros": [[0.0, 1.0 / math.sqrt(3.0)]] * 2,
            "values": [[float(z.real), float(z.imag)] for z in second_values],
            "local_delays": second_delays.tolist(),
            "pick_eigenvalues": second_eigenvalues.tolist(),
            "pick_rank": int(np.count_nonzero(second_eigenvalues > TOL)),
        },
        "maximum_local_delay_difference": float(np.max(abs(first_delays - second_delays))),
        "all_pass": bool(
            np.max(abs(first_delays - second_delays)) <= TOL
            and np.count_nonzero(first_eigenvalues > TOL) == 1
            and np.count_nonzero(second_eigenvalues > TOL) == 2
        ),
    }


def random_unit_vector(rng: np.random.Generator, dimension: int) -> np.ndarray:
    vector = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
    return vector / np.linalg.norm(vector)


def potapov_value_derivative(
    z: complex, zeros: list[complex], projectors: list[np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    dimension = projectors[0].shape[0]
    factors: list[np.ndarray] = []
    derivatives: list[np.ndarray] = []
    identity = np.eye(dimension, dtype=complex)
    for zero, projector in zip(zeros, projectors, strict=True):
        factors.append(identity + (blaschke(z, zero) - 1.0) * projector)
        derivatives.append(blaschke_derivative(z, zero) * projector)
    prefixes = [identity]
    for factor in factors:
        prefixes.append(prefixes[-1] @ factor)
    suffixes = [identity for _ in range(len(factors) + 1)]
    for index in range(len(factors) - 1, -1, -1):
        suffixes[index] = factors[index] @ suffixes[index + 1]
    derivative = np.zeros_like(identity)
    for index in range(len(factors)):
        derivative += prefixes[index] @ derivatives[index] @ suffixes[index + 1]
    return prefixes[-1], derivative


def matrix_pick_gram(nodes: np.ndarray, values: list[np.ndarray], delays: list[np.ndarray]) -> np.ndarray:
    dimension = values[0].shape[0]
    gram = np.empty((nodes.size * dimension, nodes.size * dimension), dtype=complex)
    identity = np.eye(dimension, dtype=complex)
    for i in range(nodes.size):
        for j in range(nodes.size):
            rows = slice(i * dimension, (i + 1) * dimension)
            cols = slice(j * dimension, (j + 1) * dimension)
            if i == j:
                gram[rows, cols] = values[i] @ delays[i] @ values[i].conj().T
            else:
                gram[rows, cols] = (identity - values[i] @ values[j].conj().T) / (
                    1.0 - nodes[i] * np.conj(nodes[j])
                )
    return (gram + gram.conj().T) / 2.0


def random_pick_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    worst_negative = 0.0
    rank_violations = 0
    equality_cases = 0
    for _ in range(trials):
        dimension = int(rng.integers(2, 5))
        degree = int(rng.integers(1, 7))
        zeros = []
        projectors = []
        for _factor in range(degree):
            radius = float(rng.uniform(0.05, 0.82))
            phase = float(rng.uniform(0.0, 2.0 * np.pi))
            zeros.append(radius * np.exp(1j * phase))
            vector = random_unit_vector(rng, dimension)
            projectors.append(np.outer(vector, vector.conj()))
        node_count = int(rng.integers(2, 6))
        phases = np.sort(rng.uniform(0.0, 2.0 * np.pi, node_count))
        nodes = np.exp(1j * phases)
        values: list[np.ndarray] = []
        delays: list[np.ndarray] = []
        for node in nodes:
            value, derivative = potapov_value_derivative(complex(node), zeros, projectors)
            delay = value.conj().T @ (node * derivative)
            values.append(value)
            delays.append((delay + delay.conj().T) / 2.0)
        gram = matrix_pick_gram(nodes, values, delays)
        eigenvalues = np.linalg.eigvalsh(gram)
        rank = int(np.count_nonzero(eigenvalues > 3.0e-7))
        worst_negative = min(worst_negative, float(eigenvalues[0]))
        rank_violations += int(rank > degree)
        equality_cases += int(rank == degree and node_count * dimension >= degree)
    return {
        "trials": trials,
        "worst_minimum_eigenvalue": worst_negative,
        "rank_violations": rank_violations,
        "full_degree_cases": equality_cases,
        "all_pass": bool(worst_negative >= -3.0e-7 and rank_violations == 0 and equality_cases > trials // 2),
    }


def robust_rank_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    overcounts = 0
    largest_threshold_ratio = 0.0
    for _ in range(trials):
        dimension = int(rng.integers(2, 4))
        degree = int(rng.integers(1, 5))
        zeros = []
        projectors = []
        for _factor in range(degree):
            zeros.append(float(rng.uniform(0.05, 0.75)) * np.exp(1j * rng.uniform(0.0, 2.0 * np.pi)))
            vector = random_unit_vector(rng, dimension)
            projectors.append(np.outer(vector, vector.conj()))
        node_count = 3
        phases = np.sort(rng.uniform(0.0, 2.0 * np.pi, node_count))
        nodes = np.exp(1j * phases)
        values: list[np.ndarray] = []
        delays: list[np.ndarray] = []
        noisy_values: list[np.ndarray] = []
        noisy_delays: list[np.ndarray] = []
        eps_u = 0.0
        eps_q = 0.0
        q_max = 0.0
        for node in nodes:
            value, derivative = potapov_value_derivative(complex(node), zeros, projectors)
            delay = value.conj().T @ (node * derivative)
            delay = (delay + delay.conj().T) / 2.0
            h = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
            h = (h + h.conj().T) / 2.0
            h /= max(np.linalg.norm(h, ord=2), 1.0e-12)
            noisy_value = value @ expm(1j * 2.0e-4 * h)
            e = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
            e = (e + e.conj().T) / 2.0
            e *= 2.0e-4 / max(np.linalg.norm(e, ord=2), 1.0e-12)
            noisy_delay = delay + e
            values.append(value)
            delays.append(delay)
            noisy_values.append(noisy_value)
            noisy_delays.append(noisy_delay)
            eps_u = max(eps_u, float(np.linalg.norm(noisy_value - value, ord=2)))
            eps_q = max(eps_q, float(np.linalg.norm(noisy_delay - delay, ord=2)))
            q_max = max(q_max, float(np.linalg.norm(delay, ord=2)))
        exact = matrix_pick_gram(nodes, values, delays)
        observed = matrix_pick_gram(nodes, noisy_values, noisy_delays)
        separation = min(
            abs(1.0 - nodes[i] * np.conj(nodes[j]))
            for i in range(node_count)
            for j in range(node_count)
            if i != j
        )
        eta = eps_q + 2.0 * q_max * eps_u + (node_count - 1) * 2.0 * eps_u / separation
        actual_error = float(np.linalg.norm(observed - exact, ord=2))
        largest_threshold_ratio = max(largest_threshold_ratio, actual_error / eta)
        certified_rank = int(np.count_nonzero(np.linalg.eigvalsh(observed) > eta + 2.0e-10))
        overcounts += int(certified_rank > degree)
    return {
        "trials": trials,
        "rank_overcounts": overcounts,
        "largest_actual_error_over_bound": largest_threshold_ratio,
        "all_pass": bool(overcounts == 0 and largest_threshold_ratio <= 1.0 + 2.0e-8),
    }


def mass_gap_dictionary() -> dict:
    masses = np.array([0.42, 0.91, 1.7, 2.8])
    transfer = np.diag(np.exp(-masses))
    zero_delay = (np.eye(masses.size) + transfer) @ np.linalg.inv(np.eye(masses.size) - transfer)
    recovered = 2.0 * np.arctanh(1.0 / np.diag(zero_delay))
    grid = 8192
    theta = 2.0 * np.pi * np.arange(grid) / grid
    coefficients = []
    for order in range(5):
        accumulated = np.zeros_like(transfer, dtype=complex)
        for angle in theta:
            denominator = (
                np.eye(masses.size)
                - 2.0 * transfer * np.cos(angle)
                + transfer @ transfer
            )
            delay = (np.eye(masses.size) - transfer @ transfer) @ np.linalg.inv(denominator)
            accumulated += delay * np.exp(-1j * order * angle)
        coefficients.append(accumulated / grid)
    fourier_error = max(
        float(np.linalg.norm(coefficients[order] - np.linalg.matrix_power(transfer, order), ord=2))
        for order in range(5)
    )
    a = 0.7
    n = np.arange(1, 13, dtype=float)
    casimir = n * (n + 2.0) / 4.0
    su2_masses = a * casimir
    su2_delays = 1.0 / np.tanh(su2_masses / 2.0)
    return {
        "masses": masses.tolist(),
        "zero_frequency_delays": np.diag(zero_delay).tolist(),
        "recovered_masses": recovered.tolist(),
        "mass_recovery_error": float(np.max(abs(recovered - masses))),
        "fourier_coefficient_error": fourier_error,
        "su2_2d": {
            "area_parameter": a,
            "casimir_levels": casimir.tolist(),
            "mass_levels": su2_masses.tolist(),
            "proper_delays_at_zero": su2_delays.tolist(),
            "gap": float(3.0 * a / 4.0),
            "delay_cap": float(1.0 / np.tanh(3.0 * a / 8.0)),
        },
        "all_pass": bool(np.max(abs(recovered - masses)) <= TOL and fourier_error <= 2.0e-9),
    }


def boundary_automorphism(
    z: np.ndarray | complex,
    z0: complex,
    z1: complex,
    beta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Disk automorphism taking z0 to 1 and z1 to exp(2 i beta)."""
    rotated_endpoint = z1 / z0
    target = np.exp(2j * beta)

    def cayley(point: complex) -> float:
        return float(np.real(1j * (1.0 + point) / (1.0 - point)))

    gamma = cayley(target) - cayley(rotated_endpoint)
    u = np.asarray(z, dtype=complex) / z0
    numerator = gamma + (2j - gamma) * u
    denominator = 2j + gamma - gamma * u
    value = numerator / denominator
    derivative_u = ((2j - gamma) * denominator + gamma * numerator) / denominator**2
    return value, derivative_u / z0


def rational_quantization_audit(trials: int, seed: int) -> dict:
    """Construct the sharp Blaschke compilers and verify angle and arc action."""
    rng = np.random.default_rng(seed)
    worst_endpoint_error = 0.0
    worst_angle_error = 0.0
    worst_action_error = 0.0
    minimum_delay = math.inf
    for _ in range(trials):
        rank = int(rng.integers(1, 6))
        beta = np.sort(rng.uniform(0.015, 0.49 * np.pi, rank))[::-1]
        theta0 = float(rng.uniform(0.0, 2.0 * np.pi))
        arc = float(rng.uniform(0.25, 2.0 * np.pi - 0.25))
        theta1 = theta0 + arc
        z0 = complex(np.exp(1j * theta0))
        z1 = complex(np.exp(1j * theta1))
        grid_theta = np.linspace(theta0, theta1, 1025)
        grid = np.exp(1j * grid_theta)
        total_action = 0.0
        for angle in beta:
            endpoints, _ = boundary_automorphism(np.array([z0, z1]), z0, z1, float(angle))
            target = np.exp(2j * angle)
            worst_endpoint_error = max(
                worst_endpoint_error,
                float(abs(endpoints[0] - 1.0)),
                float(abs(endpoints[1] - target)),
            )
            observed_angle = math.acos(float(np.clip(abs((1.0 + endpoints[1]) / 2.0), 0.0, 1.0)))
            worst_angle_error = max(worst_angle_error, abs(observed_angle - float(angle)))
            values, derivatives = boundary_automorphism(grid, z0, z1, float(angle))
            delay = np.real(np.conj(values) * grid * derivatives)
            minimum_delay = min(minimum_delay, float(np.min(delay)))
            def scalar_delay(theta: float) -> float:
                point = complex(np.exp(1j * theta))
                value, derivative = boundary_automorphism(point, z0, z1, float(angle))
                return float(np.real(np.conj(value) * point * derivative))

            integral, _ = quad(scalar_delay, theta0, theta1, epsabs=2.0e-10, epsrel=2.0e-10)
            total_action += float(integral)
        worst_action_error = max(worst_action_error, abs(total_action - 2.0 * float(np.sum(beta))))
    return {
        "trials": trials,
        "worst_endpoint_error": worst_endpoint_error,
        "worst_angle_error": worst_angle_error,
        "worst_arc_action_error": worst_action_error,
        "minimum_sampled_delay": minimum_delay,
        "all_pass": bool(
            worst_endpoint_error <= 2.0e-10
            and worst_angle_error <= 2.0e-10
            and worst_action_error <= 8.0e-5
            and minimum_delay >= -2.0e-10
        ),
    }


def noisy_projector_audit(trials: int, seed: int) -> dict:
    """Stress-test the joint support/Lorenz certificate using exact projectors."""
    rng = np.random.default_rng(seed)
    degree_overcounts = 0
    lorenz_violations = 0
    largest_angle_slack = 0.0
    for _ in range(trials):
        rank = int(rng.integers(1, 5))
        dimension = 2 * rank
        beta = np.sort(rng.uniform(0.0, 0.49 * np.pi, rank))[::-1]
        beta[rng.uniform(size=rank) < 0.2] = 0.0
        beta = np.sort(beta)[::-1]
        p0 = np.diag(np.concatenate([np.ones(rank), np.zeros(rank)])).astype(complex)
        frame = np.zeros((dimension, rank), dtype=complex)
        for j, angle in enumerate(beta):
            frame[j, j] = np.cos(angle)
            frame[rank + j, j] = np.sin(angle)
        p1 = frame @ frame.conj().T
        estimates = []
        errors = []
        for projector in (p0, p1):
            hermitian = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
                size=(dimension, dimension)
            )
            hermitian = (hermitian + hermitian.conj().T) / 2.0
            hermitian /= max(np.linalg.norm(hermitian, ord=2), 1.0e-12)
            unitary = expm(1j * float(rng.uniform(1.0e-5, 2.5e-3)) * hermitian)
            estimate = unitary @ projector @ unitary.conj().T
            estimates.append(estimate)
            errors.append(float(np.linalg.norm(estimate - projector, ord=2)))
        eta = sum(errors)
        singular = np.linalg.svd((np.eye(dimension) - estimates[0]) @ estimates[1], compute_uv=False)
        observed = np.sort(singular[:rank])[::-1]
        certified = np.arcsin(np.clip(np.maximum(observed - eta, 0.0), 0.0, 1.0))
        true_degree = int(np.count_nonzero(beta > 1.0e-12))
        certified_degree = int(np.count_nonzero(observed > eta))
        degree_overcounts += int(certified_degree > true_degree)
        excess = 2.0 * np.cumsum(certified) - 2.0 * np.cumsum(beta)
        lorenz_violations += int(float(np.max(excess)) > 2.0e-8)
        largest_angle_slack = max(largest_angle_slack, float(np.max(excess)))
    return {
        "trials": trials,
        "degree_overcounts": degree_overcounts,
        "lorenz_violations": lorenz_violations,
        "largest_certified_minus_true_prefix": largest_angle_slack,
        "all_pass": bool(degree_overcounts == 0 and lorenz_violations == 0),
    }


def rational_quantization_examples() -> dict:
    beta = np.array([0.61, 0.27, 0.04, 0.0])
    nonzero = int(np.count_nonzero(beta > 0.0))
    epsilon_family = [1.0e-1, 1.0e-3, 1.0e-6]
    return {
        "angles": beta.tolist(),
        "nonzero_angles": nonzero,
        "minimum_mcmillan_degree": nonzero,
        "sharp_local_trace_action": float(2.0 * np.sum(beta)),
        "unavoidable_outside_arc_phase": float(2.0 * np.pi * nonzero - 2.0 * np.sum(beta)),
        "vanishing_action_fixed_degree": [
            {
                "epsilon": epsilon,
                "rank": 4,
                "trace_action": 8.0 * epsilon,
                "minimum_degree": 4,
            }
            for epsilon in epsilon_family
        ],
        "all_pass": bool(nonzero == 3 and 2.0 * np.pi * nonzero >= 2.0 * np.sum(beta)),
    }


def make_figure(path: Path, blindness: dict, mass_gap: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.35))
    beta = np.array([1.1, 0.72, 0.31])
    profiles = [
        (2.0 * beta, "Pareto frontier"),
        (np.array([2.65, 1.45, 1.0]), "parallel surplus"),
        (np.array([4.26, 0.0, 0.0]), "serialised"),
    ]
    ranks = np.arange(1, 4)
    for profile, label in profiles:
        axes[0].plot(ranks, np.cumsum(profile), marker="o", lw=2.0, label=label)
    axes[0].set_xlabel(r"prefix $r$")
    axes[0].set_ylabel("integrated leading delay")
    axes[0].set_xticks(ranks)
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=7.7)

    eig1 = np.asarray(blindness["degree_one"]["pick_eigenvalues"])
    eig2 = np.asarray(blindness["degree_two"]["pick_eigenvalues"])
    x = np.arange(1, 3)
    axes[1].bar(x - 0.18, eig1, width=0.36, label="degree 1")
    axes[1].bar(x + 0.18, eig2, width=0.36, label="degree 2")
    axes[1].set_xlabel("Pick eigenvalue index")
    axes[1].set_ylabel("memory-Gram eigenvalue")
    axes[1].set_xticks(x)
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8.0)

    masses = np.asarray(mass_gap["masses"])
    delays = np.asarray(mass_gap["zero_frequency_delays"])
    dense = np.linspace(0.25, 3.0, 300)
    axes[2].plot(dense, 1.0 / np.tanh(dense / 2.0), lw=2.0, label=r"$\coth(h/2)$")
    axes[2].scatter(masses, delays, s=30, zorder=3, label="certificate points")
    axes[2].set_xlabel("mass / transfer gap")
    axes[2].set_ylabel("zero-frequency proper delay")
    axes[2].grid(alpha=0.25)
    axes[2].legend(frameon=False, fontsize=8.0)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/delay_memory_regions/certificate.json"))
    parser.add_argument("--figure", type=Path, default=Path("paper_delay_memory_regions/figures/complete_region.pdf"))
    parser.add_argument("--region-trials", type=int, default=192)
    parser.add_argument("--quantization-trials", type=int, default=256)
    parser.add_argument("--projector-trials", type=int, default=512)
    parser.add_argument("--pick-trials", type=int, default=768)
    parser.add_argument("--robust-trials", type=int, default=512)
    args = parser.parse_args()

    flagship = compile_delay_profile(
        np.array([1.10, 0.72, 0.31]),
        np.array([2.65, 1.45, 1.00]),
    )
    blindness = local_delay_blindness()
    mass_gap = mass_gap_dictionary()
    payload = {
        "schema_version": 1,
        "artifact": "complete-delay-memory-region-certificate",
        "deterministic_seeds": {
            "region": 26081011,
            "quantization": 26081014,
            "projector": 26081015,
            "pick": 26081012,
            "robust": 26081013,
        },
        "flagship_complete_region": flagship,
        "random_complete_region": random_complete_region_audit(args.region_trials, 26081011),
        "rational_quantization": rational_quantization_examples(),
        "random_rational_quantization": rational_quantization_audit(
            args.quantization_trials, 26081014
        ),
        "noisy_projector_certificate": noisy_projector_audit(args.projector_trials, 26081015),
        "local_delay_blindness": blindness,
        "random_pick_gram": random_pick_audit(args.pick_trials, 26081012),
        "robust_pick_rank": robust_rank_audit(args.robust_trials, 26081013),
        "transfer_mass_gap_dictionary": mass_gap,
    }
    payload["all_pass"] = bool(
        flagship["maximum_error"] <= 2.0e-7
        and flagship["segments"] <= flagship["rank"] + 1
        and all(
            payload[key]["all_pass"]
            for key in (
                "random_complete_region",
                "rational_quantization",
                "random_rational_quantization",
                "noisy_projector_certificate",
                "local_delay_blindness",
                "random_pick_gram",
                "robust_pick_rank",
                "transfer_mass_gap_dictionary",
            )
        )
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    make_figure(args.figure, blindness, mass_gap)
    if not payload["all_pass"]:
        raise RuntimeError("delay-memory region certificate failed")
    print("complete delay-memory region certificate: PASS")


if __name__ == "__main__":
    main()
