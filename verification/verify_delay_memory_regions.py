#!/usr/bin/env python3
"""Independent checks for the complete delay-memory region artifact."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "delay_memory_regions" / "certificate.json"
TOL = 5.0e-8


def close(left: float, right: float, tolerance: float = TOL) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def canonical_angles(unitary: np.ndarray, rank: int) -> np.ndarray:
    overlap = unitary[:rank, :rank]
    return np.sort(np.arccos(np.clip(np.linalg.svd(overlap, compute_uv=False), 0.0, 1.0)))[::-1]


def rebuild_compiler(record: dict) -> None:
    beta = np.asarray(record["beta"], dtype=float)
    q = np.asarray(record["target_delay_profile"], dtype=float)
    weights = np.asarray(record["weights"], dtype=float)
    vertices = np.asarray(record["signed_vertices"], dtype=float)
    if weights.size > beta.size + 1 or not close(float(weights.sum()), 1.0):
        raise RuntimeError("invalid Caratheodory decomposition")
    if not np.allclose(weights @ vertices, 2.0 * beta, rtol=TOL, atol=TOL):
        raise RuntimeError("signed-permutahedron barycenter mismatch")
    unitary = np.eye(2 * beta.size, dtype=complex)
    integrated = np.zeros(beta.size)
    generators: list[np.ndarray] = []
    for weight, vertex in zip(weights, vertices, strict=True):
        generator = np.zeros_like(unitary)
        for j, value in enumerate(vertex):
            magnitude = abs(value)
            sign = 1.0 if value >= 0.0 else -1.0
            block = 0.5 * magnitude * np.array(
                [[1.0, -1j * sign], [1j * sign, 1.0]], dtype=complex
            )
            generator[np.ix_([j, beta.size + j], [j, beta.size + j])] = block
        eig = np.linalg.eigvalsh(generator)
        if eig[0] < -TOL or not np.allclose(eig[::-1][: beta.size], q, rtol=TOL, atol=TOL):
            raise RuntimeError("compiled segment is not PSD with the requested spectrum")
        generators.append(generator)
        integrated += weight * eig[::-1][: beta.size]
        unitary = unitary @ expm(1j * weight * generator)
    for first in generators:
        for second in generators:
            if np.linalg.norm(first @ second - second @ first, ord=2) > TOL:
                raise RuntimeError("compiled generators do not commute")
    if not np.allclose(canonical_angles(unitary, beta.size), beta, rtol=TOL, atol=TOL):
        raise RuntimeError("compiled endpoint angles are incorrect")
    if not np.allclose(integrated, q, rtol=TOL, atol=TOL):
        raise RuntimeError("compiled delay profile is incorrect")


def scalar_blaschke(z: complex, zero: complex) -> complex:
    return (z - zero) / (1.0 - np.conj(zero) * z)


def strict_blindness_check() -> None:
    nodes = np.array([1.0 + 0.0j, -1.0 + 0.0j])
    families = [[0.0j], [1j / math.sqrt(3.0)] * 2]
    ranks = []
    delays_by_family = []
    for zeros in families:
        values = np.array([np.prod([scalar_blaschke(z, a) for a in zeros]) for z in nodes])
        delays = np.array(
            [sum((1.0 - abs(a) ** 2) / abs(z - a) ** 2 for a in zeros) for z in nodes]
        )
        gram = np.diag(delays.astype(complex))
        cross = (1.0 - values[0] * np.conj(values[1])) / 2.0
        gram[0, 1] = cross
        gram[1, 0] = np.conj(cross)
        eig = np.linalg.eigvalsh(gram)
        if eig[0] < -TOL:
            raise RuntimeError("strict example Pick matrix is not PSD")
        ranks.append(int(np.count_nonzero(eig > TOL)))
        delays_by_family.append(delays)
    if not np.allclose(delays_by_family[0], delays_by_family[1], rtol=TOL, atol=TOL):
        raise RuntimeError("strict example does not share local delays")
    if ranks != [1, 2]:
        raise RuntimeError("strict example does not separate hidden memory ranks")


def independent_random_scalar_pick(seed: int = 26081021, trials: int = 256) -> None:
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        degree = int(rng.integers(1, 8))
        zeros = [
            rng.uniform(0.02, 0.85) * np.exp(1j * rng.uniform(0.0, 2.0 * np.pi))
            for _factor in range(degree)
        ]
        count = int(rng.integers(2, 8))
        nodes = np.exp(1j * np.sort(rng.uniform(0.0, 2.0 * np.pi, count)))
        values = np.array([np.prod([scalar_blaschke(z, a) for a in zeros]) for z in nodes])
        delays = np.array(
            [sum((1.0 - abs(a) ** 2) / abs(z - a) ** 2 for a in zeros) for z in nodes]
        )
        gram = np.empty((count, count), dtype=complex)
        for i in range(count):
            for j in range(count):
                gram[i, j] = delays[i] if i == j else (
                    (1.0 - values[i] * np.conj(values[j]))
                    / (1.0 - nodes[i] * np.conj(nodes[j]))
                )
        eig = np.linalg.eigvalsh((gram + gram.conj().T) / 2.0)
        if eig[0] < -2.0e-7 or np.count_nonzero(eig > 2.0e-7) > degree:
            raise RuntimeError("independent scalar Pick probe failed")


def independent_quantization_compiler(seed: int = 26081022, trials: int = 64) -> None:
    rng = np.random.default_rng(seed)

    def cayley(point: complex) -> float:
        return float(np.real(1j * (1.0 + point) / (1.0 - point)))

    for _ in range(trials):
        theta0 = float(rng.uniform(0.0, 2.0 * np.pi))
        arc = float(rng.uniform(0.3, 2.0 * np.pi - 0.3))
        beta = float(rng.uniform(0.02, 0.48 * np.pi))
        z0 = complex(np.exp(1j * theta0))
        z1 = complex(np.exp(1j * (theta0 + arc)))
        target = complex(np.exp(2j * beta))
        gamma = cayley(target) - cayley(z1 / z0)

        def value_derivative(z: complex) -> tuple[complex, complex]:
            u = z / z0
            numerator = gamma + (2j - gamma) * u
            denominator = 2j + gamma - gamma * u
            value = numerator / denominator
            derivative = (
                ((2j - gamma) * denominator + gamma * numerator)
                / denominator**2
                / z0
            )
            return complex(value), complex(derivative)

        first, _ = value_derivative(z0)
        second, _ = value_derivative(z1)
        observed = math.acos(float(np.clip(abs((1.0 + second) / 2.0), 0.0, 1.0)))
        if abs(first - 1.0) > TOL or abs(second - target) > TOL or abs(observed - beta) > TOL:
            raise RuntimeError("independent Blaschke endpoint compiler failed")

        def delay(theta: float) -> float:
            point = complex(np.exp(1j * theta))
            value, derivative = value_derivative(point)
            return float(np.real(np.conj(value) * point * derivative))

        action, _ = quad(delay, theta0, theta0 + arc, epsabs=2.0e-10, epsrel=2.0e-10)
        if abs(action - 2.0 * beta) > 2.0e-8:
            raise RuntimeError("independent Blaschke action compiler failed")


def independent_projector_noise(seed: int = 26081023, trials: int = 128) -> None:
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        rank = int(rng.integers(1, 5))
        dimension = 2 * rank
        beta = np.sort(rng.uniform(0.0, 0.49 * np.pi, rank))[::-1]
        p0 = np.diag(np.concatenate([np.ones(rank), np.zeros(rank)])).astype(complex)
        frame = np.zeros((dimension, rank), dtype=complex)
        for j, angle in enumerate(beta):
            frame[j, j] = np.cos(angle)
            frame[rank + j, j] = np.sin(angle)
        p1 = frame @ frame.conj().T
        estimates = []
        eta = 0.0
        for projector in (p0, p1):
            h = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
                size=(dimension, dimension)
            )
            h = (h + h.conj().T) / 2.0
            h /= max(np.linalg.norm(h, ord=2), 1.0e-12)
            u = expm(1j * 1.0e-3 * h)
            estimate = u @ projector @ u.conj().T
            estimates.append(estimate)
            eta += float(np.linalg.norm(estimate - projector, ord=2))
        singular = np.linalg.svd((np.eye(dimension) - estimates[0]) @ estimates[1], compute_uv=False)
        certified = np.arcsin(np.clip(np.maximum(singular[:rank] - eta, 0.0), 0.0, 1.0))
        if np.max(np.cumsum(certified) - np.cumsum(beta)) > TOL:
            raise RuntimeError("independent noisy-projector Lorenz certificate failed")


def independent_matrix_pick_noise() -> None:
    dimension = 2
    zeros = [0.2 + 0.1j, -0.3 + 0.15j, 0.1 - 0.4j]
    vectors = [
        np.array([1.0, 0.0], dtype=complex),
        np.array([1.0, 1j], dtype=complex) / math.sqrt(2.0),
        np.array([1.0, 1.0], dtype=complex) / math.sqrt(2.0),
    ]
    projectors = [np.outer(vector, vector.conj()) for vector in vectors]
    nodes = np.exp(1j * np.array([0.2, 2.1, 4.4]))

    def factor_data(z: complex) -> tuple[np.ndarray, np.ndarray]:
        identity = np.eye(dimension, dtype=complex)
        factors = []
        derivatives = []
        for zero, projector in zip(zeros, projectors, strict=True):
            b = (z - zero) / (1.0 - np.conj(zero) * z)
            bp = (1.0 - abs(zero) ** 2) / (1.0 - np.conj(zero) * z) ** 2
            factors.append(identity + (b - 1.0) * projector)
            derivatives.append(bp * projector)
        value = identity.copy()
        for factor in factors:
            value = value @ factor
        derivative = np.zeros_like(identity)
        for index in range(len(factors)):
            term = identity.copy()
            for j, factor in enumerate(factors):
                term = term @ (derivatives[j] if j == index else factor)
            derivative += term
        return value, derivative

    def gram(values: list[np.ndarray], delays: list[np.ndarray]) -> np.ndarray:
        result = np.empty((nodes.size * dimension, nodes.size * dimension), dtype=complex)
        identity = np.eye(dimension, dtype=complex)
        for i in range(nodes.size):
            for j in range(nodes.size):
                rows = slice(i * dimension, (i + 1) * dimension)
                cols = slice(j * dimension, (j + 1) * dimension)
                result[rows, cols] = (
                    values[i] @ delays[i] @ values[i].conj().T
                    if i == j
                    else (identity - values[i] @ values[j].conj().T)
                    / (1.0 - nodes[i] * np.conj(nodes[j]))
                )
        return (result + result.conj().T) / 2.0

    values = []
    delays = []
    for node in nodes:
        value, derivative = factor_data(complex(node))
        delay = value.conj().T @ (node * derivative)
        values.append(value)
        delays.append((delay + delay.conj().T) / 2.0)
    exact = gram(values, delays)
    eig = np.linalg.eigvalsh(exact)
    if eig[0] < -2.0e-7 or np.count_nonzero(eig > 2.0e-7) > len(zeros):
        raise RuntimeError("independent matrix Pick factorization failed")

    noisy_values = []
    noisy_delays = []
    eps_u = 0.0
    eps_q = 0.0
    q_max = max(float(np.linalg.norm(delay, ord=2)) for delay in delays)
    perturbations = [
        np.array([[0.4, 0.2j], [-0.2j, -0.4]], dtype=complex),
        np.array([[0.1, 0.3], [0.3, -0.1]], dtype=complex),
        np.array([[-0.2, 0.1j], [-0.1j, 0.2]], dtype=complex),
    ]
    for value, delay, h in zip(values, delays, perturbations, strict=True):
        h /= np.linalg.norm(h, ord=2)
        noisy_value = value @ expm(1j * 1.5e-4 * h)
        noisy_delay = delay + 1.0e-4 * h
        noisy_values.append(noisy_value)
        noisy_delays.append(noisy_delay)
        eps_u = max(eps_u, float(np.linalg.norm(noisy_value - value, ord=2)))
        eps_q = max(eps_q, float(np.linalg.norm(noisy_delay - delay, ord=2)))
    observed = gram(noisy_values, noisy_delays)
    separation = min(
        abs(1.0 - nodes[i] * np.conj(nodes[j]))
        for i in range(nodes.size)
        for j in range(nodes.size)
        if i != j
    )
    eta = eps_q + 2.0 * q_max * eps_u + (nodes.size - 1) * 2.0 * eps_u / separation
    if np.linalg.norm(observed - exact, ord=2) > eta + TOL:
        raise RuntimeError("independent noisy matrix Pick bound failed")
    if np.count_nonzero(np.linalg.eigvalsh(observed) > eta) > len(zeros):
        raise RuntimeError("independent noisy matrix Pick rank overcount")


def mass_gap_check(record: dict) -> None:
    masses = np.asarray(record["masses"], dtype=float)
    delays = np.asarray(record["zero_frequency_delays"], dtype=float)
    recovered = 2.0 * np.arctanh(1.0 / delays)
    if not np.allclose(recovered, masses, rtol=TOL, atol=TOL):
        raise RuntimeError("mass/proper-delay inversion failed")
    su2 = record["su2_2d"]
    if not close(float(su2["gap"]), 3.0 * float(su2["area_parameter"]) / 4.0):
        raise RuntimeError("SU(2) gap arithmetic failed")
    if not close(float(su2["delay_cap"]), 1.0 / math.tanh(float(su2["gap"]) / 2.0)):
        raise RuntimeError("SU(2) delay cap arithmetic failed")


def main() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or not payload.get("all_pass"):
        raise RuntimeError("frozen certificate is not marked passing")
    rebuild_compiler(payload["flagship_complete_region"])
    strict_blindness_check()
    independent_random_scalar_pick()
    independent_quantization_compiler()
    independent_projector_noise()
    independent_matrix_pick_noise()
    mass_gap_check(payload["transfer_mass_gap_dictionary"])
    minimums = {
        "random_complete_region": 190,
        "random_rational_quantization": 250,
        "noisy_projector_certificate": 500,
        "random_pick_gram": 750,
        "robust_pick_rank": 500,
    }
    for key, minimum in minimums.items():
        if int(payload[key]["trials"]) < minimum or not payload[key]["all_pass"]:
            raise RuntimeError(f"campaign {key} is too small or failed")
    if payload["local_delay_blindness"]["maximum_local_delay_difference"] > TOL:
        raise RuntimeError("local-delay blindness arithmetic changed")
    if payload["robust_pick_rank"]["rank_overcounts"] != 0:
        raise RuntimeError("robust rank certificate overcounted memory")
    print("complete delay-memory region verifier: PASS")


if __name__ == "__main__":
    main()
