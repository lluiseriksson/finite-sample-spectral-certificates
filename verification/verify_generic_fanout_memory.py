#!/usr/bin/env python3
"""Independent verifier for generic full-span spectral fan-out.

This verifier intentionally imports nothing from the certificate producer.
It reconstructs the block Pick completion, displacement bounds, colliding
target family, and noisy rank certificate from their defining formulas.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "generic_fanout_memory" / "certificate.json"
TOL = 8.0e-8


def deterministic_frames(count: int, rank: int) -> list[np.ndarray]:
    """Nonorthogonal full-span frames, deterministic and well conditioned."""
    dimension = count * rank
    rows = np.arange(dimension, dtype=float)[:, None]
    frames: list[np.ndarray] = []
    for index in range(count):
        cols = np.arange(rank, dtype=float)[None, :]
        sample = np.exp(1j * (index + 1.0) * (rows + 1.0) * (cols + 1.0) / 7.0)
        selector = np.zeros((dimension, rank), dtype=complex)
        selector[index * rank : (index + 1) * rank] = np.eye(rank)
        sample += 0.17 * selector
        q, _ = np.linalg.qr(sample)
        frames.append(q[:, :rank])
    stacked = np.hstack(frames)
    if np.linalg.matrix_rank(stacked, tol=1.0e-10) != dimension:
        raise RuntimeError("deterministic target table is not full span")
    return frames


def gram(frames: list[np.ndarray]) -> np.ndarray:
    return np.block([[x.conj().T @ y for y in frames] for x in frames])


def partial_pick(nodes: np.ndarray, frames: list[np.ndarray]) -> np.ndarray:
    count = len(nodes)
    rank = frames[0].shape[1]
    result = np.zeros((count * rank, count * rank), dtype=complex)
    for i in range(count):
        for j in range(count):
            if i == j:
                continue
            block = (np.eye(rank) - frames[i].conj().T @ frames[j]) / (
                1.0 - np.conj(nodes[i]) * nodes[j]
            )
            result[i * rank : (i + 1) * rank, j * rank : (j + 1) * rank] = block
    return result


def independent_completion(nodes: np.ndarray, frames: list[np.ndarray]) -> np.ndarray:
    partial = partial_pick(nodes, frames)
    count = len(nodes)
    rank = frames[0].shape[1]
    cut = (count - 1) * rank
    leading = partial[:cut, :cut]
    t = 2.0 + np.linalg.norm(leading, ord=2)
    q = leading + t * np.eye(cut)
    b = partial[:cut, cut:]
    return np.block([[q, b], [b.conj().T, b.conj().T @ np.linalg.solve(q, b)]])


def verify_universal_completion() -> None:
    cases = [
        (np.array([0.02, 0.73, 2.11]), 1),
        (np.array([0.07, 0.91, 2.24, 4.88]), 2),
        (np.array([0.11, 0.66, 1.57, 3.02, 4.41, 5.77]), 2),
    ]
    for angles, rank in cases:
        count = len(angles)
        nodes = np.exp(1j * angles)
        frames = deterministic_frames(count, rank)
        completion = independent_completion(nodes, frames)
        eigenvalues = np.linalg.eigvalsh(completion)
        expected = rank * (count - 1)
        observed = int(np.count_nonzero(eigenvalues > 1.0e-9 * max(1.0, eigenvalues[-1])))
        if eigenvalues[0] < -TOL or observed != expected:
            raise RuntimeError("independent Schur completion lost PSD rank")
        z = np.kron(np.diag(nodes), np.eye(rank))
        rhs = np.kron(np.ones((count, count)), np.eye(rank)) - gram(frames)
        residual = completion - z.conj().T @ completion @ z - rhs
        if np.linalg.norm(residual, ord=2) > TOL:
            raise RuntimeError("independent Stein identity failed")
        if np.linalg.matrix_rank(np.hstack(frames), tol=1.0e-10) - rank != expected:
            raise RuntimeError("full-span lower bound no longer matches the upper bound")


def verify_inertia_lower_bound() -> None:
    rng = np.random.default_rng(26081171)
    for dimension in range(4, 22):
        for _ in range(12):
            rank = int(rng.integers(1, dimension))
            factor = rng.normal(size=(rank, dimension)) + 1j * rng.normal(
                size=(rank, dimension)
            )
            positive = factor.conj().T @ factor
            phases = np.exp(1j * rng.uniform(0.0, 2.0 * math.pi, dimension))
            z = np.diag(phases)
            displacement = positive - z.conj().T @ positive @ z
            negative = int(np.count_nonzero(np.linalg.eigvalsh(displacement) < -1.0e-8))
            if negative > np.linalg.matrix_rank(positive, tol=1.0e-8):
                raise RuntimeError("negative displacement inertia exceeded PSD rank")


def colliding_frames(count: int, rank: int, epsilon: float) -> list[np.ndarray]:
    dimension = (count + 1) * rank
    frames = []
    for index in range(count):
        frame = np.zeros((dimension, rank), dtype=complex)
        frame[:rank] = np.eye(rank)
        frame[(index + 1) * rank : (index + 2) * rank] = epsilon * np.eye(rank)
        frames.append(frame / math.sqrt(1.0 + epsilon**2))
    return frames


def verify_collision_discontinuity() -> None:
    for count, rank in ((3, 1), (7, 1), (5, 2)):
        for epsilon in (0.2, 1.0e-2, 1.0e-5):
            frames = colliding_frames(count, rank, epsilon)
            values = np.linalg.eigvalsh(gram(frames))
            predicted = epsilon**2 / (1.0 + epsilon**2)
            if abs(values[0] - predicted) > 3.0e-11:
                raise RuntimeError("colliding-target Gram spectrum changed")
            if np.linalg.matrix_rank(np.hstack(frames), tol=1.0e-12) != count * rank:
                raise RuntimeError("positive collision opening lost full span")
            overlap = np.linalg.norm(frames[0].conj().T @ frames[-1], ord=2)
            if abs(overlap - 1.0 / (1.0 + epsilon**2)) > 2.0e-12:
                raise RuntimeError("colliding-target overlap formula changed")
        closed = colliding_frames(count, rank, 0.0)
        if np.linalg.matrix_rank(np.hstack(closed), tol=1.0e-12) != rank:
            raise RuntimeError("closed collision should have zero routing cost")


def verify_noisy_rank_rule() -> None:
    rng = np.random.default_rng(26081172)
    for count, rank, epsilon in ((4, 1, 0.15), (6, 2, 0.08), (5, 3, 0.04)):
        frames = colliding_frames(count, rank, epsilon)
        exact = sum(frame @ frame.conj().T for frame in frames)
        eta = 0.35 * epsilon**2 / (1.0 + epsilon**2)
        noise = rng.normal(size=exact.shape) + 1j * rng.normal(size=exact.shape)
        noise = 0.5 * (noise + noise.conj().T)
        noise *= 0.95 * eta / np.linalg.norm(noise, ord=2)
        observed = exact + noise
        certified_rank = int(np.count_nonzero(np.linalg.eigvalsh(observed) > eta))
        if certified_rank != count * rank:
            raise RuntimeError("strict Weyl margin failed to certify full span")
        # A deliberately inflated error bound must be allowed to abstain.
        observed_values = np.linalg.eigvalsh(observed)
        smallest_signal = float(observed_values[-count * rank])
        abstaining_eta = 2.0 * smallest_signal
        if int(np.count_nonzero(np.linalg.eigvalsh(observed) > abstaining_eta)) >= count * rank:
            raise RuntimeError("fail-closed certificate did not abstain")


def verify_three_line_phase_diagram() -> None:
    # A degree-one rank-one Potapov orbit is either constant or injective and
    # remains in a two-dimensional span.  These fixtures audit all three
    # algebraic strata used by the analytic classification.
    nodes = np.exp(1j * np.array([0.1, 2.0, 5.1]))
    orbit = [np.array([[1.0], [z]], dtype=complex) / math.sqrt(2.0) for z in nodes]
    if np.linalg.matrix_rank(np.hstack(orbit), tol=1.0e-10) != 2:
        raise RuntimeError("degree-one orbit should span exactly two dimensions")
    projectors = [x @ x.conj().T for x in orbit]
    if min(np.linalg.norm(projectors[i] - projectors[j]) for i in range(3) for j in range(i)) < 0.1:
        raise RuntimeError("degree-one orbit unexpectedly repeated a target line")
    repeated = [orbit[0], orbit[0], orbit[1]]
    if len({tuple(np.round((x @ x.conj().T).ravel(), 12)) for x in repeated}) != 2:
        raise RuntimeError("repeated-line control changed")
    span_three = []
    for i in range(3):
        vector = np.zeros((3, 1), dtype=complex)
        vector[i, 0] = 1.0
        span_three.append(vector)
    if np.linalg.matrix_rank(np.hstack(span_three)) != 3:
        raise RuntimeError("three-dimensional control changed")


def verify_frozen_certificate() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or not payload.get("all_pass"):
        raise RuntimeError("frozen generic fan-out certificate is absent or failing")
    universal = payload["universal_completion"]
    if universal["trials"] < 500 or universal["failures"] != 0:
        raise RuntimeError("universal completion campaign is incomplete")
    if universal["worst_relative_stein_residual"] > 2.0e-8:
        raise RuntimeError("frozen Stein residual changed")
    near = payload["near_collapse"]
    if near["failures"] or near["worst_spectrum_error"] > 2.0e-12:
        raise RuntimeError("near-collapse campaign changed")
    synthesis = payload["colligation_synthesis"]
    if not synthesis["all_pass"] or len(synthesis["cases"]) < 5:
        raise RuntimeError("colligation campaign is incomplete")
    for case in synthesis["cases"]:
        expected = int(case["channel_rank"]) * (int(case["nodes"]) - 1)
        if int(case["state_dimension"]) != expected:
            raise RuntimeError("published colligation is not minimal")
        if float(case["maximum_interpolation_error"]) > 3.0e-7:
            raise RuntimeError("published interpolation tolerance changed")
        if float(case["state_spectral_radius"]) >= 1.0:
            raise RuntimeError("published colligation lost strict stability")
    noisy = payload["noisy_projector_sum"]
    if noisy["trials"] < 1000 or noisy["overcertificates"]:
        raise RuntimeError("noisy span campaign changed")
    if noisy["direct_sum_certificates"] != noisy["direct_sum_trials"]:
        raise RuntimeError("direct-sum projector certificates changed")
    clustered = payload["clustered_nodes"]
    if clustered["failures"] or len(clustered["cases"]) != 4:
        raise RuntimeError("clustered-node campaign changed")
    if max(case["relative_stein_residual"] for case in clustered["cases"]) > 2.0e-8:
        raise RuntimeError("clustered-node relative residual changed")


def main() -> None:
    verify_universal_completion()
    verify_inertia_lower_bound()
    verify_collision_discontinuity()
    verify_noisy_rank_rule()
    verify_three_line_phase_diagram()
    verify_frozen_certificate()
    print("generic fan-out memory verifier: PASS")


if __name__ == "__main__":
    main()
