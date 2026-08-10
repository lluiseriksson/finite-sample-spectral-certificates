#!/usr/bin/env python3
"""Independent verifier for the global fan-out memory artifact.

This file deliberately imports no functions from the producer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "global_fanout_memory" / "certificate.json"
TOL = 5.0e-8


def cauchy_off_diagonal(nodes: np.ndarray) -> np.ndarray:
    count = len(nodes)
    answer = np.zeros((count, count), dtype=complex)
    for row in range(count):
        for column in range(count):
            if row != column:
                answer[row, column] = 1.0 / (
                    1.0 - np.conj(nodes[row]) * nodes[column]
                )
    return answer


def verify_regular_spectra() -> None:
    for count in range(2, 18):
        nodes = np.exp(2j * np.pi * np.arange(count) / count)
        cauchy = cauchy_off_diagonal(nodes)
        expected = np.arange(count, dtype=float) - (count - 1.0) / 2.0
        if np.max(abs(np.linalg.eigvalsh(cauchy) - expected)) > 2.0e-11:
            raise RuntimeError("regular-polygon Cauchy spectrum changed")
        completion = cauchy + (count - 1.0) / 2.0 * np.eye(count)
        if np.linalg.matrix_rank(completion, tol=1.0e-8) != count - 1:
            raise RuntimeError("regular-polygon completion has the wrong rank")


def verify_arbitrary_nodes() -> None:
    angles = np.array([0.03, 0.61, 1.74, 2.20, 4.35, 5.91])
    nodes = np.exp(1j * angles)
    cauchy = cauchy_off_diagonal(nodes)
    smallest = np.linalg.eigvalsh(cauchy)[0]
    scalar = cauchy - smallest * np.eye(len(nodes))
    completion = np.kron(scalar, np.eye(2))
    block_nodes = np.kron(np.diag(nodes), np.eye(2))
    rhs = np.kron(np.ones((len(nodes), len(nodes))) - np.eye(len(nodes)), np.eye(2))
    residual = completion - block_nodes.conj().T @ completion @ block_nodes - rhs
    if np.linalg.eigvalsh(completion)[0] < -TOL:
        raise RuntimeError("arbitrary-node completion is not positive semidefinite")
    if np.linalg.matrix_rank(completion, tol=1.0e-8) != 2 * (len(nodes) - 1):
        raise RuntimeError("arbitrary-node sharp completion has the wrong rank")
    if np.linalg.norm(residual, ord=2) > TOL:
        raise RuntimeError("arbitrary-node Stein identity failed")


def deterministic_complement(q: np.ndarray) -> np.ndarray:
    rows, columns = q.shape
    vectors: list[np.ndarray] = []
    for index in range(rows):
        candidate = np.zeros(rows, dtype=complex)
        candidate[index] = 1.0
        candidate -= q @ (q.conj().T @ candidate)
        for vector in vectors:
            candidate -= vector * np.vdot(vector, candidate)
        length = np.linalg.norm(candidate)
        if length <= 1.0e-11:
            continue
        candidate /= length
        pivot = int(np.argmax(abs(candidate)))
        candidate *= np.exp(-1j * np.angle(candidate[pivot]))
        vectors.append(candidate)
        if len(vectors) == rows - columns:
            break
    if len(vectors) != rows - columns:
        raise RuntimeError("independent colligation complement failed")
    return np.column_stack(vectors) if vectors else np.zeros((rows, 0), dtype=complex)


def verify_one_colligation(count: int, channel_rank: int) -> None:
    nodes = np.exp(2j * np.pi * np.arange(count) / count)
    scalar = cauchy_off_diagonal(nodes) + (count - 1.0) / 2.0 * np.eye(count)
    pick = np.kron(scalar, np.eye(channel_rank))
    values, vectors = np.linalg.eigh(pick)
    positive = values > 1.0e-9
    feature_matrix = np.diag(np.sqrt(values[positive])) @ vectors[:, positive].conj().T
    state_dimension = feature_matrix.shape[0]
    expected = channel_rank * (count - 1)
    if state_dimension != expected:
        raise RuntimeError("independent Pick factor has the wrong state dimension")

    output_dimension = count * channel_rank
    input_frame = np.eye(output_dimension, channel_rank, dtype=complex)
    source_columns = []
    target_columns = []
    targets = []
    for index, node in enumerate(nodes):
        target_frame = np.zeros((output_dimension, channel_rank), dtype=complex)
        target_frame[index * channel_rank : (index + 1) * channel_rank, :] = np.eye(channel_rank)
        targets.append(target_frame)
        feature = feature_matrix[:, index * channel_rank : (index + 1) * channel_rank]
        source_columns.append(np.vstack([node * feature, input_frame]))
        target_columns.append(np.vstack([feature, target_frame]))
    source = np.hstack(source_columns)
    target = np.hstack(target_columns)
    if np.linalg.norm(source.conj().T @ source - target.conj().T @ target, ord=2) > TOL:
        raise RuntimeError("independent lurking-isometry Gram identity failed")

    sample_gram = source.conj().T @ source
    gram_values, gram_vectors = np.linalg.eigh(sample_gram)
    inverse_sqrt = gram_vectors @ np.diag(1.0 / np.sqrt(gram_values)) @ gram_vectors.conj().T
    q_source = source @ inverse_sqrt
    q_target = target @ inverse_sqrt
    source_full = np.hstack([q_source, deterministic_complement(q_source)])
    target_full = np.hstack([q_target, deterministic_complement(q_target)])
    colligation = target_full @ source_full.conj().T
    if np.linalg.norm(colligation.conj().T @ colligation - np.eye(colligation.shape[0]), ord=2) > TOL:
        raise RuntimeError("independently synthesized colligation is not unitary")

    r = state_dimension
    a = colligation[:r, :r]
    b = colligation[:r, r:]
    c = colligation[r:, :r]
    d = colligation[r:, r:]
    for index, node in enumerate(nodes):
        transfer = d + node * c @ np.linalg.solve(np.eye(r) - node * a, b)
        mapped = transfer @ input_frame
        error = np.linalg.norm(
            mapped @ mapped.conj().T - targets[index] @ targets[index].conj().T,
            ord=2,
        )
        if error > 3.0e-7:
            raise RuntimeError("independent colligation missed a fan-out target")
    controllability = np.hstack([np.linalg.matrix_power(a, power) @ b for power in range(r)])
    observability = np.vstack([c @ np.linalg.matrix_power(a, power) for power in range(r)])
    if np.linalg.matrix_rank(controllability, tol=1.0e-8) != r:
        raise RuntimeError("independent colligation is not controllable")
    if np.linalg.matrix_rank(observability, tol=1.0e-8) != r:
        raise RuntimeError("independent colligation is not observable")


def verify_colligations() -> None:
    for count, channel_rank in ((3, 1), (4, 1), (4, 2)):
        verify_one_colligation(count, channel_rank)


def verify_clustered_nodes() -> None:
    for gap in (1.0e-3, 1.0e-5, 1.0e-7, 1.0e-9):
        nodes = np.exp(1j * np.array([0.0, gap, 1.5, 3.0, 5.0]))
        cauchy = cauchy_off_diagonal(nodes)
        completion = cauchy - np.linalg.eigvalsh(cauchy)[0] * np.eye(5)
        eigenvalues = np.linalg.eigvalsh(completion)
        scale = max(1.0, float(np.max(abs(eigenvalues))))
        rank = int(np.count_nonzero(eigenvalues > 1.0e-9 * scale))
        if rank != 4 or eigenvalues[0] < -1.0e-9 * scale:
            raise RuntimeError("clustered-node completion lost its relative rank certificate")


def verify_inertia_lemma() -> None:
    rng = np.random.default_rng(26081044)
    for dimension in range(3, 16):
        for _ in range(20):
            rank = int(rng.integers(1, dimension))
            factor = rng.normal(size=(rank, dimension)) + 1j * rng.normal(
                size=(rank, dimension)
            )
            positive = factor.conj().T @ factor
            phases = np.exp(1j * rng.uniform(0.0, 2.0 * math.pi, dimension))
            unitary = np.diag(phases)
            displacement = positive - unitary.conj().T @ positive @ unitary
            negative = int(np.count_nonzero(np.linalg.eigvalsh(displacement) < -1.0e-8))
            observed_rank = int(np.linalg.matrix_rank(positive, tol=1.0e-8))
            if negative > observed_rank:
                raise RuntimeError("negative inertia exceeded the PSD rank")


def verify_robust_threshold() -> None:
    rng = np.random.default_rng(26081045)
    for count, channel_rank in ((3, 1), (5, 2), (8, 3)):
        dimension = count * channel_rank
        exact = np.kron(
            np.ones((count, count)) - np.eye(count), np.eye(channel_rank)
        )
        noise = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
            size=(dimension, dimension)
        )
        noise = 0.5 * (noise + noise.conj().T)
        # Keep a numerical margin below the declared eta=0.35.  A control
        # exactly on the strict threshold can cross it by one ulp on a
        # different BLAS/LAPACK implementation.
        noise *= 0.34 / np.linalg.norm(noise, ord=2)
        observed = exact + noise
        certified = int(np.count_nonzero(np.linalg.eigvalsh(observed) < -0.35))
        if certified != channel_rank * (count - 1):
            raise RuntimeError("robust inertia threshold lost a compulsory state")

        null_noise = noise.copy()
        false_certificate = int(
            np.count_nonzero(np.linalg.eigvalsh(null_noise) < -0.35)
        )
        if false_certificate:
            raise RuntimeError("fail-closed inertia rule produced a false positive")

    # Deterministic fixtures approach the strict boundary from below.  At
    # equality two target planes coincide and exactly k negative directions
    # disappear, providing a negative control for the open condition.
    for count, channel_rank in ((3, 1), (5, 1), (4, 2)):
        dimension = count * channel_rank
        frames = []
        for index in range(count):
            frame = np.zeros((dimension, channel_rank), dtype=complex)
            frame[index * channel_rank : (index + 1) * channel_rank, :] = np.eye(channel_rank)
            frames.append(frame)
        for overlap in (0.9, 0.99, 0.999, 0.999999, 1.0):
            perturbed = [frame.copy() for frame in frames]
            perturbed[1] = (
                overlap * frames[0]
                + math.sqrt(max(0.0, 1.0 - overlap**2)) * frames[1]
            )
            gram = np.block(
                [[left.conj().T @ right for right in perturbed] for left in perturbed]
            )
            deviation = float(np.linalg.norm(gram - np.eye(dimension), ord=2))
            rhs = np.kron(np.ones((count, count)), np.eye(channel_rank)) - gram
            negative = int(np.count_nonzero(np.linalg.eigvalsh(rhs) < -1.0e-8))
            expected = channel_rank * (count - 1 if overlap < 1.0 else count - 2)
            if abs(deviation - overlap) > 2.0e-10 or negative != expected:
                raise RuntimeError("robustness boundary fixture changed")


def verify_frozen_certificate() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or not payload.get("all_pass"):
        raise RuntimeError("frozen global fan-out certificate is absent or failing")
    arbitrary = payload["arbitrary_nodes"]
    if int(arbitrary["trials"]) < 350 or int(arbitrary["rank_failures"]) != 0:
        raise RuntimeError("arbitrary-node campaign is incomplete")
    if float(arbitrary["worst_stein_error"]) > 2.0e-10:
        raise RuntimeError("published Stein residual changed")
    polygon = payload["regular_polygon"]
    if len(polygon["cases"]) < 20 or float(polygon["worst_spectrum_error"]) > 2.0e-11:
        raise RuntimeError("regular-polygon campaign is incomplete")
    synthesis = payload["colligation_synthesis"]
    if not synthesis.get("all_pass") or len(synthesis.get("cases", [])) < 5:
        raise RuntimeError("colligation synthesis campaign is incomplete")
    for case in synthesis["cases"]:
        expected = int(case["channel_rank"]) * (int(case["nodes"]) - 1)
        if int(case["global_degree"]) != expected:
            raise RuntimeError("frozen global degree changed")
        if int(case["state_dimension"]) != expected:
            raise RuntimeError("colligation used a nonminimal state dimension")
        if int(case["controllability_rank"]) != expected:
            raise RuntimeError("frozen colligation is not controllable")
        if int(case["observability_rank"]) != expected:
            raise RuntimeError("frozen colligation is not observable")
        if float(case["maximum_interpolation_error"]) > 3.0e-7:
            raise RuntimeError("published interpolation tolerance changed")
        if float(case["colligation_unitarity_error"]) > 2.0e-8:
            raise RuntimeError("published colligation tolerance changed")
    noisy = payload["noisy_inertia"]
    if int(noisy["trials"]) < 1000:
        raise RuntimeError("noisy inertia campaign is too small")
    if int(noisy["false_certificates"]) or int(noisy["missed_full_certificates"]):
        raise RuntimeError("noisy fail-closed campaign changed")
    near = payload["near_orthogonal"]
    if int(near["trials"]) < 500 or int(near["failures"]):
        raise RuntimeError("near-orthogonal robustness campaign changed")
    if float(near["largest_gram_deviation"]) >= 1.0:
        raise RuntimeError("near-orthogonal fixtures left the robust theorem")
    clustered = payload["clustered_nodes"]
    if int(clustered["failures"]) or len(clustered["cases"]) < 4:
        raise RuntimeError("clustered-node conditioning campaign changed")
    boundary = payload["robust_boundary"]
    if int(boundary["failures"]) or len(boundary["cases"]) != 15:
        raise RuntimeError("robustness boundary campaign changed")
    if max(float(case["gram_deviation"]) for case in boundary["cases"] if float(case["overlap"]) < 1.0) < 0.99999:
        raise RuntimeError("robustness fixtures no longer approach the strict threshold")


def main() -> None:
    verify_regular_spectra()
    verify_arbitrary_nodes()
    verify_clustered_nodes()
    verify_inertia_lemma()
    verify_robust_threshold()
    verify_colligations()
    verify_frozen_certificate()
    print("global fan-out memory verifier: PASS")


if __name__ == "__main__":
    main()
