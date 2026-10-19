#!/usr/bin/env python3
"""Deterministic certificates for global spectral fan-out memory laws.

The analytic results live in ``paper_global_fanout_memory/main.tex``.  This
producer builds the sharp Pick completions, synthesizes representative
minimal unitary colligations, and stress-tests the fail-closed inertia bound.
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


TOL = 2.0e-8


def encode_complex_matrix(matrix: np.ndarray) -> list[list[list[float]]]:
    return [
        [[float(value.real), float(value.imag)] for value in row]
        for row in np.asarray(matrix, dtype=complex)
    ]


def regular_nodes(count: int) -> np.ndarray:
    return np.exp(2j * math.pi * np.arange(count) / count)


def phase_cauchy_matrix(nodes: np.ndarray) -> np.ndarray:
    nodes = np.asarray(nodes, dtype=complex)
    count = len(nodes)
    matrix = np.zeros((count, count), dtype=complex)
    for i in range(count):
        for j in range(count):
            if i != j:
                matrix[i, j] = 1.0 / (1.0 - np.conj(nodes[i]) * nodes[j])
    return matrix


def sharp_completion(nodes: np.ndarray, channel_rank: int) -> tuple[np.ndarray, float]:
    cauchy = phase_cauchy_matrix(nodes)
    smallest = float(np.linalg.eigvalsh(cauchy)[0])
    scalar = cauchy - smallest * np.eye(len(nodes))
    return np.kron(scalar, np.eye(channel_rank)), smallest


def orthogonal_frames(count: int, channel_rank: int) -> tuple[np.ndarray, list[np.ndarray]]:
    dimension = count * channel_rank
    input_frame = np.eye(dimension, channel_rank, dtype=complex)
    targets = []
    for index in range(count):
        frame = np.zeros((dimension, channel_rank), dtype=complex)
        frame[index * channel_rank : (index + 1) * channel_rank, :] = np.eye(channel_rank)
        targets.append(frame)
    return input_frame, targets


def block_gram(frames: list[np.ndarray]) -> np.ndarray:
    count = len(frames)
    channel_rank = frames[0].shape[1]
    gram = np.zeros((count * channel_rank, count * channel_rank), dtype=complex)
    for i in range(count):
        for j in range(count):
            gram[
                i * channel_rank : (i + 1) * channel_rank,
                j * channel_rank : (j + 1) * channel_rank,
            ] = frames[i].conj().T @ frames[j]
    return gram


def stein_rhs(frames: list[np.ndarray]) -> np.ndarray:
    count = len(frames)
    channel_rank = frames[0].shape[1]
    return np.kron(np.ones((count, count)), np.eye(channel_rank)) - block_gram(frames)


def inertia(matrix: np.ndarray, tolerance: float = 1.0e-9) -> tuple[int, int, int]:
    eigenvalues = np.linalg.eigvalsh(np.asarray(matrix, dtype=complex))
    negative = int(np.count_nonzero(eigenvalues < -tolerance))
    zero = int(np.count_nonzero(abs(eigenvalues) <= tolerance))
    positive = int(np.count_nonzero(eigenvalues > tolerance))
    return negative, zero, positive


def canonical_complement(isometry: np.ndarray, tolerance: float = 1.0e-11) -> np.ndarray:
    """Complete an isometry by deterministic projected coordinate vectors."""
    rows, columns = isometry.shape
    needed = rows - columns
    accepted: list[np.ndarray] = []
    for index in range(rows):
        vector = np.zeros(rows, dtype=complex)
        vector[index] = 1.0
        vector -= isometry @ (isometry.conj().T @ vector)
        for previous in accepted:
            vector -= previous * np.vdot(previous, vector)
        length = float(np.linalg.norm(vector))
        if length <= tolerance:
            continue
        vector /= length
        pivot = int(np.argmax(abs(vector)))
        vector *= np.exp(-1j * np.angle(vector[pivot]))
        accepted.append(vector)
        if len(accepted) == needed:
            break
    if len(accepted) != needed:
        raise RuntimeError("failed to construct a canonical orthogonal complement")
    return np.column_stack(accepted) if accepted else np.zeros((rows, 0), dtype=complex)


def unitary_colligation(
    nodes: np.ndarray,
    pick: np.ndarray,
    input_frame: np.ndarray,
    targets: list[np.ndarray],
) -> dict:
    """Build the lurking-isometry colligation associated with a singular Pick matrix."""
    nodes = np.asarray(nodes, dtype=complex)
    count = len(nodes)
    channel_rank = input_frame.shape[1]
    output_dimension = input_frame.shape[0]
    values, vectors = np.linalg.eigh(pick)
    positive = values > 1.0e-9
    gram_factor = np.diag(np.sqrt(values[positive])) @ vectors[:, positive].conj().T
    state_dimension = int(gram_factor.shape[0])

    source_columns: list[np.ndarray] = []
    target_columns: list[np.ndarray] = []
    for i, node in enumerate(nodes):
        feature = gram_factor[:, i * channel_rank : (i + 1) * channel_rank]
        source_columns.append(np.vstack([node * feature, input_frame]))
        target_columns.append(np.vstack([feature, targets[i]]))
    source = np.hstack(source_columns)
    target = np.hstack(target_columns)
    gram_error = float(np.linalg.norm(source.conj().T @ source - target.conj().T @ target, ord=2))

    sample_gram = source.conj().T @ source
    gram_values, gram_vectors = np.linalg.eigh(sample_gram)
    if gram_values[0] <= 1.0e-10:
        raise RuntimeError("sample vectors are not independent")
    inverse_sqrt = gram_vectors @ np.diag(1.0 / np.sqrt(gram_values)) @ gram_vectors.conj().T
    q_source = source @ inverse_sqrt
    q_target = target @ inverse_sqrt
    source_complement = canonical_complement(q_source)
    target_complement = canonical_complement(q_target)
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

    interpolation_errors = []
    boundary_unitarity_errors = []
    for i, node in enumerate(nodes):
        transfer = d + node * c @ np.linalg.solve(np.eye(r) - node * a, b)
        mapped = transfer @ input_frame
        mapped_projector = mapped @ mapped.conj().T
        target_projector = targets[i] @ targets[i].conj().T
        interpolation_errors.append(float(np.linalg.norm(mapped_projector - target_projector, ord=2)))
        boundary_unitarity_errors.append(
            float(np.linalg.norm(transfer.conj().T @ transfer - np.eye(output_dimension), ord=2))
        )

    grid_unitarity_errors = []
    for theta in np.linspace(0.013, 2.0 * math.pi - 0.017, 97):
        point = np.exp(1j * theta)
        transfer = d + point * c @ np.linalg.solve(np.eye(r) - point * a, b)
        grid_unitarity_errors.append(
            float(np.linalg.norm(transfer.conj().T @ transfer - np.eye(output_dimension), ord=2))
        )

    return {
        "state_dimension": state_dimension,
        "sample_gram_error": gram_error,
        "colligation_unitarity_error": float(
            np.linalg.norm(colligation.conj().T @ colligation - np.eye(r + output_dimension), ord=2)
        ),
        "maximum_interpolation_error": max(interpolation_errors),
        "maximum_boundary_unitarity_error": max(boundary_unitarity_errors + grid_unitarity_errors),
        "state_spectral_radius": float(max(abs(np.linalg.eigvals(a)))) if r else 0.0,
        "controllability_rank": int(np.linalg.matrix_rank(controllability, tol=1.0e-8)),
        "observability_rank": int(np.linalg.matrix_rank(observability, tol=1.0e-8)),
    }


def regular_polygon_audit(max_nodes: int) -> dict:
    cases = []
    worst_spectrum_error = 0.0
    worst_stein_error = 0.0
    for count in range(2, max_nodes + 1):
        nodes = regular_nodes(count)
        cauchy = phase_cauchy_matrix(nodes)
        expected_cauchy = np.arange(count, dtype=float) - (count - 1.0) / 2.0
        observed_cauchy = np.linalg.eigvalsh(cauchy)
        completion, _ = sharp_completion(nodes, 1)
        z = np.diag(nodes)
        rhs = np.ones((count, count)) - np.eye(count)
        spectrum_error = float(np.max(abs(observed_cauchy - expected_cauchy)))
        stein_error = float(np.linalg.norm(completion - z.conj().T @ completion @ z - rhs, ord=2))
        worst_spectrum_error = max(worst_spectrum_error, spectrum_error)
        worst_stein_error = max(worst_stein_error, stein_error)
        cases.append(
            {
                "nodes": count,
                "completion_rank": int(np.linalg.matrix_rank(completion, tol=1.0e-8)),
                "minimum_degree": count - 1,
                "completion_eigenvalues": np.linalg.eigvalsh(completion).tolist(),
            }
        )
    return {
        "cases": cases,
        "worst_spectrum_error": worst_spectrum_error,
        "worst_stein_error": worst_stein_error,
        "all_pass": bool(worst_spectrum_error < 2.0e-11 and worst_stein_error < 2.0e-11),
    }


def arbitrary_node_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    worst_psd_violation = 0.0
    worst_stein_error = 0.0
    rank_failures = 0
    minimum_node_separation = math.inf
    samples = []
    for trial in range(trials):
        count = int(rng.integers(3, 13))
        while True:
            angles = np.sort(rng.uniform(0.0, 2.0 * math.pi, count))
            gaps = np.diff(np.r_[angles, angles[0] + 2.0 * math.pi])
            if float(np.min(gaps)) > 0.025:
                break
        nodes = np.exp(1j * angles)
        completion, smallest = sharp_completion(nodes, 1)
        eigenvalues = np.linalg.eigvalsh(completion)
        z = np.diag(nodes)
        rhs = np.ones((count, count)) - np.eye(count)
        stein_error = float(np.linalg.norm(completion - z.conj().T @ completion @ z - rhs, ord=2))
        rank_tolerance = 1.0e-9 * max(1.0, float(np.max(abs(eigenvalues))))
        rank = int(np.count_nonzero(eigenvalues > rank_tolerance))
        worst_psd_violation = max(worst_psd_violation, float(max(0.0, -eigenvalues[0])))
        worst_stein_error = max(worst_stein_error, stein_error)
        rank_failures += int(rank != count - 1)
        minimum_node_separation = min(minimum_node_separation, float(np.min(gaps)))
        if trial < 8:
            samples.append(
                {
                    "nodes": count,
                    "lambda_min_cauchy": smallest,
                    "completion_rank": rank,
                    "minimum_completion_eigenvalue": float(eigenvalues[0]),
                    "second_completion_eigenvalue": float(eigenvalues[1]),
                }
            )
    return {
        "trials": trials,
        "rank_failures": rank_failures,
        "worst_psd_violation": worst_psd_violation,
        "worst_stein_error": worst_stein_error,
        "minimum_node_separation": minimum_node_separation,
        "sample": samples,
        "all_pass": bool(rank_failures == 0 and worst_psd_violation < TOL and worst_stein_error < 2.0e-10),
    }


def clustered_node_audit() -> dict:
    cases = []
    failures = 0
    for gap in (1.0e-3, 1.0e-5, 1.0e-7, 1.0e-9):
        angles = np.array([0.0, gap, 1.5, 3.0, 5.0])
        nodes = np.exp(1j * angles)
        completion, _ = sharp_completion(nodes, 1)
        eigenvalues = np.linalg.eigvalsh(completion)
        scale = max(1.0, float(np.max(abs(eigenvalues))))
        rank_tolerance = 1.0e-9 * scale
        rank = int(np.count_nonzero(eigenvalues > rank_tolerance))
        z = np.diag(nodes)
        rhs = np.ones((5, 5)) - np.eye(5)
        residual = float(np.linalg.norm(completion - z.conj().T @ completion @ z - rhs, ord=2))
        relative_residual = residual / scale
        passed = bool(rank == 4 and eigenvalues[0] >= -rank_tolerance and relative_residual < 2.0e-14)
        failures += int(not passed)
        cases.append(
            {
                "minimum_gap": gap,
                "completion_scale": scale,
                "rank_tolerance": rank_tolerance,
                "completion_rank": rank,
                "relative_stein_residual": relative_residual,
                "all_pass": passed,
            }
        )
    return {"cases": cases, "failures": failures, "all_pass": failures == 0}


def colligation_campaign() -> dict:
    specifications = [(3, 1), (4, 1), (5, 1), (4, 2), (5, 2)]
    cases = []
    all_pass = True
    for count, channel_rank in specifications:
        nodes = regular_nodes(count)
        pick, _ = sharp_completion(nodes, channel_rank)
        input_frame, targets = orthogonal_frames(count, channel_rank)
        synthesis = unitary_colligation(nodes, pick, input_frame, targets)
        expected = channel_rank * (count - 1)
        passed = bool(
            synthesis["state_dimension"] == expected
            and synthesis["controllability_rank"] == expected
            and synthesis["observability_rank"] == expected
            and synthesis["sample_gram_error"] < 2.0e-8
            and synthesis["colligation_unitarity_error"] < 2.0e-8
            and synthesis["maximum_interpolation_error"] < 3.0e-7
            and synthesis["maximum_boundary_unitarity_error"] < 3.0e-7
            and synthesis["state_spectral_radius"] < 1.0 - 1.0e-8
        )
        synthesis.update(
            {
                "nodes": count,
                "channel_rank": channel_rank,
                "pairwise_degree": channel_rank,
                "global_degree": expected,
                "local_action_degree_bound": int(math.ceil(count * channel_rank / 2.0)),
                "all_pass": passed,
            }
        )
        cases.append(synthesis)
        all_pass = all_pass and passed
    return {"cases": cases, "all_pass": all_pass}


def noisy_inertia_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    false_certificates = 0
    missed_full_certificates = 0
    smallest_margin = math.inf
    for trial in range(trials):
        count = int(rng.integers(3, 13))
        channel_rank = int(rng.integers(1, 4))
        dimension = count * channel_rank
        is_null = trial >= trials // 2
        exact = np.zeros((dimension, dimension), dtype=complex) if is_null else (
            np.kron(np.ones((count, count)) - np.eye(count), np.eye(channel_rank))
        )
        error_bound = float(rng.uniform(0.01, 0.45))
        noise = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))
        noise = 0.5 * (noise + noise.conj().T)
        noise *= error_bound * rng.uniform(0.1, 0.98) / np.linalg.norm(noise, ord=2)
        observed = exact + noise
        eigenvalues = np.linalg.eigvalsh(observed)
        certificate = int(np.count_nonzero(eigenvalues < -error_bound))
        if is_null:
            false_certificates += int(certificate > 0)
        else:
            expected = channel_rank * (count - 1)
            missed_full_certificates += int(certificate != expected)
            negative = np.sort(eigenvalues)[:expected]
            smallest_margin = min(smallest_margin, float(np.min(-error_bound - negative)))
    return {
        "trials": trials,
        "alternative_trials": trials // 2,
        "null_trials": trials - trials // 2,
        "false_certificates": false_certificates,
        "missed_full_certificates": missed_full_certificates,
        "smallest_certification_margin": smallest_margin,
        "all_pass": bool(false_certificates == 0 and missed_full_certificates == 0),
    }


def near_orthogonal_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    failures = 0
    largest_gram_deviation = 0.0
    largest_coherence_bound = 0.0
    for _ in range(trials):
        count = int(rng.integers(3, 9))
        channel_rank = int(rng.integers(1, 4))
        dimension = count * channel_rank
        _, base = orthogonal_frames(count, channel_rank)
        scale = float(rng.uniform(1.0e-4, 0.035 / count))
        frames = []
        for frame in base:
            perturbation = scale * (
                rng.normal(size=frame.shape) + 1j * rng.normal(size=frame.shape)
            )
            q, _ = np.linalg.qr(frame + perturbation)
            frames.append(q[:, :channel_rank])
        gram = block_gram(frames)
        deviation = float(np.linalg.norm(gram - np.eye(dimension), ord=2))
        coherence = max(
            float(np.linalg.norm(frames[i].conj().T @ frames[j], ord=2))
            for i in range(count)
            for j in range(i + 1, count)
        )
        negative, _, _ = inertia(stein_rhs(frames), 1.0e-8)
        expected = channel_rank * (count - 1)
        if deviation >= 1.0 or negative != expected:
            failures += 1
        largest_gram_deviation = max(largest_gram_deviation, deviation)
        largest_coherence_bound = max(largest_coherence_bound, (count - 1) * coherence)
    return {
        "trials": trials,
        "failures": failures,
        "largest_gram_deviation": largest_gram_deviation,
        "largest_gershgorin_bound": largest_coherence_bound,
        "all_pass": bool(failures == 0 and largest_gram_deviation < 1.0),
    }


def robust_boundary_audit() -> dict:
    """Probe the strict norm threshold and its equality failure deterministically."""
    cases = []
    failures = 0
    for count, channel_rank in ((3, 1), (5, 1), (4, 2)):
        dimension = count * channel_rank
        _, base = orthogonal_frames(count, channel_rank)
        for overlap in (0.9, 0.99, 0.999, 0.999999, 1.0):
            frames = [frame.copy() for frame in base]
            frames[1] = (
                overlap * base[0]
                + math.sqrt(max(0.0, 1.0 - overlap**2)) * base[1]
            )
            gram = block_gram(frames)
            deviation = float(np.linalg.norm(gram - np.eye(dimension), ord=2))
            negative, _, _ = inertia(stein_rhs(frames), 1.0e-8)
            expected = channel_rank * (count - 1 if overlap < 1.0 else count - 2)
            passed = bool(abs(deviation - overlap) < 2.0e-10 and negative == expected)
            failures += int(not passed)
            cases.append(
                {
                    "nodes": count,
                    "channel_rank": channel_rank,
                    "overlap": overlap,
                    "gram_deviation": deviation,
                    "negative_inertia": negative,
                    "expected_negative_inertia": expected,
                    "all_pass": passed,
                }
            )
    return {"cases": cases, "failures": failures, "all_pass": failures == 0}


def make_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = np.arange(2, 13)
    exact = counts - 1
    local = np.ceil(counts / 2.0)
    pairwise = np.ones_like(counts)
    spectrum_count = 9
    completion, _ = sharp_completion(regular_nodes(spectrum_count), 1)
    spectrum = np.linalg.eigvalsh(completion)

    figure, axes = plt.subplots(1, 2, figsize=(8.5, 3.05))
    axes[0].plot(counts, exact, "o-", lw=2.0, label=r"global exact $L-1$")
    axes[0].plot(counts, local, "s--", lw=1.7, label=r"endpoint bound $\lceil L/2\rceil$")
    axes[0].plot(counts, pairwise, ":", lw=2.0, label="largest pairwise cost")
    axes[0].set_xlabel("number of frequency nodes $L$")
    axes[0].set_ylabel("minimum degree / lower bound ($k=1$)")
    axes[0].grid(alpha=0.24)
    axes[0].legend(frameon=False, fontsize=7.8)

    axes[1].bar(np.arange(spectrum_count), spectrum, color="#277da1")
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set_xlabel("eigenvalue index")
    axes[1].set_ylabel("canonical Pick-completion eigenvalue")
    axes[1].set_title(r"regular $9$-gon: $0,1,\ldots,8$")
    axes[1].grid(axis="y", alpha=0.24)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/global_fanout_memory/certificate.json"),
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=Path("paper_global_fanout_memory/figures/global_fanout_memory.pdf"),
    )
    parser.add_argument("--node-trials", type=int, default=384)
    parser.add_argument("--noise-trials", type=int, default=1024)
    parser.add_argument("--near-orthogonal-trials", type=int, default=512)
    args = parser.parse_args()

    payload = {
        "schema_version": 1,
        "artifact": "global-fanout-memory-certificate",
        "deterministic_seeds": {"nodes": 26081041, "noise": 26081042, "near_orthogonal": 26081043},
        "regular_polygon": regular_polygon_audit(24),
        "arbitrary_nodes": arbitrary_node_audit(args.node_trials, 26081041),
        "clustered_nodes": clustered_node_audit(),
        "colligation_synthesis": colligation_campaign(),
        "noisy_inertia": noisy_inertia_audit(args.noise_trials, 26081042),
        "near_orthogonal": near_orthogonal_audit(args.near_orthogonal_trials, 26081043),
        "robust_boundary": robust_boundary_audit(),
    }
    payload["all_pass"] = bool(
        payload["regular_polygon"]["all_pass"]
        and payload["arbitrary_nodes"]["all_pass"]
        and payload["clustered_nodes"]["all_pass"]
        and payload["colligation_synthesis"]["all_pass"]
        and payload["noisy_inertia"]["all_pass"]
        and payload["near_orthogonal"]["all_pass"]
        and payload["robust_boundary"]["all_pass"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_figure(args.figure)
    if not payload["all_pass"]:
        raise SystemExit("global fan-out memory certificate: FAIL")
    print("global fan-out memory certificate: PASS")


if __name__ == "__main__":
    main()
