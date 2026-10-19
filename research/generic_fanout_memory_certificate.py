#!/usr/bin/env python3
"""Deterministic certificates for generic spectral fan-out memory.

The producer audits the universal block-Schur Pick completion, the exact
full-span degree law, near-coincident direct-sum targets, representative
minimal conservative realizations, and fail-closed noisy span counts.
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


TOL = 3.0e-8


def random_nodes(rng: np.random.Generator, count: int) -> np.ndarray:
    while True:
        angles = np.sort(rng.uniform(0.0, 2.0 * math.pi, count))
        gaps = np.diff(np.r_[angles, angles[0] + 2.0 * math.pi])
        if float(np.min(gaps)) > 0.04:
            return np.exp(1j * angles)


def random_frames(
    rng: np.random.Generator, count: int, channel_rank: int
) -> list[np.ndarray]:
    dimension = count * channel_rank
    frames = []
    for _ in range(count):
        sample = rng.normal(size=(dimension, channel_rank)) + 1j * rng.normal(
            size=(dimension, channel_rank)
        )
        q, _ = np.linalg.qr(sample)
        frames.append(q[:, :channel_rank])
    return frames


def near_coincident_frames(
    count: int, channel_rank: int, epsilon: float
) -> list[np.ndarray]:
    """L direct-sum k-planes converging to one common k-plane."""
    dimension = (count + 1) * channel_rank
    frames = []
    scale = math.sqrt(1.0 + epsilon**2)
    for index in range(count):
        frame = np.zeros((dimension, channel_rank), dtype=complex)
        frame[:channel_rank, :] = np.eye(channel_rank)
        frame[(index + 1) * channel_rank : (index + 2) * channel_rank, :] = (
            epsilon * np.eye(channel_rank)
        )
        frames.append(frame / scale)
    return frames


def block_gram(frames: list[np.ndarray]) -> np.ndarray:
    return np.block(
        [[left.conj().T @ right for right in frames] for left in frames]
    )


def pick_off_diagonal(nodes: np.ndarray, frames: list[np.ndarray]) -> np.ndarray:
    count = len(nodes)
    channel_rank = frames[0].shape[1]
    answer = np.zeros((count * channel_rank, count * channel_rank), dtype=complex)
    for row in range(count):
        for column in range(count):
            if row == column:
                continue
            numerator = np.eye(channel_rank) - frames[row].conj().T @ frames[column]
            denominator = 1.0 - np.conj(nodes[row]) * nodes[column]
            answer[
                row * channel_rank : (row + 1) * channel_rank,
                column * channel_rank : (column + 1) * channel_rank,
            ] = numerator / denominator
    return answer


def block_schur_completion(
    nodes: np.ndarray, frames: list[np.ndarray]
) -> tuple[np.ndarray, float]:
    """Complete all free diagonal blocks with rank at most k(L-1)."""
    partial = pick_off_diagonal(nodes, frames)
    count = len(nodes)
    channel_rank = frames[0].shape[1]
    leading_dimension = (count - 1) * channel_rank
    leading = partial[:leading_dimension, :leading_dimension]
    shift = max(1.0, -float(np.linalg.eigvalsh(leading)[0]) + 1.0)
    q = leading + shift * np.eye(leading_dimension)
    b = partial[:leading_dimension, leading_dimension:]
    last = b.conj().T @ np.linalg.solve(q, b)
    completion = np.block([[q, b], [b.conj().T, last]])
    return completion, shift


def stein_residual(
    nodes: np.ndarray, frames: list[np.ndarray], completion: np.ndarray
) -> np.ndarray:
    channel_rank = frames[0].shape[1]
    z = np.kron(np.diag(nodes), np.eye(channel_rank))
    rhs = np.kron(np.ones((len(nodes), len(nodes))), np.eye(channel_rank)) - block_gram(frames)
    return completion - z.conj().T @ completion @ z - rhs


def canonical_complement(isometry: np.ndarray, tolerance: float = 1.0e-11) -> np.ndarray:
    rows, columns = isometry.shape
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
        if len(accepted) == rows - columns:
            break
    if len(accepted) != rows - columns:
        raise RuntimeError("failed to construct deterministic complement")
    return np.column_stack(accepted) if accepted else np.zeros((rows, 0), dtype=complex)


def synthesize_colligation(
    nodes: np.ndarray, frames: list[np.ndarray], completion: np.ndarray
) -> dict:
    count = len(nodes)
    channel_rank = frames[0].shape[1]
    output_dimension = frames[0].shape[0]
    input_frame = np.eye(output_dimension, channel_rank, dtype=complex)
    values, vectors = np.linalg.eigh(completion)
    scale = max(1.0, float(np.max(abs(values))))
    positive = values > 1.0e-9 * scale
    factor = np.diag(np.sqrt(values[positive])) @ vectors[:, positive].conj().T
    state_dimension = int(factor.shape[0])

    source_columns = []
    target_columns = []
    for index, node in enumerate(nodes):
        feature = factor[:, index * channel_rank : (index + 1) * channel_rank]
        source_columns.append(np.vstack([node * feature, input_frame]))
        target_columns.append(np.vstack([feature, frames[index]]))
    source = np.hstack(source_columns)
    target = np.hstack(target_columns)
    gram_error = float(np.linalg.norm(source.conj().T @ source - target.conj().T @ target, ord=2))

    sample_gram = source.conj().T @ source
    gram_values, gram_vectors = np.linalg.eigh(sample_gram)
    inverse_sqrt = gram_vectors @ np.diag(1.0 / np.sqrt(gram_values)) @ gram_vectors.conj().T
    q_source = source @ inverse_sqrt
    q_target = target @ inverse_sqrt
    source_full = np.hstack([q_source, canonical_complement(q_source)])
    target_full = np.hstack([q_target, canonical_complement(q_target)])
    colligation = target_full @ source_full.conj().T

    r = state_dimension
    a = colligation[:r, :r]
    b = colligation[:r, r:]
    c = colligation[r:, :r]
    d = colligation[r:, r:]
    controllability = np.hstack([np.linalg.matrix_power(a, power) @ b for power in range(r)])
    observability = np.vstack([c @ np.linalg.matrix_power(a, power) for power in range(r)])
    interpolation_errors = []
    unitarity_errors = []
    for index, node in enumerate(nodes):
        transfer = d + node * c @ np.linalg.solve(np.eye(r) - node * a, b)
        mapped = transfer @ input_frame
        interpolation_errors.append(
            float(
                np.linalg.norm(
                    mapped @ mapped.conj().T
                    - frames[index] @ frames[index].conj().T,
                    ord=2,
                )
            )
        )
        unitarity_errors.append(
            float(np.linalg.norm(transfer.conj().T @ transfer - np.eye(output_dimension), ord=2))
        )
    grid_unitarity_errors = []
    for theta in np.linspace(0.0, 2.0 * math.pi, 257, endpoint=False):
        point = np.exp(1j * theta)
        transfer = d + point * c @ np.linalg.solve(np.eye(r) - point * a, b)
        grid_unitarity_errors.append(
            float(np.linalg.norm(transfer.conj().T @ transfer - np.eye(output_dimension), ord=2))
        )
    return {
        "state_dimension": state_dimension,
        "sample_gram_error": gram_error,
        "colligation_unitarity_error": float(
            np.linalg.norm(
                colligation.conj().T @ colligation - np.eye(colligation.shape[0]), ord=2
            )
        ),
        "maximum_interpolation_error": max(interpolation_errors),
        "maximum_boundary_unitarity_error": max(unitarity_errors),
        "maximum_grid_unitarity_error": max(grid_unitarity_errors),
        "state_spectral_radius": float(max(abs(np.linalg.eigvals(a)))) if r else 0.0,
        "controllability_rank": int(np.linalg.matrix_rank(controllability, tol=1.0e-8)),
        "observability_rank": int(np.linalg.matrix_rank(observability, tol=1.0e-8)),
    }


def universal_completion_campaign(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    failures = 0
    worst_psd_violation = 0.0
    worst_stein_residual = 0.0
    smallest_frame_singular_value = math.inf
    samples = []
    for trial in range(trials):
        count = int(rng.integers(3, 9))
        channel_rank = int(rng.integers(1, 4))
        nodes = random_nodes(rng, count)
        frames = random_frames(rng, count, channel_rank)
        completion, shift = block_schur_completion(nodes, frames)
        eigenvalues = np.linalg.eigvalsh(completion)
        scale = max(1.0, float(np.max(abs(eigenvalues))))
        observed_rank = int(np.count_nonzero(eigenvalues > 1.0e-9 * scale))
        expected_rank = channel_rank * (count - 1)
        residual = float(np.linalg.norm(stein_residual(nodes, frames, completion), ord=2))
        stacked = np.hstack(frames)
        sigma_min = float(np.linalg.svd(stacked, compute_uv=False)[-1])
        passed = bool(
            observed_rank == expected_rank
            and eigenvalues[0] >= -2.0e-8 * scale
            and residual < 2.0e-8 * scale
            and sigma_min > 1.0e-8
        )
        failures += int(not passed)
        worst_psd_violation = max(worst_psd_violation, float(max(0.0, -eigenvalues[0] / scale)))
        worst_stein_residual = max(worst_stein_residual, residual / scale)
        smallest_frame_singular_value = min(smallest_frame_singular_value, sigma_min)
        if trial < 10:
            samples.append(
                {
                    "nodes": count,
                    "channel_rank": channel_rank,
                    "completion_rank": observed_rank,
                    "exact_degree": expected_rank,
                    "diagonal_shift": shift,
                    "minimum_stacked_frame_singular_value": sigma_min,
                }
            )
    return {
        "trials": trials,
        "failures": failures,
        "worst_relative_psd_violation": worst_psd_violation,
        "worst_relative_stein_residual": worst_stein_residual,
        "smallest_stacked_frame_singular_value": smallest_frame_singular_value,
        "sample": samples,
        "all_pass": failures == 0,
    }


def near_collapse_campaign() -> dict:
    cases = []
    failures = 0
    worst_spectrum_error = 0.0
    for count, channel_rank in ((3, 1), (5, 1), (4, 2), (6, 2)):
        nodes = np.exp(2j * math.pi * np.arange(count) / count)
        for epsilon in (1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3):
            frames = near_coincident_frames(count, channel_rank, epsilon)
            gram = block_gram(frames)
            gram_values = np.linalg.eigvalsh(gram)
            predicted_minimum = epsilon**2 / (1.0 + epsilon**2)
            spectrum_error = float(abs(gram_values[0] - predicted_minimum))
            completion, _ = block_schur_completion(nodes, frames)
            completion_values = np.linalg.eigvalsh(completion)
            scale = max(1.0, float(np.max(abs(completion_values))))
            completion_rank = int(
                np.count_nonzero(completion_values > 1.0e-9 * scale)
            )
            expected = channel_rank * (count - 1)
            maximum_overlap = max(
                float(np.linalg.norm(frames[i].conj().T @ frames[j], ord=2))
                for i in range(count)
                for j in range(i + 1, count)
            )
            passed = bool(
                completion_rank == expected
                and gram_values[0] > 0.0
                and maximum_overlap > 0.99
                and spectrum_error < 2.0e-12
            )
            failures += int(not passed)
            worst_spectrum_error = max(worst_spectrum_error, spectrum_error)
            cases.append(
                {
                    "nodes": count,
                    "channel_rank": channel_rank,
                    "epsilon": epsilon,
                    "maximum_pairwise_overlap": maximum_overlap,
                    "minimum_gram_eigenvalue": float(gram_values[0]),
                    "predicted_minimum_gram_eigenvalue": predicted_minimum,
                    "spectrum_error": spectrum_error,
                    "completion_rank": completion_rank,
                    "exact_degree": expected,
                    "all_pass": passed,
                }
            )
    return {
        "cases": cases,
        "failures": failures,
        "worst_spectrum_error": worst_spectrum_error,
        "all_pass": failures == 0,
    }


def colligation_campaign(seed: int) -> dict:
    rng = np.random.default_rng(seed)
    cases = []
    all_pass = True
    for count, channel_rank in ((3, 1), (4, 1), (3, 2), (4, 2), (5, 2)):
        nodes = random_nodes(rng, count)
        frames = random_frames(rng, count, channel_rank)
        completion, _ = block_schur_completion(nodes, frames)
        synthesis = synthesize_colligation(nodes, frames, completion)
        expected = channel_rank * (count - 1)
        passed = bool(
            synthesis["state_dimension"] == expected
            and synthesis["controllability_rank"] == expected
            and synthesis["observability_rank"] == expected
            and synthesis["sample_gram_error"] < 3.0e-8
            and synthesis["colligation_unitarity_error"] < 3.0e-8
            and synthesis["maximum_interpolation_error"] < 3.0e-7
            and synthesis["maximum_boundary_unitarity_error"] < 3.0e-7
            and synthesis["maximum_grid_unitarity_error"] < 3.0e-7
            and synthesis["state_spectral_radius"] < 1.0
        )
        synthesis.update(
            {
                "nodes": count,
                "channel_rank": channel_rank,
                "pairwise_degree": channel_rank,
                "global_degree": expected,
                "all_pass": passed,
            }
        )
        cases.append(synthesis)
        all_pass = all_pass and passed
    return {"cases": cases, "all_pass": all_pass}


def noisy_projector_sum_campaign(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    overcertificates = 0
    direct_sum_certificates = 0
    direct_sum_trials = 0
    for trial in range(trials):
        count = int(rng.integers(3, 9))
        channel_rank = int(rng.integers(1, 4))
        full = trial < trials // 2
        if full:
            frames = random_frames(rng, count, channel_rank)
        else:
            ambient = max(channel_rank, (count - 1) * channel_rank)
            frames = []
            for _ in range(count):
                sample = rng.normal(size=(ambient, channel_rank)) + 1j * rng.normal(
                    size=(ambient, channel_rank)
                )
                q, _ = np.linalg.qr(sample)
                padded = np.zeros((count * channel_rank, channel_rank), dtype=complex)
                padded[:ambient, :] = q[:, :channel_rank]
                frames.append(padded)
        projector_sum = sum(frame @ frame.conj().T for frame in frames)
        true_rank = int(np.linalg.matrix_rank(projector_sum, tol=1.0e-9))
        estimated_projectors = []
        eta = 0.0
        for frame in frames:
            perturbation = 1.0e-8 * (
                rng.normal(size=frame.shape) + 1j * rng.normal(size=frame.shape)
            )
            estimated_frame, _ = np.linalg.qr(frame + perturbation)
            estimated_frame = estimated_frame[:, :channel_rank]
            exact_projector = frame @ frame.conj().T
            estimated_projector = estimated_frame @ estimated_frame.conj().T
            estimated_projectors.append(estimated_projector)
            eta += float(np.linalg.norm(estimated_projector - exact_projector, ord=2))
        observed = sum(estimated_projectors)
        certified_rank = int(np.count_nonzero(np.linalg.eigvalsh(observed) > eta))
        certified_degree = max(0, certified_rank - channel_rank)
        true_span_bound = max(0, true_rank - channel_rank)
        overcertificates += int(certified_degree > true_span_bound)
        if full:
            direct_sum_trials += 1
            direct_sum_certificates += int(certified_degree == channel_rank * (count - 1))
    return {
        "trials": trials,
        "direct_sum_trials": direct_sum_trials,
        "direct_sum_certificates": direct_sum_certificates,
        "overcertificates": overcertificates,
        "all_pass": bool(
            overcertificates == 0
            and direct_sum_certificates == direct_sum_trials
        ),
    }


def clustered_node_campaign(seed: int) -> dict:
    rng = np.random.default_rng(seed)
    frames = random_frames(rng, 4, 1)
    cases = []
    failures = 0
    for gap in (1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8):
        nodes = np.exp(1j * np.array([0.0, gap, 1.7, 4.2]))
        completion, _ = block_schur_completion(nodes, frames)
        values = np.linalg.eigvalsh(completion)
        scale = max(1.0, float(np.max(abs(values))))
        residual = float(np.linalg.norm(stein_residual(nodes, frames, completion), ord=2)) / scale
        numeric_rank = int(np.count_nonzero(values > 1.0e-9 * scale))
        leading = completion[:3, :3]
        condition = float(np.linalg.cond(leading))
        passed = bool(values[0] >= -2.0e-8 * scale and residual < 2.0e-8)
        failures += int(not passed)
        cases.append(
            {
                "minimum_node_gap": gap,
                "completion_scale": scale,
                "leading_condition_number": condition,
                "relative_stein_residual": residual,
                "analytic_rank": 3,
                "floating_rank_at_1e-9_relative": numeric_rank,
                "all_pass": passed,
            }
        )
    return {"cases": cases, "failures": failures, "all_pass": failures == 0}


def make_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    epsilons = np.logspace(-4, -0.25, 120)
    count = 6
    minimum_eigenvalues = []
    overlaps = []
    for epsilon in epsilons:
        frames = near_coincident_frames(count, 1, float(epsilon))
        minimum_eigenvalues.append(float(np.linalg.eigvalsh(block_gram(frames))[0]))
        overlaps.append(
            max(float(abs(frames[i].conj().T @ frames[j]).item()) for i in range(count) for j in range(i + 1, count))
        )

    figure, axes = plt.subplots(1, 2, figsize=(8.6, 3.15))
    left = axes[0]
    left.semilogx(epsilons, overlaps, color="#277da1", lw=2.1)
    left.set_xlabel(r"direct-sum opening $\epsilon$")
    left.set_ylabel("largest target overlap", color="#277da1")
    left.tick_params(axis="y", labelcolor="#277da1")
    left.grid(alpha=0.23)
    twin = left.twinx()
    twin.loglog(epsilons, minimum_eigenvalues, color="#f3722c", lw=1.8)
    twin.set_ylabel(r"smallest Gram eigenvalue", color="#f3722c")
    twin.tick_params(axis="y", labelcolor="#f3722c")
    left.set_title(r"targets coalesce, exact degree stays $L-1=5$")

    counts = np.arange(2, 14)
    axes[1].plot(counts, counts - 1, "o-", lw=2.1, label="generic global degree")
    axes[1].plot(counts, np.ones_like(counts), "s--", lw=1.7, label="every pair")
    axes[1].set_xlabel("number of frequency nodes $L$")
    axes[1].set_ylabel(r"minimum degree ($k=1$)")
    axes[1].set_title("direct-sum pairwise/global separation")
    axes[1].grid(alpha=0.23)
    axes[1].legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/generic_fanout_memory/certificate.json"),
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=Path("paper_generic_fanout_memory/figures/generic_fanout_memory.pdf"),
    )
    args = parser.parse_args()

    payload = {
        "schema_version": 1,
        "theorem": {
            "universal_upper": "degree <= k(L-1)",
            "direct_sum_exact": "degree = k(L-1)",
            "span_lower": "degree >= dim_span(Y_1,...,Y_L)-k",
        },
        "universal_completion": universal_completion_campaign(512, 26081101),
        "near_collapse": near_collapse_campaign(),
        "colligation_synthesis": colligation_campaign(26081102),
        "noisy_projector_sum": noisy_projector_sum_campaign(1024, 26081103),
        "clustered_nodes": clustered_node_campaign(26081104),
    }
    payload["all_pass"] = all(
        payload[key]["all_pass"]
        for key in (
            "universal_completion",
            "near_collapse",
            "colligation_synthesis",
            "noisy_projector_sum",
            "clustered_nodes",
        )
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_figure(args.figure)
    if not payload["all_pass"]:
        raise RuntimeError("generic fan-out certificate failed")
    print("generic fan-out memory certificate: PASS")


if __name__ == "__main__":
    main()
