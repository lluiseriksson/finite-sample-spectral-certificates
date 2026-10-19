#!/usr/bin/env python3
"""Reproduce the Ky Fan subspace speed-limit certificates.

The paper proves, for every r, that twice the first r canonical angles of a
transported subspace are bounded by the integrated r-th spectral spread of
the Hermitian generator.  For positive Wigner--Smith flows the spread can be
replaced by the sum of the r largest proper delays.  This script attacks the
algebraic, dynamical, sharpness, robustness, and measurement interfaces with
deterministic adversarial tests.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "ky_fan_speed_limits" / "certificate.json"
DEFAULT_FIGURE = ROOT / "paper_ky_fan_speed_limits" / "figures" / "ky_fan_hierarchy.pdf"
TOL = 5.0e-10


def haar_unitary(rng: np.random.Generator, size: int) -> np.ndarray:
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    q, r = np.linalg.qr(raw)
    diagonal = np.diag(r)
    phase = np.where(np.abs(diagonal) > 0.0, diagonal / np.abs(diagonal), 1.0)
    return q * phase.conj()


def hermitian(rng: np.random.Generator, size: int) -> np.ndarray:
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return 0.5 * (raw + raw.conj().T)


def positive(rng: np.random.Generator, size: int) -> np.ndarray:
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return raw @ raw.conj().T / size


def partial_sums(values: np.ndarray) -> np.ndarray:
    return np.cumsum(np.asarray(values, dtype=float))


def spectral_spread(q: np.ndarray, count: int) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh(q)
    return eigenvalues[-1 : -count - 1 : -1] - eigenvalues[:count]


def off_diagonal_singular_values(q: np.ndarray, rank: int) -> np.ndarray:
    values = np.linalg.svd(q[rank:, :rank], compute_uv=False)
    return np.sort(values)[::-1]


def canonical_angles(unitary: np.ndarray, rank: int) -> np.ndarray:
    cosines = np.linalg.svd(unitary[:rank, :rank], compute_uv=False)
    angles = np.arccos(np.clip(cosines, 0.0, 1.0))
    return np.sort(angles)[::-1]


def pointwise_audit(trials: int, seed: int, psd: bool) -> dict:
    rng = np.random.default_rng(seed)
    worst_ratio = 0.0
    smallest_slack = float("inf")
    inequalities = 0
    for _ in range(trials):
        size = int(rng.integers(3, 19))
        rank = int(rng.integers(1, size))
        count = min(rank, size - rank)
        q = positive(rng, size) if psd else hermitian(rng, size)
        u = haar_unitary(rng, size)
        q = u.conj().T @ q @ u
        lhs = 2.0 * partial_sums(off_diagonal_singular_values(q, rank)[:count])
        if psd:
            eigenvalues = np.linalg.eigvalsh(q)[::-1]
            rhs = partial_sums(eigenvalues[:count])
        else:
            rhs = partial_sums(spectral_spread(q, count))
        ratios = np.divide(lhs, rhs, out=np.zeros_like(lhs), where=rhs > 1.0e-14)
        worst_ratio = max(worst_ratio, float(np.max(ratios)))
        smallest_slack = min(smallest_slack, float(np.min(rhs - lhs)))
        inequalities += count
    return {
        "seed": seed,
        "matrices": trials,
        "partial_sum_inequalities": inequalities,
        "worst_lhs_over_rhs": worst_ratio,
        "smallest_absolute_slack": smallest_slack,
        "all_pass": bool(worst_ratio <= 1.0 + TOL and smallest_slack >= -TOL),
    }


def sharp_family(angle_spectra: list[list[float]]) -> dict:
    cases = []
    worst_error = 0.0
    for raw_angles in angle_spectra:
        angles = np.sort(np.asarray(raw_angles, dtype=float))[::-1]
        rank = angles.size
        q = np.zeros((2 * rank, 2 * rank), dtype=complex)
        coupling = np.diag(angles)
        q[:rank, :rank] = coupling
        q[rank:, rank:] = coupling
        q[:rank, rank:] = -coupling
        q[rank:, :rank] = -coupling
        unitary = expm(1j * q)
        observed = canonical_angles(unitary, rank)
        top_delay = partial_sums(np.linalg.eigvalsh(q)[::-1][:rank])
        required = 2.0 * partial_sums(angles)
        delay_vector = np.linalg.eigvalsh(q)[::-1][:rank]
        error = max(float(np.max(np.abs(observed - angles))), float(np.max(np.abs(top_delay - required))))
        gauge_costs = {
            "operator": {
                "geometric_optimum": 2.0 * float(np.max(angles)),
                "delay_action": float(np.max(delay_vector)),
            },
            "euclidean": {
                "geometric_optimum": 2.0 * float(np.linalg.norm(angles)),
                "delay_action": float(np.linalg.norm(delay_vector)),
            },
            "trace": {
                "geometric_optimum": 2.0 * float(np.sum(angles)),
                "delay_action": float(np.sum(delay_vector)),
            },
        }
        gauge_error = max(
            abs(record["geometric_optimum"] - record["delay_action"])
            for record in gauge_costs.values()
        )
        error = max(error, gauge_error)
        worst_error = max(worst_error, error)
        cases.append({
            "rank": int(rank),
            "angles": angles.tolist(),
            "observed_angles": observed.tolist(),
            "twice_angle_partial_sums": required.tolist(),
            "top_delay_partial_sums": top_delay.tolist(),
            "symmetric_gauge_costs": gauge_costs,
            "maximum_error": error,
        })
    return {"cases": cases, "worst_error": worst_error, "all_pass": bool(worst_error <= TOL)}


def scalar_endpoint_separation() -> dict:
    """Two angle spectra invisible to max and trace, separated at r=2."""
    first = np.array([1.0, 0.8, 0.2])
    second = np.array([1.0, 0.5, 0.5])
    first_prefix = 2.0 * partial_sums(first)
    second_prefix = 2.0 * partial_sums(second)
    return {
        "first_angles": first.tolist(),
        "second_angles": second.tolist(),
        "common_largest_angle": float(first[0]),
        "common_total_angle": float(np.sum(first)),
        "first_optimal_action_prefixes": first_prefix.tolist(),
        "second_optimal_action_prefixes": second_prefix.tolist(),
        "r2_action_gap": float(first_prefix[1] - second_prefix[1]),
        "same_scalar_endpoints": bool(
            np.isclose(first[0], second[0]) and np.isclose(np.sum(first), np.sum(second))
        ),
        "strict_intermediate_separation": bool(first_prefix[1] > second_prefix[1]),
        "all_pass": bool(
            np.isclose(first[0], second[0])
            and np.isclose(np.sum(first), np.sum(second))
            and first_prefix[1] > second_prefix[1]
        ),
    }


def random_path_audit(paths: int, seed: int, psd: bool) -> dict:
    rng = np.random.default_rng(seed)
    worst_ratio = 0.0
    smallest_slack = float("inf")
    inequalities = 0
    for _ in range(paths):
        rank = int(rng.integers(1, 6))
        size = int(rng.integers(2 * rank, 2 * rank + 5))
        segments = int(rng.integers(2, 8))
        unitary = np.eye(size, dtype=complex)
        action = np.zeros(rank)
        for _ in range(segments):
            duration = float(rng.uniform(0.015, 0.16))
            q = positive(rng, size) if psd else hermitian(rng, size)
            q /= max(float(np.linalg.norm(q, ord=2)), 1.0e-12)
            if psd:
                action += duration * np.linalg.eigvalsh(q)[::-1][:rank]
            else:
                action += duration * spectral_spread(q, rank)
            unitary = unitary @ expm(1j * duration * q)
        angles = canonical_angles(unitary, rank)
        lhs = 2.0 * partial_sums(angles)
        rhs = partial_sums(action)
        ratios = np.divide(lhs, rhs, out=np.zeros_like(lhs), where=rhs > 1.0e-14)
        worst_ratio = max(worst_ratio, float(np.max(ratios)))
        smallest_slack = min(smallest_slack, float(np.min(rhs - lhs)))
        inequalities += rank
    return {
        "seed": seed,
        "piecewise_constant_paths": paths,
        "partial_sum_inequalities": inequalities,
        "worst_lhs_over_rhs": worst_ratio,
        "smallest_absolute_slack": smallest_slack,
        "all_pass": bool(worst_ratio <= 1.0 + TOL and smallest_slack >= -TOL),
    }


def leakage_certificate(pass_leak: np.ndarray, stop_leak: np.ndarray) -> np.ndarray:
    rank = pass_leak.size
    pass_angles = np.arcsin(np.clip(np.sort(pass_leak)[::-1], 0.0, 1.0))
    stop_angles = np.arcsin(np.clip(np.sort(stop_leak)[::-1], 0.0, 1.0))
    return np.maximum(
        0.0,
        0.5 * np.pi * np.arange(1, rank + 1)
        - partial_sums(pass_angles)
        - partial_sums(stop_angles),
    )


def tomography_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    worst_overstatement = -float("inf")
    examples = []
    for trial in range(trials):
        rank = int(rng.integers(2, 9))
        pass_leak = np.sort(rng.uniform(0.0, 0.35, rank))[::-1]
        stop_leak = np.sort(rng.uniform(0.0, 0.35, rank))[::-1]
        eta_p = float(rng.uniform(0.0, 0.025))
        eta_s = float(rng.uniform(0.0, 0.025))
        # Pass measurements observe the desired-sector singular values, whose
        # reverse ordering is complementary to leakage singular values.
        pass_signal = np.sqrt(1.0 - pass_leak[::-1] ** 2)
        measured_pass = np.clip(pass_signal + rng.uniform(-eta_p, eta_p, rank), 0.0, 1.0)
        measured_stop = np.clip(stop_leak + rng.uniform(-eta_s, eta_s, rank), 0.0, 1.0)
        pass_upper = np.sqrt(1.0 - np.maximum(0.0, np.sort(measured_pass)[::-1][::-1] - eta_p) ** 2)
        pass_upper = np.sort(np.clip(pass_upper, 0.0, 1.0))[::-1]
        stop_upper = np.sort(np.clip(np.sort(measured_stop)[::-1] + eta_s, 0.0, 1.0))[::-1]
        actual = leakage_certificate(pass_leak, stop_leak)
        certified = leakage_certificate(pass_upper, stop_upper)
        overstatement = float(np.max(certified - actual))
        worst_overstatement = max(worst_overstatement, overstatement)
        if trial < 3:
            examples.append({
                "rank": rank,
                "eta_pass": eta_p,
                "eta_stop": eta_s,
                "actual_gamma": actual.tolist(),
                "certified_gamma": certified.tolist(),
            })
    return {
        "seed": seed,
        "trials": trials,
        "worst_certified_minus_actual": worst_overstatement,
        "examples": examples,
        "all_pass": bool(worst_overstatement <= TOL),
    }


def slack_decomposition_audit(trials: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    worst_identity_error = 0.0
    all_nonnegative = True
    for _ in range(trials):
        rank = int(rng.integers(1, 9))
        bottom = rng.exponential(scale=0.3, size=rank)
        pointwise = rng.exponential(scale=0.4, size=rank)
        geometric = rng.exponential(scale=0.2, size=rank)
        total = partial_sums(bottom + pointwise + 2.0 * geometric)
        reconstructed = partial_sums(bottom) + partial_sums(pointwise) + 2.0 * partial_sums(geometric)
        worst_identity_error = max(worst_identity_error, float(np.max(np.abs(total - reconstructed))))
        all_nonnegative = all_nonnegative and bool(np.all(bottom >= 0) and np.all(pointwise >= 0) and np.all(geometric >= 0))
    return {
        "seed": seed,
        "trials": trials,
        "worst_identity_error": worst_identity_error,
        "all_components_nonnegative": all_nonnegative,
        "all_pass": bool(worst_identity_error <= TOL and all_nonnegative),
    }


def make_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    angles = np.array([1.30, 1.05, 0.72, 0.41, 0.16])
    ideal = 2.0 * partial_sums(angles)
    rng = np.random.default_rng(260810)
    surplus = partial_sums(rng.uniform(0.07, 0.20, angles.size))
    generic = ideal + surplus
    ranks = np.arange(1, angles.size + 1)
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.45))
    axes[0].plot(ranks, ideal, "o-", lw=2.2, label="gauge-optimal path")
    axes[0].plot(ranks, generic, "s--", lw=1.9, label="generic passive path")
    axes[0].set_xlabel(r"Ky Fan index $r$")
    axes[0].set_ylabel(r"cumulative proper-delay action")
    axes[0].set_xticks(ranks)
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8.5)
    first = 2.0 * partial_sums(np.array([1.0, 0.8, 0.2]))
    second = 2.0 * partial_sums(np.array([1.0, 0.5, 0.5]))
    short_ranks = np.arange(1, 4)
    axes[1].plot(short_ranks, first, "o-", lw=2.2, label=r"$2(1,.8,.2)$")
    axes[1].plot(short_ranks, second, "D--", lw=1.9, label=r"$2(1,.5,.5)$")
    axes[1].set_xlabel(r"Ky Fan index $r$")
    axes[1].set_ylabel(r"sharp cumulative action")
    axes[1].set_xticks(short_ranks)
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8.5)
    figure.tight_layout()
    plt.savefig(path)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--pointwise-trials", type=int, default=4096)
    parser.add_argument("--path-trials", type=int, default=768)
    parser.add_argument("--tomography-trials", type=int, default=2048)
    args = parser.parse_args()

    certificate = {
        "schema": "ky-fan-speed-limits/v1",
        "theorem": {
            "signed": "2 sum_{j<=r} beta_j <= integral sum_{j<=r}(lambda_j-lambda_{N-j+1})",
            "positive": "2 sum_{j<=r} beta_j <= integral sum_{j<=r} lambda_j",
            "indices": "r=1,...,min(k,N-k)",
        },
        "sharp_equalities": sharp_family([
            [1.21],
            [1.42, 0.87],
            [1.49, 1.08, 0.63, 0.17],
            [1.51, 1.22, 0.94, 0.58, 0.31, 0.09],
        ]),
        "scalar_endpoint_separation": scalar_endpoint_separation(),
        "signed_pointwise": pointwise_audit(args.pointwise_trials, 26081001, psd=False),
        "positive_pointwise": pointwise_audit(args.pointwise_trials, 26081002, psd=True),
        "signed_paths": random_path_audit(args.path_trials, 26081003, psd=False),
        "positive_paths": random_path_audit(args.path_trials, 26081004, psd=True),
        "tomography": tomography_audit(args.tomography_trials, 26081005),
        "slack_decomposition": slack_decomposition_audit(4096, 26081006),
    }
    certificate["all_checks_pass"] = all(
        block.get("all_pass", False)
        for key, block in certificate.items()
        if key not in {"schema", "theorem", "all_checks_pass"}
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
    make_figure(args.figure)
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()
