#!/usr/bin/env python3
"""Reproducible certificate for the direct-sum occupancy law.

The paper proves that a routing word w among rank-k target planes in
direct-sum position has exact passive memory

    d_min = k * (len(w) - min_a count_w(a)).

This producer checks the two load-bearing algebraic mechanisms without
optimizing a transfer matrix: dual block detectors for the lower bound and
the polynomial/matrix Fejer--Riesz precursor for the upper bound.  For
scalar targets it additionally computes the outer spectral factor and the
resulting rational-inner column.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "unified-routing-table-memory-certificate-v1"
SEED = 26081117


def degree_law(word: str, rank: int = 1) -> int:
    counts = Counter(word)
    if not counts or rank < 1:
        raise ValueError("nonempty word and positive rank required")
    return rank * (len(word) - min(counts.values()))


def node_polynomials(word: str, nodes: np.ndarray) -> tuple[list[str], list[np.ndarray]]:
    alphabet = sorted(set(word))
    output: list[np.ndarray] = []
    for symbol in alphabet:
        polynomial = np.array([1.0 + 0.0j])
        for label, node in zip(word, nodes, strict=True):
            if label != symbol:
                polynomial = np.polynomial.polynomial.polymul(
                    polynomial, np.array([-node, 1.0 + 0.0j])
                )
        output.append(np.asarray(polynomial, dtype=complex))
    return alphabet, output


def _pad(polynomial: np.ndarray, degree: int) -> np.ndarray:
    return np.pad(polynomial, (0, degree + 1 - polynomial.size))


def scalar_spectral_factor(component_polynomials: list[np.ndarray], degree: int) -> tuple[np.ndarray, float]:
    """Outer scalar factor of sum_j |q_j|^2 by reciprocal-root selection."""
    padded = [_pad(p, degree) for p in component_polynomials]
    palindromic = np.zeros(2 * degree + 1, dtype=complex)
    for polynomial in padded:
        palindromic += np.convolve(polynomial, np.conj(polynomial[::-1]))
    if degree:
        roots = np.roots(palindromic[::-1])
        exterior = roots[np.argsort(np.abs(roots))[::-1][:degree]]
        if np.min(np.abs(exterior)) <= 1.0 - 2.0e-8:
            raise RuntimeError("spectral root split failed")
        factor = np.poly(exterior)[::-1].astype(complex)
    else:
        factor = np.array([1.0 + 0.0j])
    grid = np.exp(1j * np.linspace(-math.pi, math.pi, 8192, endpoint=False))
    density = sum(
        np.abs(np.polynomial.polynomial.polyval(grid, p)) ** 2 for p in padded
    )
    values = np.polynomial.polynomial.polyval(grid, factor)
    factor *= math.sqrt(float(np.median(density / np.abs(values) ** 2)))
    values = np.polynomial.polynomial.polyval(grid, factor)
    residual = float(np.max(np.abs(np.abs(values) ** 2 - density)) / np.max(density))
    return factor, residual


def random_direct_sum_frames(
    rng: np.random.Generator, symbols: int, rank: int, ambient_padding: int = 2
) -> list[np.ndarray]:
    ambient = symbols * rank + ambient_padding
    frames: list[np.ndarray] = []
    while True:
        frames.clear()
        for _ in range(symbols):
            raw = rng.normal(size=(ambient, rank)) + 1j * rng.normal(size=(ambient, rank))
            q, _ = np.linalg.qr(raw)
            frames.append(q[:, :rank])
        stacked = np.hstack(frames)
        if np.linalg.matrix_rank(stacked, tol=1e-11) == symbols * rank:
            return list(frames)


def compile_case(word: str, angles: np.ndarray, frames: list[np.ndarray]) -> dict:
    nodes = np.exp(1j * np.asarray(angles, dtype=float))
    alphabet, polynomials = node_polynomials(word, nodes)
    if len(frames) != len(alphabet):
        raise ValueError("one target frame per used symbol required")
    rank = frames[0].shape[1]
    degree = len(word) - min(Counter(word).values())
    stacked = np.hstack(frames)
    left_inverse = np.linalg.solve(stacked.conj().T @ stacked, stacked.conj().T)
    dual_residual = float(np.linalg.norm(left_inverse @ stacked - np.eye(len(alphabet) * rank), 2))

    grid = np.exp(1j * np.linspace(-math.pi, math.pi, 1024, endpoint=False))
    min_density_eigenvalue = math.inf
    max_density_eigenvalue = 0.0
    for z in grid:
        f0 = sum(
            np.polynomial.polynomial.polyval(z, polynomial) * frame
            for polynomial, frame in zip(polynomials, frames, strict=True)
        )
        eigenvalues = np.linalg.eigvalsh(f0.conj().T @ f0)
        min_density_eigenvalue = min(min_density_eigenvalue, float(eigenvalues[0]))
        max_density_eigenvalue = max(max_density_eigenvalue, float(eigenvalues[-1]))

    range_error = 0.0
    wrong_dual_minor = 0.0
    desired_dual_minor_floor = math.inf
    for label, node in zip(word, nodes, strict=True):
        target = alphabet.index(label)
        f0 = sum(
            np.polynomial.polynomial.polyval(node, polynomial) * frame
            for polynomial, frame in zip(polynomials, frames, strict=True)
        )
        projector_f0 = f0 @ np.linalg.inv(f0.conj().T @ f0) @ f0.conj().T
        projector_target = frames[target] @ frames[target].conj().T
        range_error = max(range_error, float(np.linalg.norm(projector_f0 - projector_target, 2)))
        for detector in range(len(alphabet)):
            block = left_inverse[
                detector * rank : (detector + 1) * rank, :
            ] @ f0
            determinant = abs(np.linalg.det(block))
            if detector == target:
                desired_dual_minor_floor = min(desired_dual_minor_floor, float(determinant))
            else:
                wrong_dual_minor = max(wrong_dual_minor, float(determinant))

    output = {
        "word": word,
        "rank": rank,
        "length": len(word),
        "symbols": len(alphabet),
        "occupancies": dict(sorted(Counter(word).items())),
        "predicted_degree": rank * degree,
        "stack_smallest_singular_value": float(np.linalg.svd(stacked, compute_uv=False)[-1]),
        "dual_isolation_residual": dual_residual,
        "minimum_density_eigenvalue": min_density_eigenvalue,
        "maximum_density_eigenvalue": max_density_eigenvalue,
        "node_range_projector_error": range_error,
        "wrong_dual_minor": wrong_dual_minor,
        "desired_dual_minor_floor": desired_dual_minor_floor,
    }

    if rank == 1:
        padded = [_pad(p, degree) for p in polynomials]
        component_polynomials = []
        for row in range(stacked.shape[0]):
            component_polynomials.append(
                sum(frames[a][row, 0] * padded[a] for a in range(len(alphabet)))
            )
        factor, factor_residual = scalar_spectral_factor(component_polynomials, degree)
        inner_residual = 0.0
        scalar_range_error = 0.0
        for z in grid:
            f0 = sum(
                np.polynomial.polynomial.polyval(z, polynomial) * frame[:, 0]
                for polynomial, frame in zip(polynomials, frames, strict=True)
            )
            f = f0 / np.polynomial.polynomial.polyval(z, factor)
            inner_residual = max(inner_residual, abs(float(np.vdot(f, f).real) - 1.0))
        for label, node in zip(word, nodes, strict=True):
            f0 = sum(
                np.polynomial.polynomial.polyval(node, polynomial) * frame[:, 0]
                for polynomial, frame in zip(polynomials, frames, strict=True)
            )
            f = f0 / np.polynomial.polynomial.polyval(node, factor)
            target = frames[alphabet.index(label)][:, 0]
            f_unit = f / np.linalg.norm(f)
            scalar_range_error = max(
                scalar_range_error,
                float(np.linalg.norm(np.outer(f_unit, f_unit.conj()) - np.outer(target, target.conj()), 2)),
            )
        output.update(
            {
                "spectral_factor_relative_residual": factor_residual,
                "inner_column_residual": float(inner_residual),
                "normalized_node_projector_error": scalar_range_error,
            }
        )
    return output


def random_word(rng: np.random.Generator, length: int, symbols: int) -> str:
    data = np.concatenate([np.arange(symbols), rng.integers(0, symbols, length - symbols)])
    rng.shuffle(data)
    return "".join(chr(ord("A") + int(x)) for x in data)


def collision_family(epsilon: float, symbols: int = 4) -> list[np.ndarray]:
    frames = []
    for index in range(symbols):
        vector = np.zeros(symbols + 1, dtype=complex)
        vector[0] = 1.0
        vector[index + 1] = epsilon
        vector /= np.linalg.norm(vector)
        frames.append(vector[:, None])
    return frames


def make_certificate(scalar_trials: int = 512, block_trials: int = 256) -> dict:
    rng = np.random.default_rng(SEED)
    scalar_cases = []
    for _ in range(scalar_trials):
        length = int(rng.integers(2, 10))
        symbols = int(rng.integers(1, min(6, length + 1)))
        word = random_word(rng, length, symbols)
        angles = np.sort(rng.uniform(-math.pi, math.pi, length))
        scalar_cases.append(compile_case(word, angles, random_direct_sum_frames(rng, symbols, 1)))

    block_cases = []
    for _ in range(block_trials):
        length = int(rng.integers(2, 9))
        symbols = int(rng.integers(1, min(5, length + 1)))
        rank = int(rng.integers(2, 4))
        word = random_word(rng, length, symbols)
        angles = np.sort(rng.uniform(-math.pi, math.pi, length))
        block_cases.append(compile_case(word, angles, random_direct_sum_frames(rng, symbols, rank)))

    collision = []
    for epsilon in [1.0, 0.3, 0.1, 0.03, 0.01, 0.003]:
        case = compile_case(
            "AABCD",
            np.array([-2.5, -1.3, -0.2, 0.9, 2.4]),
            collision_family(epsilon),
        )
        case["epsilon"] = epsilon
        collision.append(case)

    all_cases = scalar_cases + block_cases + collision
    summary = {
        "scalar_trials": scalar_trials,
        "block_trials": block_trials,
        "maximum_dual_isolation_residual": max(x["dual_isolation_residual"] for x in all_cases),
        "maximum_node_range_projector_error": max(x["node_range_projector_error"] for x in all_cases),
        "maximum_wrong_dual_minor": max(x["wrong_dual_minor"] for x in all_cases),
        "minimum_density_eigenvalue": min(x["minimum_density_eigenvalue"] for x in all_cases),
        "maximum_scalar_factor_residual": max(x["spectral_factor_relative_residual"] for x in scalar_cases + collision),
        "maximum_scalar_inner_residual": max(x["inner_column_residual"] for x in scalar_cases + collision),
        "maximum_scalar_node_projector_error": max(x["normalized_node_projector_error"] for x in scalar_cases + collision),
    }
    payload = {
        "schema": SCHEMA,
        "seed": SEED,
        "theorem": "d_min = k * (L - min_a n_a) on the direct-sum target locus",
        "scope": "distinct boundary nodes; exact subspace routing; equal target rank; square finite rational-inner completion; regular node values; direct-sum target planes",
        "summary": summary,
        "collision_family": collision,
        "curated_degree_examples": {
            word: degree_law(word) for word in ["AAAA", "AAAB", "AABB", "ABAB", "AABC", "ABCD", "AABBCC"]
        },
        "claims_not_tested": [
            "floating-point campaigns do not prove minimality",
            "the general matrix Fejer--Riesz theorem and degree-preserving lossless completion are cited analytic results",
            "near-collision conditioning is not a counterexample to the exact theorem",
        ],
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["content_digest"] = hashlib.sha256(canonical).hexdigest()
    payload["all_pass"] = bool(
        summary["maximum_dual_isolation_residual"] < 5e-8
        and summary["maximum_node_range_projector_error"] < 2e-8
        and summary["maximum_wrong_dual_minor"] < 2e-8
        and summary["minimum_density_eigenvalue"] > 1e-10
        and summary["maximum_scalar_factor_residual"] < 5e-6
        and summary["maximum_scalar_inner_residual"] < 5e-5
        and summary["maximum_scalar_node_projector_error"] < 5e-8
    )
    return payload


def make_figure(certificate: dict, path: Path) -> None:
    words = ["AAAA", "AAAB", "AABB", "ABAB", "ABCD"]
    degrees = [degree_law(word) for word in words]
    collision = certificate["collision_family"]
    eps = np.array([x["epsilon"] for x in collision])
    margins = np.array([x["stack_smallest_singular_value"] ** 2 for x in collision])

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.65))
    colors = ["#8fa6b8", "#d75d57", "#4c78a8", "#59a14f", "#ad76a5"]
    axes[0].bar(words, degrees, color=colors)
    axes[0].set_ylabel(r"exact scalar degree $d_{\min}$")
    axes[0].set_title("The rarest target fixes memory")
    axes[0].set_ylim(0, 3.35)
    for i, value in enumerate(degrees):
        axes[0].text(i, value + 0.08, str(value), ha="center", fontsize=9)

    axes[1].loglog(eps, margins, "o-", color="#d75d57", label="direct-sum margin")
    axes[1].axhline(1e-16, color="black", linewidth=0.8, alpha=0.4)
    axes[1].set_xlabel(r"target opening $\epsilon$")
    axes[1].set_ylabel(r"$\sigma_{\min}(W)^2$")
    axes[1].set_title("Exact cost survives until collision")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].text(0.004, 0.18, r"$d_{\min}=4$ for every $\epsilon>0$", fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--scalar-trials", type=int, default=512)
    parser.add_argument("--block-trials", type=int, default=256)
    args = parser.parse_args()
    certificate = make_certificate(args.scalar_trials, args.block_trials)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_figure(certificate, args.figure)
    if not certificate["all_pass"]:
        raise SystemExit("unified routing-table certificate failed")
    print(json.dumps(certificate["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
