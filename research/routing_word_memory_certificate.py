#!/usr/bin/env python3
"""Deterministic certificate for the rarest-target routing-word law.

For a word w on mutually orthogonal target blocks, the paper proves

    d_min = k * (len(w) - min_a count_w(a)).

The numerical layer independently constructs the scalar inner routing column
used by the upper bound.  It deliberately does not optimize a transfer matrix:
the construction is the proof algorithm (node polynomials, Fejer--Riesz
spectral factor, then lossless completion at the theorem level).
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "routing-word-memory-certificate-v1"
SEED = 260811


def degree_law(word: str, rank: int = 1) -> int:
    counts = Counter(word)
    if not counts:
        raise ValueError("the routing word must be nonempty")
    if rank < 1:
        raise ValueError("rank must be positive")
    return rank * (len(word) - min(counts.values()))


def node_polynomials(word: str, nodes: np.ndarray) -> tuple[list[str], list[np.ndarray]]:
    """Return p_a(z)=prod_{i:w_i != a}(z-zeta_i), ascending coefficients."""
    alphabet = sorted(set(word))
    polynomials: list[np.ndarray] = []
    for symbol in alphabet:
        coefficients = np.array([1.0 + 0.0j])
        for label, node in zip(word, nodes, strict=True):
            if label != symbol:
                coefficients = np.polynomial.polynomial.polymul(
                    coefficients, np.array([-node, 1.0 + 0.0j])
                )
        polynomials.append(np.asarray(coefficients, dtype=complex))
    return alphabet, polynomials


def _pad(polynomial: np.ndarray, degree: int) -> np.ndarray:
    return np.pad(polynomial, (0, degree + 1 - len(polynomial)))


def spectral_factor(polynomials: Iterable[np.ndarray], degree: int) -> tuple[np.ndarray, dict]:
    """Outer Fejer--Riesz factor by reciprocal-root selection.

    If R(e^it)=sum_a |p_a(e^it)|^2, then z^d R(z) is self-inversive.
    Its d exterior roots form a zero-free-in-the-disk polynomial h.  The final
    positive scale is obtained from the exact boundary identity.
    """
    padded = [_pad(np.asarray(p, dtype=complex), degree) for p in polynomials]
    palindromic = np.zeros(2 * degree + 1, dtype=complex)
    for polynomial in padded:
        palindromic += np.convolve(polynomial, np.conj(polynomial[::-1]))
    roots = np.roots(palindromic[::-1]) if degree else np.array([], dtype=complex)
    if degree:
        order = np.argsort(np.abs(roots))[::-1]
        exterior = roots[order[:degree]]
        if np.min(np.abs(exterior)) <= 1.0 - 2.0e-8:
            raise RuntimeError("reciprocal-root split failed")
        h_descending = np.poly(exterior)
        h = h_descending[::-1].astype(complex)
    else:
        exterior = np.array([], dtype=complex)
        h = np.array([1.0 + 0.0j])

    grid = np.exp(1j * np.linspace(-math.pi, math.pi, 8192, endpoint=False))
    r_values = np.zeros(grid.size, dtype=float)
    for polynomial in padded:
        r_values += np.abs(np.polynomial.polynomial.polyval(grid, polynomial)) ** 2
    h_values = np.polynomial.polynomial.polyval(grid, h)
    ratios = r_values / np.abs(h_values) ** 2
    scale = math.sqrt(float(np.median(ratios)))
    h *= scale
    h_values *= scale
    relative_residual = float(
        np.max(np.abs(np.abs(h_values) ** 2 - r_values)) / np.max(r_values)
    )
    return h, {
        "palindromic_residual": float(
            np.max(np.abs(palindromic - np.conj(palindromic[::-1])))
        ),
        "spectral_factor_relative_residual": relative_residual,
        "minimum_exterior_root_modulus": (
            float(np.min(np.abs(exterior))) if degree else None
        ),
        "minimum_boundary_density": float(np.min(r_values)),
        "maximum_boundary_density": float(np.max(r_values)),
    }


def compile_word(word: str, angles: np.ndarray) -> dict:
    nodes = np.exp(1j * np.asarray(angles, dtype=float))
    if len(nodes) != len(word):
        raise ValueError("word and node count disagree")
    if len(nodes) > 1 and min(
        abs(nodes[i] - nodes[j]) for i in range(len(nodes)) for j in range(i)
    ) < 1e-12:
        raise ValueError("nodes must be distinct")
    alphabet, polynomials = node_polynomials(word, nodes)
    counts = Counter(word)
    degree = len(word) - min(counts.values())
    h, diagnostics = spectral_factor(polynomials, degree)
    padded = [_pad(polynomial, degree) for polynomial in polynomials]

    grid = np.exp(1j * np.linspace(-math.pi, math.pi, 4096, endpoint=False))
    h_grid = np.polynomial.polynomial.polyval(grid, h)
    q_grid = np.vstack(
        [np.polynomial.polynomial.polyval(grid, p) / h_grid for p in padded]
    )
    column_isometry_residual = float(np.max(np.abs(np.sum(np.abs(q_grid) ** 2, axis=0) - 1.0)))

    interpolation_error = 0.0
    desired_floor = 1.0
    wrong_leakage = 0.0
    for index, (symbol, node) in enumerate(zip(word, nodes, strict=True)):
        h_node = np.polynomial.polynomial.polyval(node, h)
        q_node = np.array(
            [np.polynomial.polynomial.polyval(node, p) / h_node for p in padded]
        )
        target = alphabet.index(symbol)
        desired_floor = min(desired_floor, float(abs(q_node[target])))
        wrong_leakage = max(
            wrong_leakage,
            float(np.linalg.norm(np.delete(q_node, target))),
        )
        interpolation_error = max(
            interpolation_error,
            float(abs(np.sum(np.abs(q_node) ** 2) - 1.0)),
        )

    return {
        "word": word,
        "length": len(word),
        "alphabet_size": len(alphabet),
        "counts": {symbol: counts[symbol] for symbol in alphabet},
        "rarest_occupancy": min(counts.values()),
        "scalar_degree": degree,
        "rank_k_degree_formula": f"k*{degree}",
        "angles": [float(value) for value in angles],
        "node_separation": float(
            min(abs(nodes[i] - nodes[j]) for i in range(len(nodes)) for j in range(i))
        ) if len(nodes) > 1 else None,
        "column_isometry_residual": column_isometry_residual,
        "node_norm_residual": interpolation_error,
        "wrong_target_leakage": wrong_leakage,
        "minimum_desired_amplitude": desired_floor,
        **diagnostics,
    }


def _random_word(rng: np.random.Generator, length: int, alphabet_size: int) -> str:
    labels = np.arange(alphabet_size)
    tail = rng.integers(0, alphabet_size, size=length - alphabet_size)
    data = np.concatenate([labels, tail])
    rng.shuffle(data)
    return "".join(chr(ord("A") + int(value)) for value in data)


def _canonical_words(length: int) -> list[str]:
    """Restricted-growth strings: one representative modulo alphabet renaming."""
    if length < 1:
        return []
    output: list[str] = []

    def extend(prefix: list[int], maximum: int) -> None:
        if len(prefix) == length:
            output.append("".join(chr(ord("A") + value) for value in prefix))
            return
        for value in range(maximum + 2):
            prefix.append(value)
            extend(prefix, max(maximum, value))
            prefix.pop()

    extend([0], 0)
    return output


def make_certificate(random_trials: int = 1024) -> dict:
    rng = np.random.default_rng(SEED)
    curated_words = ["A", "AB", "ABA", "AAB", "AAAB", "ABAB", "ABC", "ABCA", "AABBCC", "ABCD"]
    cases: list[dict] = []
    for word in curated_words:
        angles = np.sort(rng.uniform(-math.pi, math.pi, size=len(word)))
        cases.append(compile_word(word, angles))

    random_cases: list[dict] = []
    for _ in range(random_trials):
        length = int(rng.integers(2, 10))
        alphabet_size = int(rng.integers(1, min(5, length + 1)))
        word = _random_word(rng, length, alphabet_size)
        angles = np.sort(rng.uniform(-math.pi, math.pi, size=length))
        random_cases.append(compile_word(word, angles))

    exhaustive_cases: list[dict] = []
    for length in range(1, 8):
        for index, word in enumerate(_canonical_words(length)):
            offset = 0.013 * (index + 1) / (len(_canonical_words(length)) + 1)
            angles = np.linspace(-2.7, 2.7, length) + offset
            exhaustive_cases.append(compile_word(word, angles))

    permutation_checks = 256
    permutation_failures = 0
    for _ in range(permutation_checks):
        case = random_cases[int(rng.integers(0, len(random_cases)))]
        shuffled = list(case["word"])
        rng.shuffle(shuffled)
        if degree_law("".join(shuffled)) != degree_law(case["word"]):
            permutation_failures += 1

    stress: list[dict] = []
    for gap in [1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4]:
        angles = np.array([0.0, gap, 1.7, 4.2])
        result = compile_word("ABAB", angles)
        result["prescribed_gap"] = gap
        stress.append(result)

    all_cases = cases + random_cases + exhaustive_cases + stress
    max_isometry = max(case["column_isometry_residual"] for case in all_cases)
    max_leakage = max(case["wrong_target_leakage"] for case in all_cases)
    max_factor = max(case["spectral_factor_relative_residual"] for case in all_cases)
    minimum_desired = min(case["minimum_desired_amplitude"] for case in all_cases)
    payload = {
        "schema": SCHEMA,
        "seed": SEED,
        "theorem": "d_min = k * (L - min_a n_a)",
        "scope": "distinct boundary nodes; square finite rational-inner matrices; regular node values; mutually orthogonal rank-k target planes; every listed symbol occurs",
        "curated_cases": cases,
        "random_summary": {
            "trials": random_trials,
            "maximum_column_isometry_residual": max_isometry,
            "maximum_wrong_target_leakage": max_leakage,
            "maximum_spectral_factor_relative_residual": max_factor,
            "minimum_desired_amplitude": minimum_desired,
        },
        "exhaustive_modulo_renaming": {
            "maximum_length": 7,
            "word_classes": len(exhaustive_cases),
            "maximum_column_isometry_residual": max(
                case["column_isometry_residual"] for case in exhaustive_cases
            ),
            "maximum_wrong_target_leakage": max(
                case["wrong_target_leakage"] for case in exhaustive_cases
            ),
            "maximum_spectral_factor_relative_residual": max(
                case["spectral_factor_relative_residual"] for case in exhaustive_cases
            ),
        },
        "permutation_invariance": {
            "trials": permutation_checks,
            "failures": permutation_failures,
        },
        "clustered_node_stress": stress,
        "separation_examples": {
            "AAAB": {"cyclic_switches": 2, "scalar_degree": degree_law("AAAB")},
            "ABAB": {"cyclic_switches": 4, "scalar_degree": degree_law("ABAB")},
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["content_digest"] = hashlib.sha256(canonical).hexdigest()
    payload["all_pass"] = bool(
        max_isometry < 5.0e-5
        and max_leakage < 5.0e-7
        and max_factor < 5.0e-7
        and minimum_desired > 1.0 - 5.0e-5
        and permutation_failures == 0
    )
    return payload


def make_figure(path: Path) -> None:
    words = ["AAAA", "AAAB", "AABB", "ABAB", "ABCD"]
    degrees = [degree_law(word) for word in words]
    switches = [
        sum(word[i] != word[(i + 1) % len(word)] for i in range(len(word)))
        for word in words
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.7))
    colors = ["#91a4b7", "#d95f59", "#4c78a8", "#59a14f", "#b279a2"]
    axes[0].bar(words, degrees, color=colors)
    axes[0].set_ylabel(r"exact scalar degree $d_{\min}$")
    axes[0].set_title("The rarest-target law at L=4")
    axes[0].set_ylim(0, 3.35)
    for i, value in enumerate(degrees):
        axes[0].text(i, value + 0.08, str(value), ha="center", fontsize=9)
    axes[1].scatter(switches, degrees, c=colors, s=90, edgecolor="black", linewidth=0.5)
    for word, x, y in zip(words, switches, degrees, strict=True):
        axes[1].annotate(word, (x, y), xytext=(5, 5), textcoords="offset points", fontsize=9)
    axes[1].set_xlabel("cyclic switch count")
    axes[1].set_ylabel(r"exact scalar degree $d_{\min}$")
    axes[1].set_title("Switch count does not order memory")
    axes[1].set_xlim(-0.25, 4.65)
    axes[1].set_ylim(-0.15, 3.35)
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--trials", type=int, default=1024)
    args = parser.parse_args()
    certificate = make_certificate(args.trials)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_figure(args.figure)
    if not certificate["all_pass"]:
        raise SystemExit("routing-word memory certificate failed")
    print(json.dumps(certificate["random_summary"], sort_keys=True))


if __name__ == "__main__":
    main()
