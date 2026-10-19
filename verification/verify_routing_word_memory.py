#!/usr/bin/env python3
"""Independent structural checks for the routing-word memory artifact."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "routing_word_memory" / "certificate.json"


def degree_law(word: str, rank: int = 1) -> int:
    counts = Counter(word)
    return rank * (len(word) - min(counts.values()))


def node_polynomials(word: str, nodes: np.ndarray) -> tuple[list[str], list[np.ndarray]]:
    alphabet = sorted(set(word))
    output = []
    for symbol in alphabet:
        polynomial = np.array([1.0 + 0.0j])
        for label, node in zip(word, nodes, strict=True):
            if label != symbol:
                polynomial = np.polynomial.polynomial.polymul(
                    polynomial, np.array([-node, 1.0 + 0.0j])
                )
        output.append(polynomial)
    return alphabet, output


def main() -> None:
    expected = {
        "AAAA": 0,
        "AAAB": 3,
        "AABB": 2,
        "ABAB": 2,
        "AABC": 3,
        "ABCD": 3,
        "AABBCC": 4,
    }
    for word, degree in expected.items():
        if degree_law(word) != degree:
            raise RuntimeError(f"degree law failed for {word}")
        if degree_law(word, 3) != 3 * degree:
            raise RuntimeError(f"rank scaling failed for {word}")

    # The compiler polynomials have exactly the advertised zero pattern.
    word = "AABCAB"
    nodes = np.exp(1j * np.array([-2.4, -1.1, -0.2, 0.7, 1.9, 2.7]))
    alphabet, polynomials = node_polynomials(word, nodes)
    for symbol, polynomial in zip(alphabet, polynomials, strict=True):
        for label, node in zip(word, nodes, strict=True):
            value = np.polynomial.polynomial.polyval(node, polynomial)
            if label != symbol and abs(value) > 2.0e-12:
                raise RuntimeError("wrong-target node polynomial did not vanish")
            if label == symbol and abs(value) < 1.0e-8:
                raise RuntimeError("desired-target node polynomial vanished")
        if len(polynomial) - 1 != len(word) - Counter(word)[symbol]:
            raise RuntimeError("node-polynomial degree mismatch")

    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if certificate["schema"] != "routing-word-memory-certificate-v1":
        raise RuntimeError("certificate schema mismatch")
    if certificate["theorem"] != "d_min = k * (L - min_a n_a)":
        raise RuntimeError("theorem string mismatch")
    if not certificate["all_pass"]:
        raise RuntimeError("frozen certificate is not passing")
    body = dict(certificate)
    digest = body.pop("content_digest")
    body.pop("all_pass")
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    if hashlib.sha256(canonical).hexdigest() != digest:
        raise RuntimeError("certificate content digest mismatch")
    if certificate["permutation_invariance"]["failures"] != 0:
        raise RuntimeError("permutation-invariance campaign failed")
    examples = certificate["separation_examples"]
    if not (
        examples["AAAB"]["cyclic_switches"] < examples["ABAB"]["cyclic_switches"]
        and examples["AAAB"]["scalar_degree"] > examples["ABAB"]["scalar_degree"]
    ):
        raise RuntimeError("switch-count separation missing")
    print("routing-word memory verifier: PASS")


if __name__ == "__main__":
    main()
