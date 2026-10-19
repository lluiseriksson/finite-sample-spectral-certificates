#!/usr/bin/env python3
"""Independent fail-closed verifier for the unified routing-table artifact."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "unified_routing_table_memory" / "certificate.json"


def degree_law(word: str, rank: int = 1) -> int:
    counts = Counter(word)
    return rank * (len(word) - min(counts.values()))


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
    for word, value in expected.items():
        if degree_law(word) != value or degree_law(word, 3) != 3 * value:
            raise RuntimeError(f"degree law mismatch for {word}")

    # Independent direct-sum dual-detector check.
    rng = np.random.default_rng(99173)
    blocks = []
    for _ in range(4):
        raw = rng.normal(size=(11, 2)) + 1j * rng.normal(size=(11, 2))
        q, _ = np.linalg.qr(raw)
        blocks.append(q[:, :2])
    stack = np.hstack(blocks)
    if np.linalg.matrix_rank(stack) != 8:
        raise RuntimeError("independent frame fixture lost direct-sum rank")
    dual = np.linalg.inv(stack.conj().T @ stack) @ stack.conj().T
    if np.linalg.norm(dual @ stack - np.eye(8), 2) > 1e-11:
        raise RuntimeError("dual block detectors failed")

    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if certificate["schema"] != "unified-routing-table-memory-certificate-v1":
        raise RuntimeError("certificate schema mismatch")
    if certificate["theorem"] != "d_min = k * (L - min_a n_a) on the direct-sum target locus":
        raise RuntimeError("certificate theorem mismatch")
    if certificate["curated_degree_examples"] != expected:
        raise RuntimeError("curated degree table mismatch")
    if not certificate["all_pass"]:
        raise RuntimeError("frozen certificate is not passing")
    body = dict(certificate)
    digest = body.pop("content_digest")
    body.pop("all_pass")
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    if hashlib.sha256(canonical).hexdigest() != digest:
        raise RuntimeError("certificate digest mismatch")

    summary = certificate["summary"]
    thresholds = {
        "maximum_dual_isolation_residual": 5e-8,
        "maximum_node_range_projector_error": 2e-8,
        "maximum_wrong_dual_minor": 2e-8,
        "maximum_scalar_factor_residual": 5e-6,
        "maximum_scalar_inner_residual": 5e-5,
        "maximum_scalar_node_projector_error": 5e-8,
    }
    for key, threshold in thresholds.items():
        if not summary[key] < threshold:
            raise RuntimeError(f"{key} exceeds threshold")
    if not summary["minimum_density_eigenvalue"] > 1e-10:
        raise RuntimeError("matrix density lost strict positivity")
    if len(certificate["collision_family"]) != 6:
        raise RuntimeError("collision sweep length changed")
    if any(case["predicted_degree"] != 4 for case in certificate["collision_family"]):
        raise RuntimeError("positive-opening collision law changed")
    print("unified routing-table verifier: PASS")


if __name__ == "__main__":
    main()
