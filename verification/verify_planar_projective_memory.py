#!/usr/bin/env python3
"""Independent fail-closed verifier for the frozen planar certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results" / "planar_projective_memory" / "certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    claimed = data.pop("content_sha256_without_hash_field")
    canonical = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    actual = hashlib.sha256(canonical).hexdigest()
    if actual != claimed:
        raise RuntimeError(f"certificate hash mismatch: {actual} != {claimed}")
    if data["schema"] != "planar-projective-memory-certificate-v1":
        raise RuntimeError("unexpected schema")
    strict = data["strict_cross_ratio_gap"]
    expected = {
        "node_cross_ratio": "2",
        "target_cross_ratio": "1/2",
        "M1_rank": 4,
        "M1_shape": [4, 4],
        "constant_detector_bound": 1,
    }
    for key, value in expected.items():
        if strict[key] != value:
            raise RuntimeError(f"strict fixture mismatch at {key}: {strict[key]} != {value}")
    if strict["certificate"]["degree"] != 2:
        raise RuntimeError("strict fixture must have exact degree two")
    phase = data["four_line_phase_fixtures"]
    expected_phase = {"all_same": 0, "all_distinct_compatible": 1, "2+2": 2, "2+1+1": 2, "3+1": 3}
    for key, degree in expected_phase.items():
        if phase[key]["expected_degree"] != degree or phase[key]["certificate"]["degree"] != degree:
            raise RuntimeError(f"phase fixture mismatch: {key}")
    rows = data["generic_planar_campaign"]
    if [r["L"] for r in rows] != list(range(3, 17)):
        raise RuntimeError("generic campaign does not cover L=3,...,16")
    for row in rows:
        expected_degree = (row["L"] - 1 + 1) // 2
        if row["generic_degree"] != expected_degree:
            raise RuntimeError(f"wrong generic degree at L={row['L']}")
        if row["detector_bound"] != 1:
            raise RuntimeError(f"wrong detector bound at L={row['L']}")
        if row["lower_matrix_rank"] != row["lower_matrix_columns"]:
            raise RuntimeError(f"lower matrix is not full-column-rank at L={row['L']}")
    binary = data["binary_collision_campaign"]
    if len(binary) != sum(L - 1 for L in range(3, 17)):
        raise RuntimeError("binary campaign coverage mismatch")
    for row in binary:
        expected = max(row["occupancies"])
        if row["expected_degree"] != expected or row["certificate"]["degree"] != expected:
            raise RuntimeError(f"binary collision mismatch at {row['occupancies']}")
    print(json.dumps({"status": "PASS", "certificate": str(CERT), "sha256": actual, "phase_cases": len(phase), "generic_cases": len(rows), "binary_cases": len(binary)}, sort_keys=True))


if __name__ == "__main__":
    main()
