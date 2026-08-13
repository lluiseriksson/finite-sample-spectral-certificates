#!/usr/bin/env python3
"""Independent fail-closed verifier for the frozen planar certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results" / "planar_projective_memory" / "certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    claimed = data.pop("content_sha256_without_hash_field")
    canonical = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    actual = hashlib.sha256(canonical).hexdigest()
    if actual != claimed:
        raise RuntimeError(f"certificate hash mismatch: {actual} != {claimed}")
    if data["schema"] != "planar-projective-memory-certificate-v3":
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
    approximation = data["quantitative_approximation_gap"]
    I = sp.I
    nodes = [sp.Integer(1), I, sp.Integer(-1), -I]
    targets = [(0, 1), (1, 1), (1, 0), (2, 1)]
    expected_error_squared = {0: sp.Rational(1, 5), 1: sp.Rational(3, 67)}
    for item in approximation["certified_degree_caps"]:
        d = item["degree_cap"]
        size = 2 * (d + 1)
        gram = sp.zeros(size, size)
        for node, (a, b) in zip(nodes, targets, strict=True):
            row = sp.Matrix([[b * node**j for j in range(d + 1)] + [-a * node**j for j in range(d + 1)]])
            gram += sp.conjugate(row).T * row / (a * a + b * b)
        recorded = sp.Matrix([[sp.sympify(value) for value in row] for row in item["gram"]])
        if recorded != gram.applyfunc(sp.simplify):
            raise RuntimeError(f"quantitative Gram mismatch at degree {d}")
        alpha = sp.sympify(item["certified_lambda_lower"])
        shifted = gram - alpha * sp.eye(size)
        minors = [sp.factor(shifted[:k, :k].det()) for k in range(1, size + 1)]
        if not all(value > 0 for value in minors):
            raise RuntimeError(f"Sylvester coercivity failure at degree {d}")
        if [str(value) for value in minors] != item["shifted_leading_principal_minors"]:
            raise RuntimeError(f"principal-minor record mismatch at degree {d}")
        error_squared = sp.factor(alpha / (len(nodes) * (d + 1) + alpha))
        if error_squared != expected_error_squared[d]:
            raise RuntimeError(f"unexpected quantitative lower bound at degree {d}")
        if item["certified_uniform_chordal_error_squared_lower"] != str(error_squared):
            raise RuntimeError(f"quantitative error record mismatch at degree {d}")
    zero_memory = data["exact_zero_memory_family"]
    expected_t = [sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 2), sp.Integer(1)]
    if len(zero_memory["fixtures"]) != len(expected_t):
        raise RuntimeError("unexpected zero-memory fixture coverage")
    for item, t in zip(zero_memory["fixtures"], expected_t, strict=True):
        if sp.sympify(item["t"]) != t:
            raise RuntimeError("zero-memory parameter mismatch")
        overlap = sp.factor((1 - t**2) / (1 + t**2))
        exact_error_squared = sp.factor((1 - overlap) / 2)
        sigma_squared = sp.factor(1 - overlap)
        theorem_bound_squared = sp.factor(sigma_squared / (2 + sigma_squared))
        ratio_squared = sp.factor(theorem_bound_squared / exact_error_squared)
        expected_values = {
            "target_overlap": overlap,
            "exact_zero_memory_error_squared": exact_error_squared,
            "interpolation_sigma_squared": sigma_squared,
            "theorem_3_2_bound_squared": theorem_bound_squared,
            "bound_to_exact_ratio_squared": ratio_squared,
        }
        for key, expected_value in expected_values.items():
            if sp.sympify(item[key]) != expected_value:
                raise RuntimeError(f"zero-memory mismatch at t={t}, field={key}")
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
    print(json.dumps({"status": "PASS", "certificate": str(CERT), "sha256": actual, "phase_cases": len(phase), "generic_cases": len(rows), "binary_cases": len(binary), "quantitative_degree_caps": len(approximation["certified_degree_caps"]), "exact_zero_memory_fixtures": len(zero_memory["fixtures"])}, sort_keys=True))


if __name__ == "__main__":
    main()
