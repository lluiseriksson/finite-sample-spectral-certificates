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
    if data["schema"] != "planar-projective-memory-certificate-v6":
        raise RuntimeError("unexpected schema")
    law = data["base_point_deletion_law"]
    if law != {
        "border_memory": "min over B subset [L] of |B| + delta(B^c)",
        "linear_rank_formula": "min d with rank(M_d) < 2(d+1)",
        "minimum_modulus_dichotomy": "E_d=0 iff inf_{||c||=1} ||M_d c||=0",
        "dense_elimination_field_operation_bound": "O(L^4)",
        "three_regimes": ["positive error", "zero unattained infimum", "exact realization"],
        "binary_specialization": "border=min(n0,n_infinity), exact=max(n0,n_infinity)",
    }:
        raise RuntimeError("base-point deletion law metadata mismatch")
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
    closure = data["universal_closure_family"]
    if closure["L"] != 3 or closure["threshold_degree"] != 1 or closure["exact_interpolation_degree"] != 2:
        raise RuntimeError("unexpected universal-closure threshold fixture")
    expected_epsilons = [sp.Rational(1, 2), sp.Rational(1, 10), sp.Rational(1, 100), sp.Rational(1, 1000)]
    if len(closure["degenerating_degree_one_sequence"]) != len(expected_epsilons):
        raise RuntimeError("unexpected universal-closure sequence coverage")
    z = sp.Symbol("z")
    for item, epsilon in zip(closure["degenerating_degree_one_sequence"], expected_epsilons, strict=True):
        p = sp.expand(epsilon * z)
        q = sp.expand(epsilon * z + 1 - epsilon)
        determinant = sp.factor(epsilon * (1 - epsilon))
        resultant = sp.factor(sp.resultant(p, q, z))
        error_squared = sp.factor(epsilon**2 / (1 + epsilon**2))
        expected_values = {
            "epsilon": epsilon,
            "p": p,
            "q": q,
            "mobius_determinant": determinant,
            "resultant": resultant,
            "worst_node_chordal_error_squared": error_squared,
        }
        for key, expected_value in expected_values.items():
            if sp.sympify(item[key]) != expected_value:
                raise RuntimeError(f"universal-closure mismatch at epsilon={epsilon}, field={key}")
        if item["values_at_0_1_infinity"] != ["0", str(epsilon), "1"]:
            raise RuntimeError(f"universal-closure values mismatch at epsilon={epsilon}")
    exact_witness = closure["exact_degree_two_witness"]
    exact_p = sp.expand(z * (z - 1))
    exact_q = sp.expand(z * (z - 1) + 1)
    if sp.sympify(exact_witness["p"]) != exact_p or sp.sympify(exact_witness["q"]) != exact_q:
        raise RuntimeError("universal-closure exact witness mismatch")
    if sp.sympify(exact_witness["resultant"]) != sp.resultant(exact_p, exact_q, z):
        raise RuntimeError("universal-closure exact resultant mismatch")
    phase = data["four_line_phase_fixtures"]
    expected_phase = {"all_same": 0, "all_distinct_compatible": 1, "2+2": 2, "2+1+1": 2, "3+1": 3}
    expected_border = {"all_same": 0, "all_distinct_compatible": 1, "2+2": 2, "2+1+1": 2, "3+1": 1}
    phase_nodes = [sp.cancel((1 + sp.I * sp.Rational(t)) / (1 - sp.I * sp.Rational(t))) for t in range(4)]
    phase_targets = {
        "all_same": [(0, 1)] * 4,
        "all_distinct_compatible": [(node, 1) for node in phase_nodes],
        "2+2": [(0, 1), (0, 1), (1, 0), (1, 0)],
        "2+1+1": [(1, 0), (1, 0), (0, 1), (1, 1)],
        "3+1": [(0, 1), (0, 1), (0, 1), (1, 0)],
    }
    def phase_matrix(label: str, degree: int) -> sp.Matrix:
        return sp.Matrix([
            [b * node**j for j in range(degree + 1)] + [-a * node**j for j in range(degree + 1)]
            for node, (a, b) in zip(phase_nodes, phase_targets[label], strict=True)
        ])
    for key, degree in expected_phase.items():
        if phase[key]["expected_degree"] != degree or phase[key]["certificate"]["degree"] != degree:
            raise RuntimeError(f"phase fixture mismatch: {key}")
        border = expected_border[key]
        if phase[key]["expected_border_degree"] != border:
            raise RuntimeError(f"phase border mismatch: {key}")
        if phase[key]["positive_error_degrees"] != list(range(border)):
            raise RuntimeError(f"phase positive-error mismatch: {key}")
        if phase[key]["zero_unattained_degrees"] != list(range(border, degree)):
            raise RuntimeError(f"phase nonattainment mismatch: {key}")
        at_border = phase_matrix(key, border)
        border_rank = at_border.rank()
        expected_rank = {
            "border_shape": list(at_border.shape),
            "border_rank": border_rank,
            "border_nullity": at_border.cols - border_rank,
        }
        if border > 0:
            before = phase_matrix(key, border - 1)
            expected_rank.update(
                {
                    "before_border_shape": list(before.shape),
                    "before_border_rank": before.rank(),
                }
            )
        if phase[key]["rank_certificate"] != expected_rank:
            raise RuntimeError(f"phase rank certificate mismatch: {key}")
        if expected_rank["border_nullity"] <= 0:
            raise RuntimeError(f"phase border matrix not deficient: {key}")
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
        border = min(row["occupancies"])
        if row["expected_degree"] != expected or row["certificate"]["degree"] != expected:
            raise RuntimeError(f"binary collision mismatch at {row['occupancies']}")
        if row["expected_border_degree"] != border:
            raise RuntimeError(f"binary border mismatch at {row['occupancies']}")
        if row["positive_error_degrees"] != list(range(border)):
            raise RuntimeError(f"binary positive-error phase mismatch at {row['occupancies']}")
        if row["zero_unattained_degrees"] != list(range(border, expected)):
            raise RuntimeError(f"binary nonattainment phase mismatch at {row['occupancies']}")
        if row["exact_zero_error_from_degree"] != expected:
            raise RuntimeError(f"binary exact phase mismatch at {row['occupancies']}")
        n_zero, n_inf = row["occupancies"]
        circle_nodes = [sp.cancel((1 + sp.I * sp.Rational(t)) / (1 - sp.I * sp.Rational(t))) for t in range(row["L"])]
        binary_targets = [(0, 1)] * n_zero + [(1, 0)] * n_inf
        def binary_matrix(degree: int) -> sp.Matrix:
            return sp.Matrix([
                [b * node**j for j in range(degree + 1)] + [-a * node**j for j in range(degree + 1)]
                for node, (a, b) in zip(circle_nodes, binary_targets, strict=True)
            ])
        before = binary_matrix(border - 1)
        at_border = binary_matrix(border)
        before_rank = before.rank()
        border_rank = at_border.rank()
        rank_cert = row["rank_certificate"]
        expected_rank_cert = {
            "before_border_degree": border - 1,
            "before_border_shape": list(before.shape),
            "before_border_rank": before_rank,
            "border_shape": list(at_border.shape),
            "border_rank": border_rank,
            "border_nullity": at_border.cols - border_rank,
        }
        if rank_cert != expected_rank_cert:
            raise RuntimeError(f"binary rank certificate mismatch at {row['occupancies']}")
        if before_rank != before.cols or border_rank >= at_border.cols:
            raise RuntimeError(f"binary first-deficiency law failed at {row['occupancies']}")

    border_gap = data["one_vs_rest_border_gap_campaign"]
    if [row["L"] for row in border_gap] != list(range(3, 17)):
        raise RuntimeError("one-vs-rest border campaign coverage mismatch")
    epsilons = [sp.Rational(1, 2), sp.Rational(1, 10), sp.Rational(1, 100)]
    for row in border_gap:
        L = row["L"]
        if row["occupancies"] != [L - 1, 1] or row["exact_degree"] != L - 1:
            raise RuntimeError(f"one-vs-rest exact law mismatch at L={L}")
        if row["border_degree"] != 1 or row["exact_to_border_ratio"] != L - 1:
            raise RuntimeError(f"one-vs-rest border law mismatch at L={L}")
        circle_nodes = [sp.cancel((1 + sp.I * sp.Rational(t)) / (1 - sp.I * sp.Rational(t))) for t in range(L)]
        exceptional = circle_nodes[-1]
        if len(row["degenerating_degree_one_sequence"]) != len(epsilons):
            raise RuntimeError(f"one-vs-rest sequence coverage mismatch at L={L}")
        for item, epsilon in zip(row["degenerating_degree_one_sequence"], epsilons, strict=True):
            p = epsilon
            q = sp.expand((1 - epsilon) * (z - exceptional))
            if sp.sympify(item["epsilon"]) != epsilon or sp.sympify(item["p"]) != p:
                raise RuntimeError(f"one-vs-rest parameter mismatch at L={L}")
            if sp.simplify(sp.sympify(item["q"]) - q) != 0:
                raise RuntimeError(f"one-vs-rest denominator mismatch at L={L}")
            if sp.simplify(q.subs(z, exceptional)) != 0:
                raise RuntimeError(f"one-vs-rest exceptional value mismatch at L={L}")
            errors = []
            for node in circle_nodes[:-1]:
                value = sp.cancel(p / q.subs(z, node))
                errors.append(sp.factor(value * sp.conjugate(value) / (1 + value * sp.conjugate(value))))
            if sp.sympify(item["worst_node_chordal_error_squared"]) != max(errors):
                raise RuntimeError(f"one-vs-rest error mismatch at L={L}, epsilon={epsilon}")

    print(json.dumps({"status": "PASS", "certificate": str(CERT), "sha256": actual, "phase_cases": len(phase), "generic_cases": len(rows), "binary_cases": len(binary), "independent_border_rank_cases": len(binary), "border_gap_cases": len(border_gap), "quantitative_degree_caps": len(approximation["certified_degree_caps"]), "exact_zero_memory_fixtures": len(zero_memory["fixtures"]), "universal_closure_fixtures": len(closure["degenerating_degree_one_sequence"])}, sort_keys=True))


if __name__ == "__main__":
    main()
