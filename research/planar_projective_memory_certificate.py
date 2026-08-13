#!/usr/bin/env python3
"""Exact certificates for planar projective passive-memory laws.

All algebraic checks use SymPy exact arithmetic over Gaussian rationals.
Floating point is used only for the explanatory figure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import sympy as sp

I = sp.I
z = sp.Symbol("z")


def circle_node(t: int) -> sp.Expr:
    """A Gaussian-rational point on the unit circle."""
    return sp.cancel((1 + I * sp.Rational(t)) / (1 - I * sp.Rational(t)))


def cross_ratio(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    """Ordered cross-ratio (a,b;c,d), with a small infinity chart."""
    oo = sp.oo
    if c == oo:
        return sp.cancel((b - d) / (a - d))
    if d == oo:
        return sp.cancel((a - c) / (b - c))
    if a == oo:
        return sp.cancel((b - d) / (b - c))
    if b == oo:
        return sp.cancel((a - c) / (a - d))
    return sp.cancel((a - c) * (b - d) / ((a - d) * (b - c)))


def projective_value(p: sp.Expr, q: sp.Expr, node: sp.Expr) -> sp.Expr:
    pn = sp.simplify(p.subs(z, node))
    qn = sp.simplify(q.subs(z, node))
    if qn == 0:
        if pn == 0:
            raise AssertionError("base point at interpolation node")
        return sp.oo
    return sp.cancel(pn / qn)


def interpolation_matrix(nodes: list[sp.Expr], targets: list[sp.Expr], d: int) -> sp.Matrix:
    """M_d for affine targets [a:b]=[y:1], with infinity [1:0]."""
    rows = []
    for node, target in zip(nodes, targets, strict=True):
        if target == sp.oo:
            a, b = sp.Integer(1), sp.Integer(0)
        else:
            a, b = target, sp.Integer(1)
        rows.append([b * node**j for j in range(d + 1)] + [-a * node**j for j in range(d + 1)])
    return sp.Matrix(rows)


def poly_degree(expr: sp.Expr) -> int:
    poly = sp.Poly(sp.expand(expr), z, extension=I)
    return 0 if poly.is_zero else int(poly.degree())


def certify_witness(nodes: list[sp.Expr], targets: list[sp.Expr], p: sp.Expr, q: sp.Expr, *, include_resultant: bool = True) -> dict:
    values = [projective_value(p, q, node) for node in nodes]
    assert values == targets, (values, targets)
    gcd = sp.gcd(sp.Poly(p, z, extension=I), sp.Poly(q, z, extension=I))
    assert gcd.degree() == 0
    nonconstant_pair = poly_degree(p) > 0 or poly_degree(q) > 0
    resultant = sp.simplify(sp.resultant(p, q, z)) if include_resultant and nonconstant_pair else "not needed (coprimality certified directly)"
    if include_resultant and nonconstant_pair:
        assert resultant != 0
    return {
        "degree": max(poly_degree(p), poly_degree(q)),
        "p": str(sp.expand(p)),
        "q": str(sp.expand(q)),
        "resultant": str(resultant),
        "gcd_degree": gcd.degree(),
        "values": [str(v) for v in values],
    }


def strict_cross_ratio_fixture() -> dict:
    nodes = [sp.Integer(1), I, sp.Integer(-1), -I]
    targets = [sp.Integer(0), sp.Integer(1), sp.oo, sp.Integer(2)]
    p = -((z - 1) * (3 * z - I)) / 2
    q = z + 1
    node_cr = cross_ratio(*nodes)
    target_cr = cross_ratio(*targets)
    assert node_cr == 2
    assert target_cr == sp.Rational(1, 2)
    cert = certify_witness(nodes, targets, p, q)
    assert cert["degree"] == 2
    m1 = interpolation_matrix(nodes, targets, 1)
    assert m1.rank() == 4
    return {
        "nodes": [str(v) for v in nodes],
        "targets": [str(v) for v in targets],
        "node_cross_ratio": str(node_cr),
        "target_cross_ratio": str(target_cr),
        "M1_rank": m1.rank(),
        "M1_shape": list(m1.shape),
        "constant_detector_bound": 1,
        "certificate": cert,
    }


def normalized_projective_gram(nodes: list[sp.Expr], targets: list[sp.Expr], d: int) -> sp.Matrix:
    """Return M_d^* M_d for unit homogeneous target representatives.

    Square roots from normalizing [a:b] cancel in each row outer product, so
    Gaussian-rational data produce an exact Gaussian-rational Gram matrix.
    """
    size = 2 * (d + 1)
    gram = sp.zeros(size, size)
    for node, target in zip(nodes, targets, strict=True):
        if target == sp.oo:
            a, b = sp.Integer(1), sp.Integer(0)
        else:
            a, b = target, sp.Integer(1)
        row = sp.Matrix([[b * node**j for j in range(d + 1)] + [-a * node**j for j in range(d + 1)]])
        norm_squared = sp.simplify(sp.conjugate(a) * a + sp.conjugate(b) * b)
        gram += sp.conjugate(row).T * row / norm_squared
    return gram.applyfunc(sp.simplify)


def quantitative_approximation_fixture() -> dict:
    """Exact coercivity certificates for the strict four-node fixture."""
    nodes = [sp.Integer(1), I, sp.Integer(-1), -I]
    targets = [sp.Integer(0), sp.Integer(1), sp.oo, sp.Integer(2)]
    certified = []
    for d, alpha in ((0, sp.Integer(1)), (1, sp.Rational(3, 8))):
        gram = normalized_projective_gram(nodes, targets, d)
        shifted = gram - alpha * sp.eye(gram.rows)
        leading_minors = [sp.factor(shifted[:k, :k].det()) for k in range(1, shifted.rows + 1)]
        assert all(value > 0 for value in leading_minors)
        # Sylvester's criterion proves M_d^*M_d > alpha I.  The quantitative
        # routing theorem then yields epsilon_infty^2 > alpha/(L(d+1)+alpha).
        error_squared = sp.factor(alpha / (len(nodes) * (d + 1) + alpha))
        certified.append(
            {
                "degree_cap": d,
                "gram": [[str(value) for value in gram.row(k)] for k in range(gram.rows)],
                "gram_trace": str(sp.factor(sp.trace(gram))),
                "gram_determinant": str(sp.factor(gram.det())),
                "certified_lambda_lower": str(alpha),
                "shifted_leading_principal_minors": [str(value) for value in leading_minors],
                "certified_uniform_chordal_error_squared_lower": str(error_squared),
                "certified_uniform_chordal_error_lower": str(sp.sqrt(error_squared)),
            }
        )
    assert certified[0]["certified_uniform_chordal_error_squared_lower"] == "1/5"
    assert certified[1]["certified_uniform_chordal_error_squared_lower"] == "3/67"
    return {
        "metric": "complex-projective chordal distance",
        "target_representatives": ["[0:1]", "[1:1]/sqrt(2)", "[1:0]", "[2:1]/sqrt(5)"],
        "certified_degree_caps": certified,
        "degree_two_exact_error": "0",
    }


def four_line_fixtures() -> dict:
    nodes = [circle_node(t) for t in (0, 1, 2, 3)]
    x1, x2, x3, x4 = nodes

    fixtures: dict[str, dict] = {}
    fixtures["all_same"] = {
        "expected_degree": 0,
        "reason": "constant projective map",
        "certificate": certify_witness(nodes, [sp.Integer(0)] * 4, sp.Integer(0), sp.Integer(1)),
    }
    fixtures["all_distinct_compatible"] = {
        "expected_degree": 1,
        "reason": "ordered cross-ratios agree",
        "certificate": certify_witness(nodes, nodes, z, sp.Integer(1)),
    }
    fixtures["2+2"] = {
        "expected_degree": 2,
        "reason": "two double fibers",
        "certificate": certify_witness(
            nodes,
            [sp.Integer(0), sp.Integer(0), sp.oo, sp.oo],
            (z - x1) * (z - x2),
            (z - x3) * (z - x4),
        ),
    }
    fixtures["3+1"] = {
        "expected_degree": 3,
        "reason": "a nonconstant degree-d fiber has at most d distinct preimages",
        "certificate": certify_witness(
            nodes,
            [sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.oo],
            (z - x1) * (z - x2) * (z - x3),
            z - x4,
        ),
    }
    # Use a separate normalized projective chart for the 2+1+1 construction.
    nn = [sp.Integer(0), sp.oo, sp.Integer(1), sp.Integer(2)]
    # Infinity is represented by chart evaluation of leading coefficients below.
    a = sp.Integer(3)
    A = sp.Rational(2, 1) / ((2 - 1) * (2 - a))
    p = sp.expand(A * (z - 1) * (z - a))
    q = z
    finite_values = [projective_value(p, q, nn[j]) for j in (0, 2, 3)]
    infinity_value = sp.oo if poly_degree(p) > poly_degree(q) else sp.LC(sp.Poly(p, z)) / sp.LC(sp.Poly(q, z))
    values = [finite_values[0], infinity_value, finite_values[1], finite_values[2]]
    assert values == [sp.oo, sp.oo, 0, 1]
    res = sp.simplify(sp.resultant(p, q, z))
    assert res != 0
    fixtures["2+1+1"] = {
        "expected_degree": 2,
        "reason": "explicit normalized quadratic witness",
        "certificate": {"degree": 2, "p": str(p), "q": str(q), "resultant": str(res), "values": [str(v) for v in values]},
    }
    return fixtures


def generic_planar_campaign(max_L: int) -> list[dict]:
    rows = []
    for L in range(3, max_L + 1):
        d0 = (L - 1 + 1) // 2  # ceil((L-1)/2)
        nodes = [circle_node(t) for t in range(L)]
        targets = [sp.cancel(node**d0) for node in nodes]
        lower_rank = None
        if d0 > 0:
            m_lower = interpolation_matrix(nodes, targets, d0 - 1)
            lower_rank = m_lower.rank()
            assert m_lower.nullspace() == []
        m0 = interpolation_matrix(nodes, targets, d0)
        vec = sp.Matrix([0] * d0 + [1] + [1] + [0] * d0)
        # Coefficients are p=z^d0 followed by q=1.
        assert (m0 * vec).applyfunc(sp.simplify) == sp.zeros(L, 1)
        assert sp.gcd(sp.Poly(z**d0, z), sp.Poly(1, z)).degree() == 0
        rows.append(
            {
                "L": L,
                "generic_degree": d0,
                "detector_bound": 1,
                "lower_matrix_rank": lower_rank,
                "lower_matrix_columns": 2 * d0,
                "threshold_matrix_rank": m0.rank(),
                "threshold_matrix_shape": list(m0.shape),
                "exact_witness": "[z^%d:1]" % d0,
                "witness_role": "nonempty rank-and-coprimality open set; target-distinctness is intersected as a separate nonempty Zariski open condition",
            }
        )
    return rows


def binary_collision_campaign(max_L: int) -> list[dict]:
    rows = []
    for L in range(3, max_L + 1):
        nodes = [circle_node(t) for t in range(L)]
        for n_zero in range(1, L):
            n_inf = L - n_zero
            expected = max(n_zero, n_inf)
            # The two displayed factor lists are disjoint because the nodes are
            # distinct.  They are therefore a coprime exact projective witness
            # without requiring costly expanded resultants.
            cert = {
                "degree": expected,
                "p_roots": [str(node) for node in nodes[:n_zero]],
                "q_roots": [str(node) for node in nodes[n_zero:]],
                "coprime": True,
                "construction": "p=product(z-p_root), q=product(z-q_root)",
            }
            rows.append({"L": L, "occupancies": [n_zero, n_inf], "expected_degree": expected, "certificate": cert})
    return rows


def make_figure(rows: list[dict], approximation: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    xs = [r["L"] for r in rows]
    ds = [r["generic_degree"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.25))
    ax = axes[0]
    ax.step(xs, ds, where="mid", linewidth=2.2, label=r"generic $d_{\min}=\lceil(L-1)/2\rceil$")
    ax.plot(xs, [1] * len(xs), "--", linewidth=2, label="constant-detector bound")
    ax.fill_between(xs, 1, ds, step="mid", alpha=0.16, color="#2563eb")
    ax.set_xlabel("number of distinct coplanar target lines $L$")
    ax.set_ylabel("state count and bound values")
    ax.set_xticks(xs[::2] if len(xs) > 8 else xs)
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, loc="upper left")
    ax = axes[1]
    error_bounds = [
        float(sp.sqrt(sp.Rational(item["certified_uniform_chordal_error_squared_lower"])))
        for item in approximation["certified_degree_caps"]
    ] + [0.0]
    ax.bar([0, 1, 2], error_bounds, color=["#7c3aed", "#2563eb", "#94a3b8"], width=0.68)
    ax.set_xticks([0, 1, 2])
    ax.set_xlabel("router degree cap $d$")
    ax.set_ylabel("certified worst-node error")
    ax.set_title("strict four-node fixture")
    ax.set_ylim(0, 0.5)
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--max-L", type=int, default=16)
    args = parser.parse_args()

    payload = {
        "schema": "planar-projective-memory-certificate-v2",
        "arithmetic": "SymPy exact Gaussian-rational arithmetic",
        "strict_cross_ratio_gap": strict_cross_ratio_fixture(),
        "quantitative_approximation_gap": quantitative_approximation_fixture(),
        "four_line_phase_fixtures": four_line_fixtures(),
        "generic_planar_campaign": generic_planar_campaign(args.max_L),
        "binary_collision_campaign": binary_collision_campaign(args.max_L),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["content_sha256_without_hash_field"] = hashlib.sha256(canonical).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_figure(payload["generic_planar_campaign"], payload["quantitative_approximation_gap"], args.figure)
    print(json.dumps({"output": str(args.output), "figure": str(args.figure), "cases": len(payload["generic_planar_campaign"]), "sha256": payload["content_sha256_without_hash_field"]}, sort_keys=True))


if __name__ == "__main__":
    main()
