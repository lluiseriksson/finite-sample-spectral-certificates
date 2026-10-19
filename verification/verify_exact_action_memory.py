#!/usr/bin/env python3
"""Independent verifier for the exact action--memory artifact.

This file deliberately does not import the producer in ``research``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "exact_action_memory" / "certificate.json"
TOL = 5.0e-8


def feasible(beta: np.ndarray, action: float, budget: int) -> bool:
    positive = beta[beta > 0.0]
    rank = int(positive.size)
    geometric = float(positive.sum())
    if rank == 0:
        return action == 0.0 or (budget >= 1 and 0.0 < action < 2.0 * math.pi * budget)
    return bool(
        budget >= rank
        and 2.0 * geometric - TOL <= action
        and action <= 2.0 * math.pi * budget - 2.0 * geometric + TOL
    )


def deterministic_faces() -> None:
    fixtures = [
        (np.array([math.pi / 6.0]), 1, [math.pi / 3.0, math.pi, 5.0 * math.pi / 3.0]),
        (np.array([math.pi / 3.0, math.pi / 6.0, 0.0]), 2, [math.pi, 2.0 * math.pi, 3.0 * math.pi]),
        (np.array([math.pi / 3.0, math.pi / 6.0, 0.0]), 3, [math.pi, 3.0 * math.pi, 5.0 * math.pi]),
        (np.array([math.pi / 2.0]), 1, [math.pi]),
        (np.array([math.pi / 2.0, math.pi / 2.0]), 2, [2.0 * math.pi]),
    ]
    for beta, degree, actions in fixtures:
        for action in actions:
            if not feasible(beta, action, degree):
                raise RuntimeError("a closed face of the diamond was rejected")
    if feasible(np.array([math.pi / 5.0, math.pi / 7.0, 0.0]), 2.0, 1):
        raise RuntimeError("rank obstruction was not enforced")
    if feasible(np.array([math.pi / 6.0]), math.pi / 3.0 - 1.0e-4, 1):
        raise RuntimeError("lower action face leaked")
    if feasible(np.array([math.pi / 6.0]), 5.0 * math.pi / 3.0 + 1.0e-4, 1):
        raise RuntimeError("upper action face leaked")
    if not feasible(np.array([0.0]), 0.0, 0):
        raise RuntimeError("constant case was rejected")
    if not feasible(np.array([0.0]), 2.0 * math.pi - 1.0e-6, 1):
        raise RuntimeError("open trivial-subspace region was truncated")
    if feasible(np.array([0.0]), 2.0 * math.pi, 1):
        raise RuntimeError("open upper endpoint was incorrectly included")


def independent_scalar_phase(seed: int = 26081041, trials: int = 512) -> None:
    rng = np.random.default_rng(seed)
    gamma = 0.73
    for _ in range(trials):
        delta = float(rng.uniform(2.0e-3, 2.0 * math.pi - 2.0e-3))
        tangent_delta = math.tan(delta / 4.0)
        tangent_arc = math.tan(gamma / 2.0)
        rho = (tangent_delta - tangent_arc) / (tangent_delta + tangent_arc)
        if not -1.0 < rho < 1.0:
            raise RuntimeError("inverse phase map left the disk")

        def poisson(theta: float) -> float:
            return (1.0 - rho * rho) / (1.0 - 2.0 * rho * math.cos(theta) + rho * rho)

        integral, _ = quad(poisson, -gamma, gamma, epsabs=2.0e-11, epsrel=2.0e-11, limit=300)
        formula = 4.0 * math.atan(((1.0 + rho) / (1.0 - rho)) * tangent_arc)
        if abs(integral - delta) > 2.0e-9 or abs(formula - delta) > 2.0e-11:
            raise RuntimeError("scalar Blaschke phase oracle failed")


def independent_gate_lemma(seed: int = 26081042, trials: int = 2048) -> None:
    rng = np.random.default_rng(seed)
    e1 = np.array([1.0, 0.0], dtype=complex)
    for _ in range(trials):
        alpha = float(rng.uniform(1.0e-5, 0.5 * math.pi))
        delta = float(rng.uniform(2.0 * alpha, 2.0 * math.pi - 2.0 * alpha))
        ratio = math.sin(alpha) / abs(math.sin(delta / 2.0))
        t = 0.5 * (1.0 - math.sqrt(max(0.0, 1.0 - ratio * ratio)))
        lam = (1.0 + (np.exp(1j * delta) - 1.0) * t) / math.cos(alpha)
        v1 = math.sqrt(t)
        v2 = lam * math.sin(alpha) / ((np.exp(1j * delta) - 1.0) * v1)
        v = np.array([v1, v2], dtype=complex)
        gate = np.eye(2) + (np.exp(1j * delta) - 1.0) * np.outer(v, v.conj())
        target = lam * np.array([math.cos(alpha), math.sin(alpha)], dtype=complex)
        if abs(np.linalg.norm(v) - 1.0) > TOL or np.linalg.norm(gate @ e1 - target) > TOL:
            raise RuntimeError("rank-one gate construction failed")
        observed = math.acos(float(np.clip(abs((gate @ e1)[0]), 0.0, 1.0)))
        if abs(observed - alpha) > TOL:
            raise RuntimeError("gate moved the line by the wrong angle")


def independent_coincident_compiler() -> None:
    gamma = 0.73
    for action in [1.0e-8, 2.0 * math.pi, 2.0 * math.pi + 0.1, 4.0 * math.pi - 1.0e-4]:
        degree = math.floor(action / (2.0 * math.pi)) + 1
        phases = np.full(degree, action / degree)
        if not np.all((phases > 0.0) & (phases < 2.0 * math.pi)):
            raise RuntimeError("coincident-subspace split hit a forbidden scalar endpoint")
        recovered = 0.0
        for phase in phases:
            tangent_phase = math.tan(float(phase) / 4.0)
            tangent_arc = math.tan(gamma / 2.0)
            rho = (tangent_phase - tangent_arc) / (tangent_phase + tangent_arc)
            recovered += 4.0 * math.atan(((1.0 + rho) / (1.0 - rho)) * tangent_arc)
        if abs(recovered - action) > 2.0e-10:
            raise RuntimeError("coincident-subspace phase split failed")
        if degree != math.floor(action / (2.0 * math.pi)) + 1:
            raise RuntimeError("coincident-subspace compiler used the wrong degree")


def exact_holonomy() -> None:
    omega = np.exp(2j * math.pi / 3.0)
    nodes = np.array([1.0 + 0j, omega, omega**2])
    pick = np.eye(3, dtype=complex)
    for i in range(3):
        for j in range(3):
            if i != j:
                pick[i, j] = 1.0 / (1.0 - np.conj(nodes[i]) * nodes[j])
    eigenvalues = np.linalg.eigvalsh(pick)
    cycle = pick[0, 1] * pick[1, 2] * pick[2, 0]
    if not np.allclose(eigenvalues, [0.0, 1.0, 2.0], atol=TOL, rtol=TOL):
        raise RuntimeError("three-node Pick spectrum changed")
    if abs(cycle - 1j / (3.0 * math.sqrt(3.0))) > TOL:
        raise RuntimeError("three-node cycle holonomy changed")
    for i in range(3):
        for j in range(i + 1, 3):
            cross = pick[i, j]
            pair = np.array([[abs(cross), cross], [np.conj(cross), abs(cross)]])
            pair_eigenvalues = np.linalg.eigvalsh(pair)
            if pair_eigenvalues[0] < -TOL or np.count_nonzero(pair_eigenvalues > TOL) != 1:
                raise RuntimeError("a pairwise subproblem is not rank one")


def verify_frozen_artifact() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or not payload.get("all_pass"):
        raise RuntimeError("frozen certificate is absent or not passing")
    if int(payload["random_diamond"]["trials"]) < 120:
        raise RuntimeError("random diamond campaign is too small")
    if int(payload["noisy_holonomy"]["trials"]) < 2000:
        raise RuntimeError("holonomy noise campaign is too small")
    flagship = payload["flagship_diamond"]
    beta = np.asarray(flagship["positive_angles"])
    action = float(flagship["action"])
    degree = int(flagship["degree"])
    if not feasible(beta, action, degree):
        raise RuntimeError("flagship lies outside its claimed diamond")
    if abs(sum(flagship["gate_phases"]) - action) > TOL:
        raise RuntimeError("gate phases do not add to the frozen action")
    if abs(float(flagship["circle_action_numeric"]) - 2.0 * math.pi * degree) > TOL:
        raise RuntimeError("full-circle winding budget is incorrect")
    faces = payload["deterministic_faces"]
    if not faces.get("all_pass") or len(faces.get("cases", {})) < 12:
        raise RuntimeError("deterministic boundary campaign is incomplete")
    if int(faces["cases"]["coincident_positive_action"]["degree"]) != 1:
        raise RuntimeError("coincident-subspace positive action used the wrong degree")
    for name in ("coincident_one_turn", "coincident_above_turn", "coincident_near_two_turns"):
        if int(faces["cases"][name]["degree"]) != 2 or not faces["cases"][name]["all_pass"]:
            raise RuntimeError("multi-turn coincident-subspace compiler failed")
    noncommuting = payload["noncommuting_order_control"]
    if noncommuting["projector_commutator_norm"] <= 0.49:
        raise RuntimeError("noncommuting order fixture accidentally commutes")
    if noncommuting["reversed_endpoint_projector_error"] <= 0.1:
        raise RuntimeError("reversed order did not trigger the negative control")
    for case in payload["random_diamond"]["sample"]:
        if int(case["degree"]) > int(case["requested_budget"]) or int(case["budget_slack"]) < 0:
            raise RuntimeError("constructed degree exceeded the requested budget")
    holonomy = payload["three_node_holonomy"]
    if holonomy["global_minimum_degree"] != 2 or holonomy["pairwise_minimum_degree"] != 1:
        raise RuntimeError("frozen holonomy gap changed")
    if float(holonomy["maximum_interpolation_error"]) > TOL:
        raise RuntimeError("lurking-isometry interpolant missed a target")
    if float(holonomy["colligation_unitarity_error"]) > TOL:
        raise RuntimeError("synthesized colligation is not unitary")
    if holonomy["controllability_rank"] != 2 or holonomy["observability_rank"] != 2:
        raise RuntimeError("synthesized degree-two realization is not minimal")
    noisy = payload["noisy_holonomy"]
    if int(noisy["rank_one_null_trials"]) < 1000 or int(noisy["false_certificates"]) != 0:
        raise RuntimeError("noisy cycle witness failed its rank-one null controls")
    if float(flagship["angle_error"]) >= 3.0e-16:
        raise RuntimeError("published flagship angle error changed")
    if float(flagship["arc_action_error"]) >= 2.0e-14:
        raise RuntimeError("published flagship arc-action error changed")
    if float(flagship["circle_action_error"]) >= 2.0e-14:
        raise RuntimeError("published flagship winding error changed")
    if float(payload["random_diamond"]["worst_arc_action_error"]) >= 1.0e-10:
        raise RuntimeError("published random-campaign action error changed")


def main() -> None:
    deterministic_faces()
    independent_scalar_phase()
    independent_gate_lemma()
    independent_coincident_compiler()
    exact_holonomy()
    verify_frozen_artifact()
    print("exact action-memory verifier: PASS")


if __name__ == "__main__":
    main()
