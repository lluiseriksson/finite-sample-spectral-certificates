#!/usr/bin/env python3
"""Independent arithmetic checks for the routing speed-limit certificate."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "routing_speed_limit" / "certificate.json"


def close(left: float, right: float, tolerance: float = 2.0e-11) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise RuntimeError("unexpected routing certificate schema")

    transitions = payload["finite_error_local_equality_cases"]
    if len(transitions) < 15:
        raise RuntimeError("too few finite-error equality cases")
    for index, record in enumerate(transitions):
        rank = int(record["rank"])
        order = int(record["order"])
        epsilon_pass = float(record["epsilon_pass"])
        epsilon_stop = float(record["epsilon_stop"])
        angle = 0.5 * math.pi - math.asin(epsilon_pass) - math.asin(epsilon_stop)
        rho = math.sin(math.asin(epsilon_pass) + math.asin(epsilon_stop))
        if angle <= 0.0:
            raise RuntimeError(f"inadmissible error pair in transition {index}")
        expected_width = 2.0 * angle / order
        if not close(float(record["overlap_bound_rho"]), rho):
            raise RuntimeError(f"rho mismatch in transition {index}")
        if not close(float(record["observed_endpoint_overlap"]), rho):
            raise RuntimeError(f"endpoint overlap mismatch in transition {index}")
        if not close(float(record["principal_angle"]), angle):
            raise RuntimeError(f"principal angle mismatch in transition {index}")
        if not close(float(record["frequency_arc_width"]), expected_width):
            raise RuntimeError(f"arc width mismatch in transition {index}")
        trace_bound = 2.0 * rank * angle
        peak_bound = 2.0 * angle
        if not close(float(record["integrated_trace_delay"]), trace_bound):
            raise RuntimeError(f"trace equality failed in transition {index}")
        if not close(float(record["integrated_peak_delay"]), peak_bound):
            raise RuntimeError(f"peak equality failed in transition {index}")

    cyclic = payload["exact_cyclic_degree_saturation"]
    if len(cyclic) < 4:
        raise RuntimeError("too few cyclic saturation cases")
    for index, record in enumerate(cyclic):
        rank = int(record["rank"])
        order = int(record["order"])
        expected_degree = rank * order
        expected_action = 2.0 * math.pi * expected_degree
        if int(record["alternating_transitions"]) != 2 * order:
            raise RuntimeError(f"transition count failed in cyclic case {index}")
        if int(record["mcmillan_degree"]) != expected_degree:
            raise RuntimeError(f"degree count failed in cyclic case {index}")
        if not close(float(record["degree_lower_bound"]), expected_degree):
            raise RuntimeError(f"degree bound not saturated in cyclic case {index}")
        if not close(float(record["integrated_trace_delay"]), expected_action):
            raise RuntimeError(f"action identity failed in cyclic case {index}")
        if not close(float(record["routing_action_lower_bound"]), expected_action):
            raise RuntimeError(f"routing action not saturated in cyclic case {index}")

    block = payload["positive_block_stress_test"]
    if int(block["trials"]) < 2000:
        raise RuntimeError("positive-block stress test is too small")
    if not block["all_trace_block_bounds_pass"] or not block["all_peak_block_bounds_pass"]:
        raise RuntimeError("positive-block stress test reports a violation")
    if float(block["largest_two_nuclear_offdiag_over_trace"]) > 1.0 + 2.0e-12:
        raise RuntimeError("trace block inequality violated")
    if float(block["largest_two_operator_offdiag_over_peak"]) > 1.0 + 2.0e-12:
        raise RuntimeError("peak block inequality violated")

    paths = payload["random_potapov_path_stress_test"]
    if int(paths["trials"]) < 500 or int(paths["maximum_factor_count"]) < 12:
        raise RuntimeError("random path stress test is too small")
    if paths["subspace_ranks_tested"] != [1, 2, 3]:
        raise RuntimeError("random path ranks are incomplete")
    if not paths["all_random_paths_respect_speed_limit"]:
        raise RuntimeError("random Potapov path reports a speed-limit violation")
    ratio = float(paths["largest_distance_over_half_trace_action"])
    if ratio > 1.0 + 2.0e-9:
        raise RuntimeError("random Potapov path violates the speed limit")
    if ratio < 0.97:
        raise RuntimeError("stress test did not find a near-saturating path")

    figure = ROOT / "paper_routing_speed_limit" / "figures" / "routing_law.pdf"
    if not figure.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("routing figure is not a PDF")
    print("routing speed-limit certificate: PASS")


if __name__ == "__main__":
    main()
