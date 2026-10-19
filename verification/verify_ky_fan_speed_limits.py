#!/usr/bin/env python3
"""Independent checks for the frozen Ky Fan speed-limit certificate."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "ky_fan_speed_limits" / "certificate.json"
TOL = 5.0e-9


def close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=TOL, abs_tol=TOL)


def direct_pointwise_probe() -> None:
    """Use a separate block construction to probe every partial sum."""
    rng = np.random.default_rng(781229)
    for signed in (False, True):
        for _ in range(512):
            size = int(rng.integers(3, 15))
            rank = int(rng.integers(1, size))
            count = min(rank, size - rank)
            raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
            if signed:
                q = raw + raw.conj().T
            else:
                q = raw @ raw.conj().T
            singular = np.linalg.svd(q[rank:, :rank], compute_uv=False)[:count]
            eigenvalues = np.linalg.eigvalsh(q)
            if signed:
                budget = eigenvalues[::-1][:count] - eigenvalues[:count]
            else:
                budget = eigenvalues[::-1][:count]
            if np.any(2.0 * np.cumsum(singular) > np.cumsum(budget) + TOL):
                raise RuntimeError("independent pointwise Ky Fan probe failed")


def main() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if payload.get("schema") != "ky-fan-speed-limits/v1":
        raise RuntimeError("unexpected certificate schema")
    if not payload.get("all_checks_pass"):
        raise RuntimeError("frozen certificate reports a failure")

    sharp = payload["sharp_equalities"]
    if len(sharp["cases"]) < 4 or not sharp["all_pass"]:
        raise RuntimeError("sharpness suite is incomplete")
    for case in sharp["cases"]:
        angles = np.asarray(case["angles"], dtype=float)
        observed = np.asarray(case["observed_angles"], dtype=float)
        required = 2.0 * np.cumsum(angles)
        delays = np.asarray(case["top_delay_partial_sums"], dtype=float)
        if not np.allclose(angles, observed, rtol=TOL, atol=TOL):
            raise RuntimeError("canonical-angle equality family failed")
        if not np.allclose(required, delays, rtol=TOL, atol=TOL):
            raise RuntimeError("proper-delay equality family failed")
        for name, gauge in case["symmetric_gauge_costs"].items():
            if not close(float(gauge["geometric_optimum"]), float(gauge["delay_action"])):
                raise RuntimeError(f"symmetric-gauge variational equality failed: {name}")

    separation = payload["scalar_endpoint_separation"]
    if not separation["same_scalar_endpoints"] or not separation["strict_intermediate_separation"]:
        raise RuntimeError("scalar endpoints failed to hide a strict intermediate-prefix gap")
    first = np.asarray(separation["first_angles"], dtype=float)
    second = np.asarray(separation["second_angles"], dtype=float)
    if not close(float(first[0]), float(second[0])) or not close(float(np.sum(first)), float(np.sum(second))):
        raise RuntimeError("separation example does not share max and total angle")
    if not close(float(separation["r2_action_gap"]), 2.0 * float(first[:2].sum() - second[:2].sum())):
        raise RuntimeError("intermediate-prefix separation arithmetic failed")

    for key, minimum in (
        ("signed_pointwise", 4000),
        ("positive_pointwise", 4000),
        ("signed_paths", 700),
        ("positive_paths", 700),
        ("tomography", 2000),
        ("slack_decomposition", 4000),
    ):
        block = payload[key]
        count = block.get("matrices", block.get("piecewise_constant_paths", block.get("trials", 0)))
        if int(count) < minimum or not block.get("all_pass"):
            raise RuntimeError(f"insufficient or failed suite: {key}")

    for key in ("signed_pointwise", "positive_pointwise", "signed_paths", "positive_paths"):
        block = payload[key]
        if float(block["worst_lhs_over_rhs"]) > 1.0 + TOL:
            raise RuntimeError(f"hierarchy violation in {key}")
        if float(block["smallest_absolute_slack"]) < -TOL:
            raise RuntimeError(f"negative slack in {key}")
    if float(payload["tomography"]["worst_certified_minus_actual"]) > TOL:
        raise RuntimeError("measurement certificate overstates the true guarantee")
    if not close(float(payload["slack_decomposition"]["worst_identity_error"]), 0.0):
        raise RuntimeError("slack identity failed")

    direct_pointwise_probe()
    figure = ROOT / "paper_ky_fan_speed_limits" / "figures" / "ky_fan_hierarchy.pdf"
    if not figure.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("hierarchy figure is not a PDF")
    print("Ky Fan speed-limit certificate: PASS")


if __name__ == "__main__":
    main()
