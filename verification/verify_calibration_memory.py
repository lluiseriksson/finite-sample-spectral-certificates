#!/usr/bin/env python3
"""Fail-closed verification of the calibration-memory certificate."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results" / "calibration_memory" / "certificate.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    payload = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    require(payload["schema_version"] == 1, "unexpected schema")
    require(len(payload["records"]) >= 6, "scaling ledger is incomplete")
    for record in payload["records"]:
        order = record["order"]
        require(record["calibration_multiplicity"] == 3 * order, f"wrong calibration count at S={order}")
        require(record["full_network_mcmillan_degree"] == 3 * order + 6, f"wrong degree at S={order}")
        require(record["degree_overhead"] == 6, f"construction is not near-optimal at S={order}")
        require(record["loss_minor_calibration_zeros"] == 3 * order, f"loss zero count failed at S={order}")
        require(record["maximum_boundary_root_error"] < 1.0e-14, f"root audit failed at S={order}")
        require(record["maximum_loss_minor_zero_residual"] < 2.0e-10, f"loss zero residual failed at S={order}")
        require(record["strict_stop_loss_determinant"] > 0.5, f"loss minor may be trivial at S={order}")
        require(abs(record["wigner_smith_integral_over_2pi"] - (3 * order + 6)) < 2.0e-7, f"delay sum rule failed at S={order}")
        require(record["topological_lower_bound_satisfied"], f"topological bound failed at S={order}")
    rouche = payload["rouche_passive_perturbation"]
    require(rouche["maximum_rouche_ratio"] < 1.0, "Rouche inequality is not strict")
    require(rouche["all_contours_retain_one_zero"], "perturbed root left a certified contour")
    require(all(count == 1 for count in rouche["root_counts_per_contour"]), "wrong contour root multiplicity")
    random_minors = payload["random_potapov_minor_stress_test"]
    require(random_minors["trials"] >= 100, "random minor stress test is too small")
    require(random_minors["minor_sizes_tested"] == [1, 2, 3, 4, 5, 6], "not all compound ranks were tested")
    require(random_minors["all_minors_fit_the_factor_count_budget"], "a random minor exceeded its Potapov degree budget")
    require(random_minors["maximum_relative_interpolation_residual"] < 1.0e-8, "random minor interpolation residual is too large")
    print("calibration-memory certificate: PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"calibration-memory certificate: FAIL: {exc}", file=sys.stderr)
        raise
