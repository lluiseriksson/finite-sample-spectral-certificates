#!/usr/bin/env python3
"""Independent structural verifier for the detector-rank hierarchy fixtures."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "results/unified_routing_table_memory/detector_rank_certificate.json"


def main() -> None:
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    if certificate["schema"] != "detector-rank-hierarchy-certificate-v1":
        raise RuntimeError("detector-rank schema mismatch")
    if certificate["theorem"] != "d_min >= max_D sum_a n_a * (k - rank(D|Y_a))":
        raise RuntimeError("detector-rank theorem mismatch")
    body = dict(certificate)
    digest = body.pop("content_digest")
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    if hashlib.sha256(canonical).hexdigest() != digest:
        raise RuntimeError("detector-rank digest mismatch")

    scalar_expected = {
        "all_coincident": (1, 0),
        "three_distinct_coplanar": (2, 1),
        "two_distinct_one_repeated": (2, 2),
        "three_independent": (3, 2),
        "rarest_target_word_AAAB": (2, 3),
        "four_generic_lines_in_a_plane": (2, 1),
    }
    scalar = certificate["scalar_hyperplane_cases"]
    if set(scalar) != set(scalar_expected):
        raise RuntimeError("scalar fixture set mismatch")
    for name, (span_rank, bound) in scalar_expected.items():
        if scalar[name]["span_rank"] != span_rank:
            raise RuntimeError(f"span rank mismatch for {name}")
        if scalar[name]["separator_bound"] != bound:
            raise RuntimeError(f"separator bound mismatch for {name}")

    block = certificate["rank_two_detector_witnesses"]
    if block["direct_sum_AAAB"] != {
        "detector_ranks": {"A": 0, "B": 2},
        "weighted_rank_deficiency": 6,
    }:
        raise RuntimeError("direct-sum block witness mismatch")
    if block["three_planes_with_common_line"] != {
        "detector_ranks": {"A": 1, "B": 2, "C": 1},
        "weighted_rank_deficiency": 4,
    }:
        raise RuntimeError("intersecting block witness mismatch")
    print("detector-rank hierarchy verifier: PASS")


if __name__ == "__main__":
    main()
