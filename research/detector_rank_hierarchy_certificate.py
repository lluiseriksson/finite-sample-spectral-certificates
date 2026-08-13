#!/usr/bin/env python3
"""Exact fixtures for the detector-rank hierarchy.

The calculations use rational Gaussian elimination only.  They audit the
finite-dimensional incidence layer of the theorem; the analytic zero budget
remains in the paper.
"""

from __future__ import annotations

from fractions import Fraction
import argparse
import hashlib
import itertools
import json
from pathlib import Path


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def rank(columns: tuple[Vector, ...]) -> int:
    if not columns:
        return 0
    rows = [[Fraction(column[row]) for column in columns] for row in range(len(columns[0]))]
    pivot_row = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(pivot_row, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [entry / scale for entry in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [left - factor * right for left, right in zip(rows[row], rows[pivot_row], strict=True)]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def in_span(vector: Vector, columns: tuple[Vector, ...]) -> bool:
    return rank(columns + (vector,)) == rank(columns)


def scalar_separator(vectors: dict[str, Vector], occupancies: dict[str, int]) -> dict:
    labels = tuple(sorted(vectors))
    ambient_rank = rank(tuple(vectors[label] for label in labels))
    best_weight = -1
    best_closures: set[tuple[str, ...]] = set()
    for size in range(len(labels) + 1):
        for subset in itertools.combinations(labels, size):
            span = tuple(vectors[label] for label in subset)
            if rank(span) >= ambient_rank:
                continue
            closure = tuple(label for label in labels if in_span(vectors[label], span))
            weight = sum(occupancies[label] for label in closure)
            if weight > best_weight:
                best_weight = weight
                best_closures = {closure}
            elif weight == best_weight:
                best_closures.add(closure)
    return {
        "span_rank": ambient_rank,
        "separator_bound": best_weight,
        "maximizing_closures": [list(closure) for closure in sorted(best_closures)],
    }


def matmul(left: Matrix, right_columns: tuple[Vector, ...]) -> tuple[Vector, ...]:
    return tuple(
        tuple(sum(left_row[index] * column[index] for index in range(len(column))) for left_row in left)
        for column in right_columns
    )


def detector_weight(
    frames: dict[str, tuple[Vector, ...]],
    occupancies: dict[str, int],
    detector: Matrix,
) -> dict:
    target_rank = len(detector)
    ranks = {label: rank(matmul(detector, frame)) for label, frame in frames.items()}
    if max(ranks.values()) != target_rank:
        raise ValueError("detector is not full rank on any target")
    return {
        "detector_ranks": ranks,
        "weighted_rank_deficiency": sum(
            occupancies[label] * (target_rank - ranks[label]) for label in sorted(frames)
        ),
    }


def make_certificate() -> dict:
    scalar_cases = {
        "all_coincident": (
            {"A": (1, 0, 0), "B": (1, 0, 0), "C": (1, 0, 0)},
            {"A": 1, "B": 1, "C": 1},
        ),
        "three_distinct_coplanar": (
            {"A": (1, 0, 0), "B": (0, 1, 0), "C": (1, 1, 0)},
            {"A": 1, "B": 1, "C": 1},
        ),
        "two_distinct_one_repeated": (
            {"A": (1, 0, 0), "B": (1, 0, 0), "C": (0, 1, 0)},
            {"A": 1, "B": 1, "C": 1},
        ),
        "three_independent": (
            {"A": (1, 0, 0), "B": (0, 1, 0), "C": (0, 0, 1)},
            {"A": 1, "B": 1, "C": 1},
        ),
        "rarest_target_word_AAAB": (
            {"A": (1, 0), "B": (0, 1)},
            {"A": 3, "B": 1},
        ),
        "four_generic_lines_in_a_plane": (
            {"A": (1, 0), "B": (0, 1), "C": (1, 1), "D": (1, 2)},
            {"A": 1, "B": 1, "C": 1, "D": 1},
        ),
    }
    scalar = {
        name: {"occupancies": occupancies, **scalar_separator(vectors, occupancies)}
        for name, (vectors, occupancies) in scalar_cases.items()
    }

    block_frames = {
        "A": ((1, 0, 0, 0), (0, 1, 0, 0)),
        "B": ((0, 0, 1, 0), (0, 0, 0, 1)),
    }
    block_direct_sum = detector_weight(
        block_frames,
        {"A": 3, "B": 1},
        ((0, 0, 1, 0), (0, 0, 0, 1)),
    )

    intersecting_frames = {
        "A": ((1, 0, 0, 0), (0, 1, 0, 0)),
        "B": ((1, 0, 0, 0), (0, 0, 1, 0)),
        "C": ((1, 0, 0, 0), (0, 0, 0, 1)),
    }
    block_intersection = detector_weight(
        intersecting_frames,
        {"A": 3, "B": 1, "C": 1},
        ((1, 0, 0, 0), (0, 0, 1, 0)),
    )

    payload = {
        "schema": "detector-rank-hierarchy-certificate-v1",
        "theorem": "d_min >= max_D sum_a n_a * (k - rank(D|Y_a))",
        "arithmetic": "exact rational Gaussian elimination",
        "scalar_hyperplane_cases": scalar,
        "rank_two_detector_witnesses": {
            "direct_sum_AAAB": block_direct_sum,
            "three_planes_with_common_line": block_intersection,
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["content_digest"] = hashlib.sha256(canonical).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    certificate = make_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "scalar_bounds": {
            name: case["separator_bound"]
            for name, case in certificate["scalar_hyperplane_cases"].items()
        },
        "block_bounds": {
            name: case["weighted_rank_deficiency"]
            for name, case in certificate["rank_two_detector_witnesses"].items()
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()
