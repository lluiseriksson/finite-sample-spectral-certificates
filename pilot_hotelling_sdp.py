"""Unknown-covariance pilot using the exact one-sample Hotelling ellipsoid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import f, t

from pilot_covariance_sdp import (
    covariance_model,
    localizer_matrix,
    margins,
    moments,
)


def campaign(
    rng: np.random.Generator,
    true_moments: np.ndarray,
    true_covariance: np.ndarray,
    theta: float,
    target: float,
    degree: int,
    coverage: float,
    sample_count: int,
    trials: int,
) -> dict[str, object]:
    dimension = true_covariance.shape[0]
    alpha = 1.0 - coverage
    hotelling_q_sq = (
        dimension
        * (sample_count - 1)
        / (sample_count - dimension)
        * f.ppf(coverage, dimension, sample_count - dimension)
    )
    hotelling_q = float(np.sqrt(hotelling_q_sq))
    bonferroni_t = float(t.ppf(1.0 - alpha / (2.0 * dimension), sample_count - 1))
    detected = {"fixed_hotelling": 0, "fixed_bonferroni": 0, "optimized_hotelling": 0}
    ellipsoid_misses = 0
    optimized_inside = 0
    examples = []

    for trial in range(trials):
        noise = rng.multivariate_normal(np.zeros(dimension), true_covariance, size=sample_count)
        sample_mean = noise.mean(axis=0)
        sample_covariance = np.cov(noise, rowvar=False, ddof=1)
        observed = true_moments.copy()
        observed[1:] += sample_mean
        covariance_of_mean_metric = sample_covariance / sample_count
        result = margins(
            observed,
            covariance_of_mean_metric,
            theta,
            target,
            degree,
            coverage,
            ellipsoid_q=hotelling_q,
            box_multiplier=bonferroni_t,
        )
        detected["fixed_hotelling"] += result["fixed_chebyshev_ellipsoid_upper"] < 0
        detected["fixed_bonferroni"] += result["fixed_chebyshev_box_upper"] < 0
        detected["optimized_hotelling"] += result["optimized_ellipsoid_upper"] < 0

        hotelling_t_sq = float(
            sample_count * sample_mean @ np.linalg.solve(sample_covariance, sample_mean)
        )
        inside = hotelling_t_sq <= hotelling_q_sq
        ellipsoid_misses += not inside
        optimized_inside += inside and result["optimized_ellipsoid_upper"] < -1e-7
        if trial < 3:
            x_value = np.asarray(result["optimized_lift"])
            result["true_localizer_at_optimized_lift"] = float(
                np.trace(localizer_matrix(true_moments, theta, degree) @ x_value)
            )
            result["hotelling_t_sq"] = hotelling_t_sq
            result["hotelling_q_sq"] = hotelling_q_sq
            result["inside_hotelling_ellipsoid"] = inside
            examples.append(result)

    return {
        "sample_count": sample_count,
        "trials": trials,
        "hotelling_q_sq": hotelling_q_sq,
        "detections": detected,
        "rates": {key: value / trials for key, value in detected.items()},
        "ellipsoid_miss_rate": ellipsoid_misses / trials,
        "optimized_detections_inside_ellipsoid": optimized_inside,
        "first_three": examples,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=5)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--trials", type=int, default=500)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--sigma", type=float, default=5e-6)
    parser.add_argument("--correlation", type=float, default=0.92)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("pilot_hotelling_sdp.json"))
    args = parser.parse_args()
    if args.samples <= 2 * args.degree + 1:
        raise SystemExit("Hotelling requires samples > moment dimension")

    theta, delta, gamma = 0.8, 0.1, 0.01
    low_atoms = np.array([0.15, 0.45, 0.78])
    low_weights = np.array([0.2, 0.3, 0.5])
    null_moments = moments(np.array([theta]), np.array([1.0]), 2 * args.degree + 2)
    alternative_moments = moments(
        np.append(low_atoms, theta + delta),
        np.append((1.0 - gamma) * low_weights, gamma),
        2 * args.degree + 2,
    )
    covariance = covariance_model(2 * args.degree + 1, args.sigma, args.correlation)
    rng = np.random.default_rng(args.seed)
    payload = {
        "setup": vars(args) | {"output": str(args.output), "theta": theta, "delta": delta, "gamma": gamma},
        "null_boundary": campaign(
            rng,
            null_moments,
            covariance,
            theta,
            theta + delta,
            args.degree,
            args.coverage,
            args.samples,
            args.trials,
        ),
        "visible_alternative": campaign(
            rng,
            alternative_moments,
            covariance,
            theta,
            theta + delta,
            args.degree,
            args.coverage,
            args.samples,
            args.trials,
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(payload["null_boundary"]["rates"], indent=2))
    print(json.dumps(payload["visible_alternative"]["rates"], indent=2))


if __name__ == "__main__":
    main()
