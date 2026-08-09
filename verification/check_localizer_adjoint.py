"""Independent randomized check of the scaled block-localizer adjoint identity."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def localizer(q: np.ndarray, degree: int, channels: int, theta: float, scales: np.ndarray) -> np.ndarray:
    d = 2 * degree + 2
    k = np.diag(scales) @ q @ np.diag(scales)
    blocks = [k[t * channels : (t + 1) * channels, t * channels : (t + 1) * channels] for t in range(d)]
    return np.block(
        [
            [theta * blocks[i + j] - blocks[i + j + 1] for j in range(degree + 1)]
            for i in range(degree + 1)
        ]
    )


def adjoint(y: np.ndarray, degree: int, channels: int, theta: float, scales: np.ndarray) -> np.ndarray:
    d = 2 * degree + 2
    answer = np.zeros((d * channels, d * channels))
    for k in range(d):
        block = np.zeros((channels, channels))
        for i in range(degree + 1):
            for j in range(degree + 1):
                coefficient = theta * (i + j == k) - (i + j + 1 == k)
                block += coefficient * y[
                    i * channels : (i + 1) * channels,
                    j * channels : (j + 1) * channels,
                ]
        indices = slice(k * channels, (k + 1) * channels)
        diagonal = np.diag(scales[indices])
        answer[indices, indices] = diagonal @ block @ diagonal
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260816)
    parser.add_argument("--output", type=Path, default=Path("results/verification/localizer_adjoint.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    errors = []
    relative_errors = []
    cases = []
    for trial in range(args.trials):
        degree = int(rng.integers(1, 5))
        channels = int(rng.integers(1, 4))
        theta = float(rng.uniform(0.1, 0.95))
        size = (2 * degree + 2) * channels
        localizer_size = (degree + 1) * channels
        q = rng.normal(size=(size, size))
        q = (q + q.T) / 2
        y = rng.normal(size=(localizer_size, localizer_size))
        y = (y + y.T) / 2
        scales = np.exp(rng.uniform(-3.0, 3.0, size=size))
        lhs = float(np.trace(y @ localizer(q, degree, channels, theta, scales)))
        rhs = float(np.trace(adjoint(y, degree, channels, theta, scales) @ q))
        error = abs(lhs - rhs)
        relative = error / max(1.0, abs(lhs), abs(rhs))
        errors.append(error)
        relative_errors.append(relative)
        if trial < 5:
            cases.append({"degree": degree, "channels": channels, "theta": theta, "lhs": lhs, "rhs": rhs})
    payload = {
        "schema_version": 1,
        "setup": vars(args) | {"output": str(args.output)},
        "maximum_absolute_error": max(errors),
        "maximum_relative_error": max(relative_errors),
        "first_five": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["maximum_relative_error"] > 2e-12:
        raise SystemExit("adjoint identity failed")


if __name__ == "__main__":
    main()
