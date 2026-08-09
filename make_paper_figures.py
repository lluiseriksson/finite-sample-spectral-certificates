"""Freeze production summaries and generate publication figures."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter


ROOT = Path(__file__).resolve().parent
PRODUCTION = ROOT / "results" / "production"
FIGURES = ROOT / "paper" / "figures"


def load(name: str) -> dict:
    return json.loads((PRODUCTION / name).read_text(encoding="utf-8"))


def digest(name: str) -> str:
    return hashlib.sha256((PRODUCTION / name).read_bytes()).hexdigest()


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - radius, center + radius


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    main_campaign = load("wishart_campaign_colab_repaired.json")
    ablation = load("degree_ablation_colab.json")
    misspecification = load("primitive_misspecification_colab.json")

    summary: dict[str, object] = {
        "schema_version": 1,
        "hashes": {
            name: digest(name)
            for name in (
                "wishart_campaign_colab_repaired.json",
                "degree_ablation_colab.json",
                "primitive_misspecification_colab.json",
            )
        },
    }
    main_cells = []
    null_successes = alternative_35_successes = 0
    null_total = alternative_35_total = 0
    maximum_residual_bound = 0.0
    minimum_psd_eigenvalue = math.inf
    for campaign in main_campaign["campaigns"]:
        length = campaign["setup"]["length"]
        for row in campaign["rows"]:
            successes = round(row["detection_rate"] * campaign["setup"]["ensembles"])
            total = campaign["setup"]["ensembles"]
            if row["gap_factor"] == 1.0:
                null_successes += successes
                null_total += total
            if row["gap_factor"] == 1.35:
                alternative_35_successes += successes
                alternative_35_total += total
            maximum_residual_bound = max(maximum_residual_bound, row["maximum_safe_stationarity_bound"])
            minimum_psd_eigenvalue = min(minimum_psd_eigenvalue, row["minimum_safe_psd_eigenvalue"])
            if row["gap_factor"] == 1.15:
                low, high = wilson(successes, total)
                main_cells.append(
                    {
                        "length": length,
                        "sample_count": row["sample_count"],
                        "successes": successes,
                        "total": total,
                        "rate": row["detection_rate"],
                        "wilson_95": [low, high],
                    }
                )
    summary["main_campaign"] = {
        "null_rejections": [null_successes, null_total],
        "inflated_35_rejections": [alternative_35_successes, alternative_35_total],
        "maximum_safe_stationarity_bound": maximum_residual_bound,
        "minimum_safe_psd_eigenvalue": minimum_psd_eigenvalue,
        "inflated_15_cells": main_cells,
        "finite_volume_gaps": [
            {"length": c["setup"]["length"], "gap": c["model"]["true_gap"]}
            for c in main_campaign["campaigns"]
        ],
    }

    ablation_cells = []
    ablation_null_successes = ablation_null_total = 0
    for campaign in ablation["campaigns"]:
        length, degree = campaign["setup"]["length"], campaign["setup"]["degree"]
        total = campaign["setup"]["ensembles"]
        for row in campaign["rows"]:
            successes = round(row["detection_rate"] * total)
            if row["gap_factor"] == 1.0:
                ablation_null_successes += successes
                ablation_null_total += total
            else:
                low, high = wilson(successes, total)
                ablation_cells.append(
                    {
                        "length": length,
                        "degree": degree,
                        "sample_count": row["sample_count"],
                        "successes": successes,
                        "total": total,
                        "rate": row["detection_rate"],
                        "wilson_95": [low, high],
                    }
                )
    summary["degree_ablation"] = {
        "null_rejections": [ablation_null_successes, ablation_null_total],
        "inflated_15_cells": ablation_cells,
    }
    summary["misspecification"] = misspecification["rows"]
    (PRODUCTION / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    colors = plt.get_cmap("viridis")(np.linspace(0.08, 0.9, 6))
    fig = plt.figure(figsize=(7.1, 6.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[1.18, 1.0])
    axis_main = fig.add_subplot(grid[0, :])
    for color, length in zip(colors, sorted({cell["length"] for cell in main_cells})):
        cells = sorted((cell for cell in main_cells if cell["length"] == length), key=lambda x: x["sample_count"])
        x = np.array([cell["sample_count"] for cell in cells])
        y = np.array([cell["rate"] for cell in cells])
        lower = np.maximum(0.0, y - np.array([cell["wilson_95"][0] for cell in cells]))
        upper = np.maximum(0.0, np.array([cell["wilson_95"][1] for cell in cells]) - y)
        axis_main.errorbar(x, y, yerr=[lower, upper], marker="o", linewidth=1.45, capsize=2.5, color=color, label=f"L={length}")
    axis_main.set_xscale("log")
    axis_main.xaxis.set_minor_formatter(NullFormatter())
    axis_main.set_xticks([500, 1000, 2000, 5000], labels=["500", "1000", "2000", "5000"])
    axis_main.set_ylim(-0.04, 1.06)
    axis_main.set_ylabel("rejection probability")
    axis_main.set_xlabel("Gaussian sketches n")
    axis_main.set_title("(a) 15% inflated visible-gap claim, moment degree N=2", loc="left", fontsize=10)
    axis_main.grid(alpha=0.22)
    axis_main.legend(ncol=3, frameon=False, fontsize=8)

    degree_colors = ["#0072B2", "#D55E00", "#009E73"]
    for column, length in enumerate([8, 12, 16]):
        axis = fig.add_subplot(grid[1, column])
        for color, degree in zip(degree_colors, [1, 2, 3]):
            cells = sorted(
                (cell for cell in ablation_cells if cell["length"] == length and cell["degree"] == degree),
                key=lambda x: x["sample_count"],
            )
            x = [cell["sample_count"] for cell in cells]
            y = [cell["rate"] for cell in cells]
            axis.plot(x, y, marker="o", linewidth=1.4, color=color, label=f"N={degree}")
        axis.set_xscale("log")
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.set_xticks([500, 1000, 2000, 5000], labels=["0.5k", "1k", "2k", "5k"])
        axis.set_ylim(-0.04, 1.06)
        axis.set_title(f"({chr(98 + column)}) L={length}", loc="left", fontsize=9)
        axis.set_xlabel("n")
        axis.grid(alpha=0.22)
        if column == 0:
            axis.set_ylabel("rejection probability")
        else:
            axis.set_yticklabels([])
        if column == 2:
            axis.legend(frameon=False, fontsize=8, loc="lower right")
    for extension in ("pdf", "png"):
        fig.savefig(FIGURES / f"power_and_degree.{extension}", dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6.7, 3.25), constrained_layout=True)
    laws = ["gaussian", "student", "rademacher"]
    labels = ["Gaussian\n(in model)", "Student-$t_5$\n(out of model)", "Rademacher\n(out of model)"]
    positions = np.arange(len(laws))
    width = 0.34
    for offset, count, color in [(-width / 2, 500, "#56B4E9"), (width / 2, 1000, "#E69F00")]:
        rates = []
        errors_low, errors_high = [], []
        for law in laws:
            rows = [row for row in misspecification["rows"] if row["law"] == law and row["sample_count"] == count]
            successes = round(sum(row["coverage_rate"] for row in rows) * 200)
            total = len(rows) * 200
            rate = successes / total
            low, high = wilson(successes, total)
            rates.append(rate)
            errors_low.append(max(0.0, rate - low))
            errors_high.append(max(0.0, high - rate))
        axis.bar(positions + offset, rates, width, color=color, label=f"n={count}", yerr=[errors_low, errors_high], capsize=3)
    axis.axhline(0.95, color="black", linewidth=1.0, linestyle="--", label="nominal 95%")
    axis.set_xticks(positions, labels)
    axis.set_ylim(0.64, 1.025)
    axis.set_ylabel("simultaneous-band coverage")
    axis.set_title("Primitive-law stress test (400 ensembles per bar)", loc="left", fontsize=10)
    axis.legend(frameon=False, ncol=3, fontsize=8, loc="lower left")
    axis.grid(axis="y", alpha=0.22)
    for extension in ("pdf", "png"):
        fig.savefig(FIGURES / f"misspecification_coverage.{extension}", dpi=240, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary["main_campaign"], indent=2))
    print(f"wrote {PRODUCTION / 'summary.json'} and figures in {FIGURES}")


if __name__ == "__main__":
    main()
