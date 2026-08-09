"""ANNNI multichannel correlators with non-Gaussian Gaussian-sketch samples.

For vectors u_{k,a}=T^(k/2) psi_a and one standard Gaussian g, define
z_{k,a}=<g,u_{k,a}>.  Then z_k z_k^T is an unbiased rank-one estimator of
B_k.  Reusing g across k produces the full physical time/channel covariance;
the outer products are Wishart-like rather than multivariate normal.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.stats import f

from pilot_local_visibility import hamiltonian, parity


def probes_from_ground(ground: np.ndarray, length: int) -> np.ndarray:
    states = np.arange(1 << length)
    center = length // 2
    x_ground = ground[states ^ (1 << center)]
    x_ground -= np.vdot(ground, x_ground).real * ground
    z = 1.0 - 2.0 * ((states >> center) & 1)
    z_ground = z * ground
    z_ground -= np.vdot(ground, z_ground).real * ground
    return np.column_stack([x_ground, z_ground])


def exact_half_moments(
    length: int, degree: int, tau: float, j1: float, j2: float, hx: float
) -> tuple[np.ndarray, dict[str, object]]:
    h = hamiltonian(length, j1, j2, hx)
    energies, vectors = eigsh(h, k=3, which="SA", tol=1e-11, maxiter=300_000)
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    ground = vectors[:, 0]
    probes = probes_from_ground(ground, length)
    max_moment = 2 * degree + 1
    generator = -(h - energies[0] * sparse.eye(h.shape[0], format="csr"))
    evolved = expm_multiply(
        generator,
        probes,
        start=0.0,
        stop=tau * max_moment,
        num=2 * max_moment + 1,
        endpoint=True,
        traceA=energies[0] * h.shape[0],
    )
    half = np.array([probes.T @ item for item in evolved])
    half = (half + half.transpose(0, 2, 1)) / 2
    metadata = {
        "length": length,
        "dimension": 1 << length,
        "parameters": {"j1": j1, "j2": j2, "hx": hx, "tau": tau},
        "energies": energies.tolist(),
        "true_gap": float(energies[1] - energies[0]),
        "parities": [parity(vectors[:, i]) for i in range(3)],
    }
    return half, metadata


def joint_sketch_covariance(half: np.ndarray, max_moment: int) -> np.ndarray:
    channels = half.shape[1]
    covariance = np.empty((max_moment + 1, channels, max_moment + 1, channels))
    for i in range(max_moment + 1):
        for j in range(max_moment + 1):
            covariance[i, :, j, :] = half[i + j]
    covariance = covariance.reshape((max_moment + 1) * channels, -1)
    covariance = (covariance + covariance.T) / 2
    values, vectors = np.linalg.eigh(covariance)
    return (vectors * np.maximum(values, 0.0)) @ vectors.T


def parameter_index(max_moment: int, channels: int) -> list[tuple[int, int, int]]:
    return [
        (k, a, b)
        for k in range(max_moment + 1)
        for a in range(channels)
        for b in range(a, channels)
    ]


def flatten_moments(matrices: np.ndarray) -> np.ndarray:
    return np.array(
        [matrices[k, a, b] for k, a, b in parameter_index(len(matrices) - 1, matrices.shape[1])]
    )


def sketch_feature_covariance(joint: np.ndarray, max_moment: int, channels: int) -> np.ndarray:
    """Exact covariance of flattened products z_{k,a} z_{k,b} by Isserlis."""
    params = parameter_index(max_moment, channels)
    answer = np.empty((len(params), len(params)))
    for i, (k, a, b) in enumerate(params):
        ka, kb = k * channels + a, k * channels + b
        for j, (ell, c, d) in enumerate(params):
            lc, ld = ell * channels + c, ell * channels + d
            answer[i, j] = joint[ka, lc] * joint[kb, ld] + joint[ka, ld] * joint[kb, lc]
    return (answer + answer.T) / 2


def block_design(theta: float, degree: int, channels: int) -> np.ndarray:
    maximum = 2 * degree + 1
    params = parameter_index(maximum, channels)
    size = (degree + 1) * channels
    design = np.zeros((len(params), size, size))
    for parameter, (k, a, b) in enumerate(params):
        for i in range(degree + 1):
            for j in range(degree + 1):
                coefficient = theta * (i + j == k) - (i + j + 1 == k)
                if coefficient == 0:
                    continue
                design[parameter, i * channels + a, j * channels + b] += coefficient
                if a != b:
                    design[parameter, i * channels + b, j * channels + a] += coefficient
    return design


def robust_tests(
    observed: np.ndarray,
    covariance_mean: np.ndarray,
    design: np.ndarray,
    q: float,
) -> dict[str, object]:
    matrix = np.einsum("k,kij->ij", observed, design)
    root = np.linalg.cholesky(covariance_mean + 1e-16 * np.eye(len(observed)))

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    direction = eigenvectors[:, 0]
    rank_one = np.outer(direction, direction)
    l_rank_one = np.einsum("kij,ij->k", design, rank_one)
    rank_one_upper = float(
        eigenvalues[0] + q * np.linalg.norm(root.T @ l_rank_one)
    )

    x = cp.Variable(matrix.shape, PSD=True)
    l_expression = cp.hstack([cp.trace(item @ x) for item in design])
    objective = cp.trace(matrix @ x) + q * cp.norm(root.T @ l_expression, 2)
    problem = cp.Problem(cp.Minimize(objective), [cp.trace(x) == 1.0])
    problem.solve(solver="CLARABEL", tol_gap_abs=1e-9, tol_feas=1e-9, max_iter=500)
    if problem.status not in {"optimal", "optimal_inaccurate"}:
        raise RuntimeError(problem.status)
    raw = (x.value + x.value.T) / 2
    values, vectors = np.linalg.eigh(raw)
    feasible = (vectors * np.maximum(values, 0.0)) @ vectors.T
    feasible /= np.trace(feasible)
    l_feasible = np.einsum("kij,ij->k", design, feasible)
    measured = float(np.trace(matrix @ feasible))
    margin = q * float(np.linalg.norm(root.T @ l_feasible))
    return {
        "observed_min_eigenvalue": float(eigenvalues[0]),
        "adaptive_rank_one_upper": rank_one_upper,
        "optimized_psd_upper": measured + margin,
        "optimized_psd_measured": measured,
        "optimized_psd_margin": margin,
        "optimized_psd_rank_1e-7": int(np.count_nonzero(np.linalg.eigvalsh(feasible) > 1e-7)),
    }


def one_ensemble(
    rng: np.random.Generator,
    joint_covariance: np.ndarray,
    feature_covariance: np.ndarray,
    exact_flat: np.ndarray,
    design: np.ndarray,
    sample_count: int,
    coverage: float,
) -> dict[str, object]:
    channels = 2
    maximum = joint_covariance.shape[0] // channels - 1
    samples_z = rng.multivariate_normal(np.zeros(joint_covariance.shape[0]), joint_covariance, size=sample_count)
    samples_z = samples_z.reshape(sample_count, maximum + 1, channels)
    samples = np.empty((sample_count, len(exact_flat)))
    for r in range(sample_count):
        outer = np.einsum("ka,kb->kab", samples_z[r], samples_z[r])
        samples[r] = flatten_moments(outer)
    mean = samples.mean(axis=0)
    covariance = np.cov(samples, rowvar=False, ddof=1)
    dimension = len(mean)
    q_sq = dimension * (sample_count - 1) / (sample_count - dimension) * f.ppf(
        coverage, dimension, sample_count - dimension
    )
    error = mean - exact_flat
    t_sq = float(sample_count * error @ np.linalg.solve(covariance, error))
    result = robust_tests(mean, covariance / sample_count, design, float(np.sqrt(q_sq)))
    result["hotelling_t_sq"] = t_sq
    result["hotelling_q_sq"] = float(q_sq)
    result["inside_hotelling_ellipsoid"] = t_sq <= q_sq
    # Exact distribution-free multivariate Chebyshev ellipsoid with known
    # covariance of the sketch features.  E[Mahalanobis^2]=dimension.
    distribution_free_q_sq = dimension / (1.0 - coverage)
    covariance_of_mean = feature_covariance / sample_count
    distribution_free_t_sq = float(error @ np.linalg.solve(covariance_of_mean, error))
    distribution_free = robust_tests(
        mean, covariance_of_mean, design, float(np.sqrt(distribution_free_q_sq))
    )
    result["distribution_free_rank_one_upper"] = distribution_free["adaptive_rank_one_upper"]
    result["distribution_free_psd_upper"] = distribution_free["optimized_psd_upper"]
    result["distribution_free_t_sq"] = distribution_free_t_sq
    result["distribution_free_q_sq"] = distribution_free_q_sq
    result["inside_distribution_free_ellipsoid"] = distribution_free_t_sq <= distribution_free_q_sq
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=12)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--samples", type=int, default=120)
    parser.add_argument("--ensembles", type=int, default=300)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--gap-factor", type=float, default=1.35)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("pilot_annni_block_hotelling.json"))
    args = parser.parse_args()
    tau, j1, j2, hx = 0.2, 1.0, 0.37, 2.2
    half, metadata = exact_half_moments(args.length, args.degree, tau, j1, j2, hx)
    maximum = 2 * args.degree + 1
    exact_matrices = half[::2][: maximum + 1]
    exact_flat = flatten_moments(exact_matrices)
    joint = joint_sketch_covariance(half, maximum)
    feature_covariance = sketch_feature_covariance(joint, maximum, channels=2)
    false_gap_factor = args.gap_factor
    theta = float(np.exp(-tau * false_gap_factor * metadata["true_gap"]))
    design = block_design(theta, args.degree, channels=2)
    exact_localizer = np.einsum("k,kij->ij", exact_flat, design)
    rng = np.random.default_rng(args.seed)
    rows = [
        one_ensemble(
            rng,
            joint,
            feature_covariance,
            exact_flat,
            design,
            args.samples,
            args.coverage,
        )
        for _ in range(args.ensembles)
    ]
    payload = {
        "setup": vars(args) | {"output": str(args.output), "theta": theta, "false_gap_factor": false_gap_factor},
        "model": metadata,
        "exact": {
            "localizer_min_eigenvalue": float(np.linalg.eigvalsh(exact_localizer)[0]),
            "moment_parameter_dimension": len(exact_flat),
            "joint_sketch_covariance_rank": int(np.linalg.matrix_rank(joint, tol=1e-10)),
            "feature_covariance_rank": int(np.linalg.matrix_rank(feature_covariance, tol=1e-10)),
        },
        "summary": {
            "hotelling_coverage_rate": float(np.mean([row["inside_hotelling_ellipsoid"] for row in rows])),
            "rank_one_detection_rate": float(np.mean([row["adaptive_rank_one_upper"] < 0 for row in rows])),
            "optimized_psd_detection_rate": float(np.mean([row["optimized_psd_upper"] < 0 for row in rows])),
            "mean_optimized_rank": float(np.mean([row["optimized_psd_rank_1e-7"] for row in rows])),
            "distribution_free_coverage_rate": float(
                np.mean([row["inside_distribution_free_ellipsoid"] for row in rows])
            ),
            "distribution_free_rank_one_detection_rate": float(
                np.mean([row["distribution_free_rank_one_upper"] < 0 for row in rows])
            ),
            "distribution_free_psd_detection_rate": float(
                np.mean([row["distribution_free_psd_upper"] < 0 for row in rows])
            ),
        },
        "first_ten": rows[:10],
    }
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(payload["exact"], indent=2))
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
