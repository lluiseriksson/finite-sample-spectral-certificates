"""Pilot B: covariance-optimal finite-moment support falsifiers.

The test uses a known Gaussian confidence ellipsoid for moments m_1,...,m_K.
For a PSD lifted polynomial mixture X, its exact localizer is tr(A(m) X).
The support function of the ellipsoid gives an exact simultaneous margin,
including after X is selected from the same data by convex optimization.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cvxpy as cp
import numpy as np
from scipy.stats import chi2, norm


def moments(atoms: np.ndarray, weights: np.ndarray, count: int) -> np.ndarray:
    return np.array([np.sum(weights * atoms**k) for k in range(count)])


def basis_power_coefficients(degree: int) -> np.ndarray:
    """Rows are shifted Chebyshev basis polynomials on [0,1]."""
    from numpy.polynomial import Chebyshev, Polynomial

    answer = np.zeros((degree + 1, degree + 1))
    for j in range(degree + 1):
        coefficients = Chebyshev.basis(j, domain=[0.0, 1.0]).convert(kind=Polynomial).coef
        answer[j, : len(coefficients)] = coefficients
    return answer


def localizer_design(theta: float, degree: int) -> np.ndarray:
    """D[k] gives the coefficient multiplying power moment m_k."""
    basis = basis_power_coefficients(degree)
    design = np.zeros((2 * degree + 2, degree + 1, degree + 1))
    for i in range(degree + 1):
        for j in range(degree + 1):
            product = np.convolve(basis[i], basis[j])
            design[: len(product), i, j] += theta * product
            design[1 : len(product) + 1, i, j] -= product
    return design


def localizer_matrix(moment: np.ndarray, theta: float, degree: int) -> np.ndarray:
    return np.einsum("k,kij->ij", moment, localizer_design(theta, degree))


def chebyshev_coefficients(theta: float, target: float, degree: int) -> np.ndarray:
    from numpy.polynomial import Chebyshev, Polynomial

    polynomial = Chebyshev.basis(degree, domain=[0.0, theta]).convert(kind=Polynomial)
    power = np.pad(
        polynomial.coef / polynomial(target), (0, degree + 1 - len(polynomial.coef))
    )
    # Convert power coefficients to the shifted-Chebyshev optimization basis.
    return np.linalg.solve(basis_power_coefficients(degree).T, power)


def margins(
    observed: np.ndarray,
    covariance: np.ndarray,
    theta: float,
    target: float,
    degree: int,
    coverage: float,
    ellipsoid_q: float | None = None,
    box_multiplier: float | None = None,
) -> dict[str, float | int | list[float]]:
    all_design = localizer_design(theta, degree)
    design = all_design[1:]
    q = (
        float(np.sqrt(chi2.ppf(coverage, covariance.shape[0])))
        if ellipsoid_q is None
        else ellipsoid_q
    )
    covariance_root = np.linalg.cholesky(covariance)
    a = chebyshev_coefficients(theta, target, degree)
    fixed_x = np.outer(a, a)
    fixed_l = np.einsum("kij,ij->k", design, fixed_x)
    fixed_value = float(np.trace(localizer_matrix(observed, theta, degree) @ fixed_x))
    ellipsoid_fixed = fixed_value + q * float(np.linalg.norm(covariance_root.T @ fixed_l))

    alpha = 1.0 - coverage
    box_z = (
        float(norm.ppf(1.0 - alpha / (2.0 * covariance.shape[0])))
        if box_multiplier is None
        else box_multiplier
    )
    box_fixed = fixed_value + box_z * float(
        np.dot(np.sqrt(np.diag(covariance)), np.abs(fixed_l))
    )

    x = cp.Variable((degree + 1, degree + 1), PSD=True)
    matrix = localizer_matrix(observed, theta, degree)
    l_expression = cp.hstack([cp.trace(design[k] @ x) for k in range(design.shape[0])])
    objective = cp.trace(matrix @ x) + q * cp.norm(covariance_root.T @ l_expression, 2)
    # The robust sign objective is positively homogeneous.  Trace one is a
    # compact slice meeting every nonzero PSD ray, so it changes neither the
    # existence nor nonexistence of a negative certificate.
    problem = cp.Problem(cp.Minimize(objective), [cp.trace(x) == 1.0])
    solver_used = "CLARABEL-tight"
    try:
        problem.solve(solver="CLARABEL", tol_gap_abs=1e-9, tol_feas=1e-9, max_iter=500)
    except cp.error.SolverError:
        solver_used = "CLARABEL-default"
        try:
            problem.solve(solver="CLARABEL")
        except cp.error.SolverError:
            solver_used = "SCS"
            problem.solve(solver="SCS", eps=1e-7, max_iters=50_000)
    if problem.status not in {"optimal", "optimal_inaccurate"}:
        raise RuntimeError(f"solver status {problem.status}")
    raw_x = (x.value + x.value.T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(raw_x)
    # Convert the numerical solver output into an explicitly PSD, normalized
    # matrix before evaluating the certificate.
    x_value = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T
    normalization = float(np.trace(x_value))
    x_value /= normalization
    eigenvalues = np.linalg.eigvalsh(x_value)
    optimized_l = np.einsum("kij,ij->k", design, x_value)
    optimized_measured = float(np.trace(matrix @ x_value))
    optimized_margin = q * float(np.linalg.norm(covariance_root.T @ optimized_l))
    return {
        "fixed_chebyshev_ellipsoid_upper": ellipsoid_fixed,
        "fixed_chebyshev_box_upper": box_fixed,
        "optimized_ellipsoid_upper": optimized_measured + optimized_margin,
        "optimized_measured_localizer": optimized_measured,
        "optimized_ellipsoid_margin": optimized_margin,
        "optimized_recomputed_upper": optimized_measured + optimized_margin,
        "optimized_trace_normalization": float(np.trace(x_value)),
        "solver_used": solver_used,
        "optimized_lift_rank_1e-7": int(np.count_nonzero(eigenvalues > 1e-7)),
        "optimized_lift_eigenvalues": eigenvalues.tolist(),
        "optimized_lift": x_value.tolist(),
    }


def covariance_model(count: int, sigma: float, correlation: float) -> np.ndarray:
    indices = np.arange(count)
    correlation_matrix = correlation ** np.abs(indices[:, None] - indices[None, :])
    scales = sigma * (1.0 + 0.08 * indices)
    return correlation_matrix * scales[:, None] * scales[None, :]


def trial_campaign(
    rng: np.random.Generator,
    true_moments: np.ndarray,
    covariance: np.ndarray,
    theta: float,
    target: float,
    degree: int,
    coverage: float,
    trials: int,
) -> dict[str, object]:
    detections = {"fixed_ellipsoid": 0, "fixed_box": 0, "optimized_ellipsoid": 0}
    ellipsoid_misses = 0
    optimized_detections_inside_ellipsoid = 0
    examples = []
    for trial in range(trials):
        observed = true_moments.copy()
        observed[1:] += rng.multivariate_normal(np.zeros(covariance.shape[0]), covariance)
        result = margins(observed, covariance, theta, target, degree, coverage)
        error = observed[1:] - true_moments[1:]
        mahalanobis_sq = float(error @ np.linalg.solve(covariance, error))
        ellipsoid_q_sq = float(chi2.ppf(coverage, covariance.shape[0]))
        inside = mahalanobis_sq <= ellipsoid_q_sq
        detections["fixed_ellipsoid"] += result["fixed_chebyshev_ellipsoid_upper"] < 0
        detections["fixed_box"] += result["fixed_chebyshev_box_upper"] < 0
        detections["optimized_ellipsoid"] += result["optimized_ellipsoid_upper"] < 0
        ellipsoid_misses += not inside
        optimized_detections_inside_ellipsoid += inside and result["optimized_ellipsoid_upper"] < -1e-7
        if trial < 3:
            x_value = np.asarray(result["optimized_lift"])
            result["true_localizer_at_optimized_lift"] = float(
                np.trace(localizer_matrix(true_moments, theta, degree) @ x_value)
            )
            result["mahalanobis_sq"] = mahalanobis_sq
            result["inside_confidence_ellipsoid"] = inside
            examples.append(result)
    return {
        "trials": trials,
        "detections": detections,
        "rates": {key: value / trials for key, value in detections.items()},
        "ellipsoid_misses": ellipsoid_misses,
        "ellipsoid_miss_rate": ellipsoid_misses / trials,
        "optimized_detections_inside_ellipsoid": optimized_detections_inside_ellipsoid,
        "first_three": examples,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=5)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--sigma", type=float, default=2e-5)
    parser.add_argument("--correlation", type=float, default=0.92)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("pilot_covariance_sdp.json"))
    args = parser.parse_args()

    theta, delta, gamma = 0.8, 0.1, 0.01
    low_atoms = np.array([0.15, 0.45, 0.78])
    low_weights = np.array([0.2, 0.3, 0.5])
    # The boundary null delta_theta has zero localizer in every direction, so
    # it stress-tests simultaneous coverage rather than benefiting from slack.
    null_moments = moments(np.array([theta]), np.array([1.0]), 2 * args.degree + 2)
    alternative_moments = moments(
        np.append(low_atoms, theta + delta),
        np.append((1.0 - gamma) * low_weights, gamma),
        2 * args.degree + 2,
    )
    covariance = covariance_model(2 * args.degree + 1, args.sigma, args.correlation)
    rng = np.random.default_rng(args.seed)
    payload = {
        "setup": {
            "theta": theta,
            "delta": delta,
            "gamma": gamma,
            "degree": args.degree,
            "coverage": args.coverage,
            "sigma": args.sigma,
            "ar1_correlation": args.correlation,
            "seed": args.seed,
        },
        "null_support_inside_claim": trial_campaign(
            rng, null_moments, covariance, theta, theta + delta, args.degree, args.coverage, args.trials
        ),
        "alternative_visible_atom": trial_campaign(
            rng,
            alternative_moments,
            covariance,
            theta,
            theta + delta,
            args.degree,
            args.coverage,
            args.trials,
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["null_support_inside_claim"]["rates"], indent=2))
    print(json.dumps(payload["alternative_visible_atom"]["rates"], indent=2))


if __name__ == "__main__":
    main()
