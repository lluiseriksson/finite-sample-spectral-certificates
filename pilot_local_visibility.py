"""Pilot A: optimal local visibility of the lowest excitation.

For a contiguous block A and two normalized states, form

    X_A = Tr_{A^c} |psi0><psi1|.

The Frobenius norm is the aggregate transition strength over an orthonormal
operator basis; the nuclear norm is the exact optimum over all A-local probes
with operator norm at most one.  This pilot checks those identities and maps
their finite-size behaviour in TFIM/ANNNI regimes before any paper claim is
made.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh


REGIMES = {
    "tfim_paramagnet": {"j1": 1.0, "j2": 0.0, "hx": 1.5},
    "tfim_critical": {"j1": 1.0, "j2": 0.0, "hx": 1.0},
    "tfim_ferromagnet": {"j1": 1.0, "j2": 0.0, "hx": 0.6},
    "annni_paper": {"j1": 1.0, "j2": 0.37, "hx": 2.2},
    "annni_frustrated": {"j1": 1.0, "j2": 0.65, "hx": 1.0},
}


def hamiltonian(length: int, j1: float, j2: float, hx: float) -> sparse.csr_matrix:
    dim = 1 << length
    states = np.arange(dim, dtype=np.int64)
    diagonal = np.zeros(dim)
    z = [1.0 - 2.0 * ((states >> site) & 1) for site in range(length)]
    for site in range(length - 1):
        diagonal -= j1 * z[site] * z[site + 1]
    for site in range(length - 2):
        diagonal -= j2 * z[site] * z[site + 2]
    rows = [states]
    cols = [states]
    data = [diagonal]
    for site in range(length):
        rows.append(states)
        cols.append(states ^ (1 << site))
        data.append(np.full(dim, -hx))
    return sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(dim, dim),
    ).tocsr()


def parity(vector: np.ndarray) -> float:
    dim = vector.size
    complement = (dim - 1) ^ np.arange(dim)
    return float(np.vdot(vector, vector[complement]).real)


def transition_reduction(
    psi0: np.ndarray, psi1: np.ndarray, length: int, sites: list[int]
) -> np.ndarray:
    complement = [site for site in range(length) if site not in sites]
    order = sites + complement
    # Fortran order makes tensor axis `site` correspond to bit `site`.
    shape = [2] * length
    ground = psi0.reshape(shape, order="F").transpose(order)
    excited = psi1.reshape(shape, order="F").transpose(order)
    block_dim = 1 << len(sites)
    ground = ground.reshape(block_dim, -1)
    excited = excited.reshape(block_dim, -1)
    return ground @ excited.conj().T


def centered_sites(length: int, block_size: int) -> list[int]:
    start = (length - block_size) // 2
    return list(range(start, start + block_size))


def single_z_weight(
    psi0: np.ndarray, psi1: np.ndarray, length: int, site: int
) -> float:
    states = np.arange(1 << length)
    z = 1.0 - 2.0 * ((states >> site) & 1)
    return float(abs(np.vdot(psi1, z * psi0)) ** 2)


def run_case(length: int, regime: str, max_block: int) -> dict[str, object]:
    parameters = REGIMES[regime]
    h = hamiltonian(length, **parameters)
    energies, vectors = eigsh(h, k=3, which="SA", tol=1e-11, maxiter=300_000)
    indices = np.argsort(energies)
    energies = energies[indices]
    vectors = vectors[:, indices]
    psi0, psi1 = vectors[:, 0], vectors[:, 1]
    rows = []
    for block_size in range(1, min(max_block, length) + 1):
        sites = centered_sites(length, block_size)
        transition = transition_reduction(psi0, psi1, length, sites)
        singular_values = np.linalg.svd(transition, compute_uv=False)
        frobenius_sq = float(np.sum(singular_values**2))
        nuclear_sq = float(np.sum(singular_values) ** 2)
        rows.append(
            {
                "block_size": block_size,
                "sites": sites,
                "aggregate_hs_visibility": frobenius_sq,
                "optimal_contraction_visibility": nuclear_sq,
                "transition_rank": int(np.count_nonzero(singular_values > 1e-12)),
                "largest_singular_value": float(singular_values[0]),
            }
        )
    return {
        "regime": regime,
        "parameters": parameters,
        "length": length,
        "energies": energies.tolist(),
        "gap": float(energies[1] - energies[0]),
        "parities": [parity(vectors[:, i]) for i in range(3)],
        "center_z_weight": single_z_weight(psi0, psi1, length, length // 2),
        "blocks": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", nargs="+", type=int, default=[6, 8, 10, 12])
    parser.add_argument("--regimes", nargs="+", choices=REGIMES, default=list(REGIMES))
    parser.add_argument("--max-block", type=int, default=6)
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("pilot_local_visibility.json")
    )
    args = parser.parse_args()
    payload = {
        "definition": "X_A = Tr_{A^c} |psi0><psi1|",
        "normalizations": {
            "aggregate_hs_visibility": "||X_A||_F^2",
            "optimal_contraction_visibility": "||X_A||_*^2",
        },
        "results": [
            run_case(length, regime, args.max_block)
            for regime in args.regimes
            for length in args.lengths
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(payload['results'])} cases to {args.output}")


if __name__ == "__main__":
    main()
