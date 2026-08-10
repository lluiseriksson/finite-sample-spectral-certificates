# Finite-Sample Spectral Certificates

[![verify](https://github.com/lluiseriksson/finite-sample-spectral-certificates/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/lluiseriksson/finite-sample-spectral-certificates/actions/workflows/verify.yml)

Research artifact for finite-sample upper-support and visible-gap tests from
noisy block moments.  The central method combines a Wishart--Loewner confidence
band with truncated Hausdorff localizers.  The confidence set remains valid
after the separating matrix witness is optimized on the same data.

## New: topological calibration-memory law

The newest paper, `paper_calibration_memory/lossless_calibration_memory.pdf`,
proves an architecture-independent resource law for finite rational passive
networks. If a signal block is lossless on subspaces of total dimension `K`
across distinct boundary frequencies and is strictly contractive somewhere
else, then the full network has McMillan degree at least `K`. The same degree
is exactly its integrated Wigner--Smith trace delay. The coefficient is sharp;
an analytic Rouché condition makes the count perturbation-stable; and a
counterexample records why unstructured approximate samples alone cannot
support the claim. Applied to the irreducible six-port reservoir family, the
theorem proves universal near-optimality: `3S` calibrated directions require
at least `3S` states and the construction uses `3S+6`.

The active research branch also contains a separate new paper on
irreducible reciprocal MIMO FIR filters.  Its headline theorem gives a
fixed-degree exponential separation from every reducible three-channel
symmetric filter family, including noncommuting `1+2` block architectures and
therefore strictly more than fixed orthogonal banks of scalar filters.  A
sampled off-block defect makes the obstruction quantitative for approximately
reducible filters.  An exact five-tap witness and a second certificate grounded in the
public VBL-VA001 triaxial vibration data are replayed in standard Python.  An
FDD invariance lemma gives the directional calibrations a modal-preprocessing
meaning, and the held-out cospectral replay measures less than `0.222` degrees
of leading-direction drift.  Lean/mathlib kernel-checks the polynomial
root-count core used by the invariant-line obstruction.

The current frontier experiment studies Gaussian sketches of Euclidean
correlators in the interacting axial next-nearest-neighbour Ising (ANNNI)
chain.  Two finite-sample routes are implemented:

- a covariance-known, distribution-free Chebyshev ellipsoid;
- an unknown-covariance Wishart Loewner band followed by one semidefinite
  compatibility problem.

The main theorem gives unknown-covariance type-I control at finite sample size;
the companion power theorem gives an explicit witness-dependent sample bound.
Every numerical rejection uses a repaired dual witness and a conservative
stationarity-residual allowance.  One representative rejection is additionally
encoded and checked in exact rational arithmetic.

## Reproduce locally

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt
.venv/Scripts/python pilot_primal_dual_equivalence.py
.venv/Scripts/python pilot_wishart_loewner.py --length 8 --ensembles 10
python verification/verify_rational_witness.py
```

Linux/Colab users should replace `.venv/Scripts/python` with
`.venv/bin/python`, or use the supplied notebook.

## Current empirical boundary

The Hotelling ellipsoid is deliberately retained as a negative control: the
quadratic sketch features are Wishart-like, not Gaussian, and its nominal 95%
coverage under-covered at small sample sizes.  It must not be reported as an
exact certificate for this model.

The theorem is also deliberately not extended to arbitrary primitive laws.
The production misspecification grid shows severe undercoverage for a
variance-matched Student-$t_5$ primitive, while Rademacher sketches happen to
be conservative in the tested cells.  Neither observation changes the stated
Gaussian scope.

## Layout

- `THEOREM_DRAFT.md`: theorem statements and open audit points.
- `PILOT_LEDGER.md`: positive and negative experimental results.
- `pilot_wishart_loewner.py`: covariance-unknown finite-sample test.
- `production_annni_campaign.py`: covariance-known reference campaign.
- `colab/annni_wishart_campaign.ipynb`: scalable Colab entry point.
- `results/production/`: frozen Colab campaigns and their raw summaries.
- `results/verification/rational_witness.json`: exact-rational rejection witness.
- `verification/verify_rational_witness.py`: independent standard-library replay.
- `paper_noncommuting/`: new matrix-polynomial/MIMO FIR manuscript and PDF.
- `paper_quantum_reservoir/`: passive six-port decoherence-filter manuscript.
- `paper_calibration_memory/`: topological McMillan-degree/Wigner--Smith paper.
- `research/calibration_memory_certificate.py`: degree, delay and perturbation replay.
- `verification/verify_calibration_memory.py`: fail-closed new-paper verifier.
- `formal/`: Lean/mathlib proof of the scalar-branch root-count obstruction.
- `verification/verify_vbl_va001_witness.py`: exact measured-calibration replay.
- `data/route_c/vbl_va001_calibration.json`: public-data provenance and frozen split.
- `paper/main.tex`: full manuscript source; `output/pdf/main.pdf` is the built paper.

## What is and is not certified

A rejection proves that an asserted upper edge is incompatible with the
simultaneous confidence set, hence the asserted *visible* gap is too large.
Non-rejection is not a positive lower-gap certificate.  The finite ANNNI chain
is a controlled interacting benchmark, not a thermodynamic mass-gap proof.
