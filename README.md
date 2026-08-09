# Finite-Sample Spectral Certificates

[![verify](https://github.com/lluiseriksson/finite-sample-spectral-certificates/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/lluiseriksson/finite-sample-spectral-certificates/actions/workflows/verify.yml)

Research artifact for finite-sample upper-support and visible-gap tests from
noisy block moments.  The central method combines a Wishart--Loewner confidence
band with truncated Hausdorff localizers.  The confidence set remains valid
after the separating matrix witness is optimized on the same data.

The active research branch also contains a separate new paper on
noncommuting reciprocal MIMO FIR filters.  Its headline theorem gives a
fixed-degree exponential separation from every fixed orthogonal bank of scalar
filters; an exact five-tap witness and a second certificate grounded in the
public VBL-VA001 triaxial vibration data are replayed in standard Python.

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
- `verification/verify_vbl_va001_witness.py`: exact measured-calibration replay.
- `data/route_c/vbl_va001_calibration.json`: public-data provenance and frozen split.
- `paper/main.tex`: full manuscript source; `output/pdf/main.pdf` is the built paper.

## What is and is not certified

A rejection proves that an asserted upper edge is incompatible with the
simultaneous confidence set, hence the asserted *visible* gap is too large.
Non-rejection is not a positive lower-gap certificate.  The finite ANNNI chain
is a controlled interacting benchmark, not a thermodynamic mass-gap proof.
