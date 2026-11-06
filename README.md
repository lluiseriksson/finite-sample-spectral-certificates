# Finite-Sample Spectral Certificates

<!-- RESEARCH-UPDATE-2026-09-09 -->
> **Research update · 9 September 2026:** [verified source cards, limits and next experiments](research/updates/2026-09-09.md) · [machine-readable dossier](research/updates/2026-09-09.json).

[![verify](https://github.com/lluiseriksson/finite-sample-spectral-certificates/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/lluiseriksson/finite-sample-spectral-certificates/actions/workflows/verify.yml)

Research artifact for finite-sample upper-support and visible-gap tests from
noisy block moments.  The central method combines a Wishart--Loewner confidence
band with truncated Hausdorff localizers.  The confidence set remains valid
after the separating matrix witness is optimized on the same data.

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
- `paper/main.tex`: full manuscript source; `output/pdf/main.pdf` is the built paper.

## What is and is not certified

A rejection proves that an asserted upper edge is incompatible with the
simultaneous confidence set, hence the asserted *visible* gap is too large.
Non-rejection is not a positive lower-gap certificate.  The finite ANNNI chain
is a controlled interacting benchmark, not a thermodynamic mass-gap proof.

## Optional fixed-monomial sign test

A separate [opt-in sign test](MONOMIAL_SIGN_TEST.md) handles one fixed
monomial/channel projection of raw 12D centred elliptical sketches with
finite covariance. It uses exact binomial tails and certified rational
intervals. Its declared-model and selection requirements are narrower than
the existing adaptive SDP workflow; it does not extend the Wishart theorem
or reinterpret prior misspecification results.

The separate [separated-time entrypoint](SEPARATED_SIGN_TEST.md) also permits
one preselected certified pair `[0,3]`, `[2,5]`, or `[0,5]`. Its exact SOS
proof applies to the full-covariance null without half-Hankel identities.
These options remain single tests, not an unadjusted scan across pairs, and
the resulting family does not characterize the full PSD localizer.

