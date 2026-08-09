# Finite-Sample Spectral Certificates

Research code for rigorous upper-support and gap tests from noisy block moments.
The central method combines truncated Hausdorff localizers with confidence sets
that remain valid after the separating witness is optimized on the same data.

The current frontier experiment studies Gaussian sketches of Euclidean
correlators in the interacting axial next-nearest-neighbour Ising (ANNNI)
chain.  Two finite-sample routes are implemented:

- a covariance-known, distribution-free Chebyshev ellipsoid;
- an unknown-covariance Wishart Loewner band followed by one semidefinite
  compatibility problem.

This repository is an active research artifact.  Claims in `THEOREM_DRAFT.md`
are explicitly marked as under audit until the manuscript and independent
numerical certificates are complete.

## Reproduce locally

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt
.venv/Scripts/python pilot_primal_dual_equivalence.py
.venv/Scripts/python pilot_wishart_loewner.py --length 8 --ensembles 10
```

Linux/Colab users should replace `.venv/Scripts/python` with
`.venv/bin/python`, or use the supplied notebook.

## Current empirical boundary

The Hotelling ellipsoid is deliberately retained as a negative control: the
quadratic sketch features are Wishart-like, not Gaussian, and its nominal 95%
coverage under-covered at small sample sizes.  It must not be reported as an
exact certificate for this model.

## Layout

- `THEOREM_DRAFT.md`: theorem statements and open audit points.
- `PILOT_LEDGER.md`: positive and negative experimental results.
- `pilot_wishart_loewner.py`: covariance-unknown finite-sample test.
- `production_annni_campaign.py`: covariance-known reference campaign.
- `colab/annni_wishart_campaign.ipynb`: scalable Colab entry point.

