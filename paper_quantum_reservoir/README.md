# Passive quantum reservoir architecture paper

This directory contains a new paper connecting architecture-dependent passive
filtering to decoherence-rate and coherence-maintenance boundaries.

## Core claim

At fixed rational bidegree and under the same full-spark delayed tangential
calibrations, an explicitly realizable irreducible three-signal/three-loss-port
passive network achieves exponentially small stopband leakage. Every transfer
with a constant nontrivial reducing channel retains unit leakage. A uniformly
nondegenerate bath spectral density squares the gap at the Kossakowski-rate
level. The exact obstruction has a quantitative finite-error version, and an
explicit scalar Schur comparator proves that its conditioning must deteriorate
at least exponentially. A closed Markov Ramsey model accounts for every vacuum
port as a common measurable floor. A separate theorem proves that exponential
suppression requires large all-pass group delay, so the resource is moved into
passive memory.

## Reproduce the results

From the repository root:

```powershell
python research/quantum_reservoir_filter.py
python research/quantum_reservoir_robustness.py
python research/quantum_reservoir_robust_bound.py
python research/quantum_reservoir_closed_dephasing.py
python verification/verify_passive_quantum_filter.py
```

The first command regenerates:

- `results/quantum_reservoir/passive_filter_scaling.json`
- `paper_quantum_reservoir/figures/passive_scaling.pdf`
- `paper_quantum_reservoir/figures/passive_architecture.pdf`

The tolerance command regenerates the frozen passive-preserving component
stress test and `paper_quantum_reservoir/figures/passive_tolerance.pdf`.

The multiprecision command writes
`results/quantum_reservoir/robust_separation.json`; the closed-model command
writes `results/quantum_reservoir/closed_dephasing.json` and regenerates
`paper_quantum_reservoir/figures/closed_dephasing.pdf`.

The verification command independently checks the exact quadratic no-go, pass-node
calibrations, full-spark signatures, global Schur bound, stop envelope,
six-port boundary unitarity, rational root count, approximate-separation
ingredients, scalar robustness barrier, vacuum covariance closure, squared
rate law and delay identity.

## Compile the manuscript

```powershell
New-Item -ItemType Directory -Force paper_quantum_reservoir/build | Out-Null
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper_quantum_reservoir/build paper_quantum_reservoir/main.tex
bibtex paper_quantum_reservoir/build/main
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper_quantum_reservoir/build paper_quantum_reservoir/main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper_quantum_reservoir/build paper_quantum_reservoir/main.tex
```

The release PDF is copied only after all references resolve, all verification
gates pass and every rendered page has been inspected.

After the artifact manifest has been frozen, the top-level release gate is:

```powershell
python verification/run_quantum_reservoir_release_checks.py
```
