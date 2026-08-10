# Complete delay-memory regions

This directory contains the seventh paper in the passive-network sequence.
It proves three linked results:

1. `2 beta ≺_w q` is necessary and sufficient for an integrated leading-delay
   profile `q` to transport two rank-`k` subspaces; every feasible profile has
   a compiler with at most `k+1` commuting positive segments.
2. The minimum McMillan degree for two-frequency subspace interpolation is
   exactly the number of nonzero canonical angles.
3. Boundary proper delays are diagonal blocks of a Pick memory-Gram matrix;
   its spectrum gives robust hidden-state and reduction-order certificates.

From the repository root, regenerate and verify the artifact with:

```powershell
python research/delay_memory_region_certificate.py
python verification/verify_delay_memory_regions.py
```

Compile the manuscript from this directory with:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The proofs are analytic. The deterministic code reconstructs the majorization
compiler, audits matrix Blaschke--Potapov Pick ranks, stress-tests the robust
thresholds for Pick matrices and noisy projectors, constructs sharp two-node
Blaschke interpolants, verifies the exact local-delay-blindness example, and
checks the transfer/proper-delay Fourier dictionary.
