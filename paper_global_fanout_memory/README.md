# Orthogonal Spectral Fan-Out Costs k(L-1) States

This directory contains the ninth paper in the passive-network sequence.  It
proves that routing one rank-`k` input subspace to `L` mutually orthogonal
output subspaces at distinct boundary frequencies requires exactly
`k(L-1)` McMillan states, even though every two-node restriction requires
only `k`.

The paper also proves a general Pick--Stein negative-inertia lower bound, a
gauge-independent span bound, an operator-norm robust version, and a
fail-closed noisy certificate.  It compares the exact global cost with the
strictly weaker sum of endpoint-only Grassmann speed-limit bounds.  The full
Wigner--Smith trace action still recovers degree exactly.

From the repository root, regenerate and independently verify the artifact:

```powershell
python research/global_fanout_memory_certificate.py
python verification/verify_global_fanout_memory.py
```

Compile the manuscript from this directory:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```
