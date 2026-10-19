# The Action--Memory Diamond

This directory contains the eighth paper in the passive-network sequence. It
proves the exact joint attainable region for arc-integrated trace delay and
McMillan degree under a prescribed subspace transfer, gives an explicit
Blaschke--Potapov compiler, derives robust degree certificates, and exhibits a
strict three-node cyclic memory cost invisible to pairwise tests.

The manuscript also gives a theorem-level comparison with the preceding
modal-delay paper.  That work classifies the continuous delay spectrum and
proves minimum degree; this paper instead classifies every jointly feasible
total-action/degree pair, including the complementary-return face and the
discontinuous coincident-endpoint case.

From the repository root, regenerate and independently verify the artifact:

```powershell
python research/exact_action_memory_certificate.py
python verification/verify_exact_action_memory.py
```

Compile the manuscript from this directory:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```
