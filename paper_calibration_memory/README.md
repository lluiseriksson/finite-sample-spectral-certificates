# Lossless calibration is stored memory

This directory contains the manuscript source and generated figure for the
topological calibration-memory law.

Reproduce the numerical ledger and figure from the repository root:

```bash
python research/calibration_memory_certificate.py
python verification/verify_calibration_memory.py
```

Compile the paper with:

```bash
cd paper_calibration_memory
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The proofs of the main degree and delay identities are analytic. The JSON
ledger audits the explicit six-port construction, numerical quadrature, root
locations, and a passive perturbation governed by the Rouché criterion.
