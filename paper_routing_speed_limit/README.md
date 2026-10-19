# Every spectral switch costs memory

This directory contains the source, frozen figure and final PDF for the sharp
robust Wigner--Smith routing speed limit.

Reproduce the certificate from the repository root:

```bash
python research/routing_speed_limit_certificate.py
python verification/verify_routing_speed_limit.py
python verification/run_routing_speed_limit_release_checks.py
```

Compile the manuscript with:

```bash
cd paper_routing_speed_limit
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The principal theorems are analytic.  The numerical artifact checks explicit
finite-error equality cases, exact cyclic saturation, the two positive-block
inequalities and independent random Blaschke--Potapov paths.
