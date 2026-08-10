# Proper-delay spectra majorize subspace rotation

This directory contains the sixth paper's source, frozen figure and final PDF.
It proves a sharp Ky Fan hierarchy connecting every prefix of a subspace's
canonical-angle spectrum to the matching prefix of the generator's spectral
spread and, for positive Wigner--Smith flows, to the leading proper delays.
The definitive version lifts this hierarchy to an exact variational principle
for every symmetric gauge: the minimum gauge action is `2 Phi(beta)`, attained
by one common positive coupled-mode path for all gauges simultaneously.

From the repository root, reproduce the computational certificate with:

```bash
python research/ky_fan_speed_limit_certificate.py
python verification/verify_ky_fan_speed_limits.py
python verification/run_ky_fan_speed_limit_release_checks.py
```

Compile with:

```bash
cd paper_ky_fan_speed_limits
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The proofs are analytic. The deterministic artifact audits sharp equality and
gauge-optimal families, a strict separation invisible to max and trace,
signed and positive matrix inequalities, noncommuting piecewise-constant
paths, noisy tomography and the exact slack identity.
