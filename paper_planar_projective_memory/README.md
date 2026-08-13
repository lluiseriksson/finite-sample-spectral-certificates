# Projective bottlenecks in passive spectral routing

This directory contains the source of the paper proving that passive memory on
the coplanar line-target stratum equals the minimum projective rational
interpolation degree, and deriving a singular-value-controlled lower bound on
the worst-node chordal error of every low-state lossless router.  It also solves
the complete zero-memory minimax frontier through smallest caps on the Bloch
sphere and proves sharp collision asymptotics for two targets.

Build:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Exact artifact:

```powershell
python ..\research\planar_projective_memory_certificate.py `
  --output ..\results\planar_projective_memory\certificate.json `
  --figure figures\unbounded_detector_gap.pdf --max-L 16
python ..\verification\verify_planar_projective_memory.py
```
