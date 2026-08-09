# Noncommuting-filter paper draft

This directory is a new paper, not version 3 of the finite-sample Gaussian
correlator paper.

## Artifact identity

- Title: *Exponential Tangential Advantages of Noncommuting Matrix-Polynomial
  and Linear-Phase MIMO FIR Filters*
- Source: `paper_noncommuting/main.tex`
- Compiled artifact: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Status: research draft; the quantified proof, exact five-tap FIR certificate
  and negative graph pilot are closed, while the full-text priority audit and
  a measured/standards-derived application remain open under Gate G1.

The current PDF hash is recorded after each accepted compilation in
`programme/ARTIFACT_MANIFEST.md`.  Evaluations must quote the title and SHA-256
from that manifest; a score for any other PDF is not evidence about this draft.

## Reproduce

From this directory on a TeX installation with `pdflatex` and `bibtex`:

```powershell
New-Item -ItemType Directory -Force build | Out-Null
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
bibtex build/main
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
```

Copy `build/main.pdf` to `noncommuting_filters_draft.pdf` only after the log has
no unresolved references and all rendered pages have been visually inspected.
