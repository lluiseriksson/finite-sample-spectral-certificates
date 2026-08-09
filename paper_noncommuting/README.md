# Noncommuting-filter paper draft

This directory is a new paper, not version 3 of the finite-sample Gaussian
correlator paper.

## Artifact identity

- Title: *Noncommuting Hermitian Matrix-Polynomial Filters Can Have an
  Exponential Tangential Advantage*
- Source: `paper_noncommuting/main.tex`
- Compiled artifact: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Status: research draft; the quantified proof is closed, while the priority
  audit and natural-application gate remain open under Gate G1.

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
