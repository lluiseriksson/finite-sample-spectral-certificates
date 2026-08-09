# Noncommuting-filter paper draft

This directory is a new paper, not version 3 of the finite-sample Gaussian
correlator paper.

## Artifact identity

- Title: *Exponential Tangential Advantages of Noncommuting Matrix-Polynomial
  and Linear-Phase MIMO FIR Filters*
- Source: `paper_noncommuting/main.tex`
- Compiled artifact: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Status: independently reviewed research release; the quantified proof gives an exponential
  separation at fixed dimension `d=3` with `M=N+3`; the exact five-tap FIR certificate,
  rational approximate-calibration margin, negative graph pilot, and an exact
  rational certificate calibrated from held-out-controlled VBL-VA001 sensor
  records are closed.  The modal-task interpretation is backed by an exact FDD
  invariance lemma and a held-out cospectral replay.  A fresh independent
  adversarial review of the exact hash below assigned 7/10 to both scientific
  contribution and manuscript quality and found no proof defect.  The
  restricted full-text priority uncertainty remains disclosed; this release
  is not a certificate of novelty or named-human peer review.

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

The measured calibration replays offline from the frozen summary:

```powershell
python verification/verify_vbl_va001_witness.py
python research/route_c_vibration_calibration.py
```

To reproduce the summary from the public archive, run
`research/acquire_vbl_va001_subset.py` followed by
`research/extract_vbl_va001_calibration.py`.  Range requests retrieve only the
twelve selected CSV records, whose hashes are fixed in the manifest.
