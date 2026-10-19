# Irreducible-channel-filter paper draft

This directory is a new paper, not version 3 of the finite-sample Gaussian
correlator paper.

## Artifact identity

- Title: *Irreducible Channel Mixing and Exponential Calibration Laws for
  Matrix-Polynomial and Linear-Phase MIMO FIR Filters*
- Source: `paper_noncommuting/main.tex`
- Compiled artifact: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Status: ai.viXra release candidate; the quantified proof gives an exponential
  separation at fixed dimension `d=3` with `M=N+3` from every reducible
  symmetric coefficient family, including genuinely noncommuting `1+2`
  block architectures.  A sampled off-block defect quantifies approximate
  reducibility.  The proof replaces the old
  `exp(-N^3)` perturbation by rational `2^(-13N)` signatures, and proves a
  universal scalar barrier showing that constant calibration tolerance is
  impossible.  The new exact five-tap FIR certificate has leakage `25/32` and
  preserves strict separation for `delta < 7/1920`.  The negative graph pilot and an exact
  rational certificate calibrated from held-out-controlled VBL-VA001 sensor
  records are closed.  The modal-task interpretation is backed by an exact FDD
  invariance lemma and a held-out cospectral replay.  Reviews of the superseded
  15- and 16-page PDFs do not apply to this 18-page irreducibility release; fresh
  independent review is pending.  Restricted full-text priority uncertainty
  remains disclosed.  Lean/mathlib kernel-checks the scalar
  polynomial root-count core of the invariant-line obstruction.

The current PDF hash is recorded after each accepted compilation in
`programme/ARTIFACT_MANIFEST.md`.  Evaluations must quote the title and SHA-256
from that manifest; a score for any other PDF is not evidence about this draft.

## Reproduce

The frozen release has one top-level gate:

```powershell
python verification/run_ai_vixra_release_checks.py
```

It checks the pinned Python environment, all paper-facing replays, the
Lean/mathlib build and source boundary, and the exact PDF SHA-256.  See
`programme/AI_VIXRA_RELEASE.md` for the immutable release record.

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

The formal core is replayed separately:

```powershell
Set-Location formal
lake exe cache get
lake build
```
