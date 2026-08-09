# Artifact identity manifest

This manifest is the first gate of the 7+ programme.  No score is accepted
unless the evaluator reports the exact SHA-256 of the file it read.

## Superseded paper (the 5.31 evaluation)

- File: `finite_window_gap_identifiability.pdf`
- Title: *Noisy Euclidean Correlators Do Not Certify a Spectral Gap Without Visibility*
- Pages: 11
- SHA-256: `7baeb546e47fc8e004c32809569c3e39aa3089fc757c7f514edb44ab474eaa49`
- Status: scientifically useful prior work, but **not** the artifact to use when
  evaluating the later finite-sample Wishart paper or the new 7+ programme.

The assessment stored at
`C:/Users/lluis/.codex/attachments/e2c4ccac-e6d3-4b60-9ee8-4e675617a239/pasted-text.txt`
explicitly reports this hash and is therefore an assessment of the superseded
paper.

## Current published statistical paper

- Repository: <https://github.com/lluiseriksson/finite-sample-spectral-certificates>
- File in repository: `output/pdf/main.pdf`
- Title: *Finite-Sample Spectral-Edge Exclusion from Gaussian Correlator Sketches*
- Pages: 12
- SHA-256: `593b5b5f7d568bc546bd2f18c2cd166df9d6c4bd402888bc2075183e75ca9c7b`
- Status: published infrastructure and baseline for the new programme; it is
  **not presumed** to be a 7+ paper.

## Candidate next paper (independently reviewed robustness-law release)

- File in research branch: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Title: *Exponential Calibration Laws for Noncommuting Matrix-Polynomial and
  Linear-Phase MIMO FIR Filters*
- Pages: 15
- Local SHA-256: `d62ed80dba4a4c9cd42016aa23931277154a6bdfc070686ff19b91e8f3e7e323`
- Scientific commit: `46ed63d12d1639b949625f67b1dfbc175ad4580f`
- Public-byte SHA-256 recheck: identical to the local hash above.
- Status: **independently reviewed scientific revision**.
  Its sharpened theorem needs only
  the fixed dimension `d=3` and `M=N+3` signatures, rather than a growing
  channel dimension; it also contains
  a quantitative rational full-spark margin, an exponential robust-calibration
  lower bound for the commuting class, and a universal scalar construction
  proving that an exponential tolerance scale is necessary.  It also contains
  exact affine and full-spark quadratic witnesses, an exact
  palindromic five-tap FIR realization, a rational commuting lower bound under
  approximate calibration, the block-Krylov proposition, a
  controlled Hermitian subspace benchmark, and a deliberately negative natural
  graph pilot.  It now also contains a public VBL-VA001 triaxial calibration:
  exact rational arithmetic proves continuum leakage below `24/25` versus one
  for the exactly calibrated commuting class, and six held-out records give
  maximum directional residual `0.0038491`.  An exact FDD modal-component
  invariance lemma gives the directional constraints a task-level meaning; on
  frozen held-out cospectra the leading modal direction moves by at most
  `0.2217` degrees and the leading ordinate changes by at most `2.40e-5`
  relatively.  A fresh blind Gemini 3.1 Pro audit of the exact attached PDF
  reported the title and 15 pages, marked all eight targeted proof checks
  `PASS`, scored the scientific contribution `9.0/10` and manuscript quality
  `8.5/10`, and recommended human peer review.  It did not hash the bytes or run
  scripts; local/public hashing and GitHub CI are separate evidence.  The
  restricted full-text/FIR priority uncertainty remains open and
  explicitly prevents an unqualified novelty claim.  The audit also
  treats Kootsookos's 1991 fixed-length MIMO `H-infinity` thesis and
  Alpay--Lewkowicz's structured matrix-polynomial interpolation as direct
  collisions with broader formulations of the claim.

The scores above apply only to the current title, page count and hash.  Any PDF
edit creates a new artifact and invalidates the review gate.

The release candidate has:

1. a unique title that does not resemble either prior paper;
2. an embedded date and a separately recorded scientific commit;
3. a SHA-256 checked locally and again from the public GitHub bytes;
4. a one-page evaluator cover sheet listing the headline theorem, closest
   literature, what is new, and what is explicitly not claimed.
