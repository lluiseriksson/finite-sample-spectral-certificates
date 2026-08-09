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

## Candidate next paper (unfrozen research draft)

- File in research branch: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Title: *Exponential Tangential Advantages of Noncommuting Matrix-Polynomial
  and Linear-Phase MIMO FIR Filters*
- Pages: 12
- Local SHA-256: `d5a79b1f5a0ff627ebf052929fe0167b907b34240b702df6aa20bde49446bba8`
- Status: **unfrozen research draft**.  It contains the quantified Week-3
  theorem proof, exact affine and full-spark quadratic witnesses, an exact
  palindromic five-tap FIR realization, a rational commuting lower bound under
  approximate calibration, the block-Krylov proposition, a
  controlled Hermitian subspace benchmark, and a deliberately negative natural
  graph pilot.  It now also contains a public VBL-VA001 triaxial calibration:
  exact rational arithmetic proves continuum leakage below `24/25` versus one
  for the exactly calibrated commuting class, and six held-out records give
  maximum directional residual `0.0038491`.  Gate G1 remains open pending the
  full-text/FIR priority audit and a task-level validation that mandates the
  directional frequency calibrations and fixed delay.

This draft must not be scored as a final 7+ submission.  If it is reviewed,
the evaluator must quote the title and hash above so the feedback can be
attributed to the correct bytes.

The frozen next-paper entry will be created only after Gates G1--G3 are passed.
A release candidate must have:

1. a unique title that does not resemble either prior paper;
2. an embedded version string and Git commit;
3. a SHA-256 generated from the public GitHub bytes, not a platform-normalized
   local copy;
4. a one-page evaluator cover sheet listing the headline theorem, closest
   literature, what is new, and what is explicitly not claimed.
