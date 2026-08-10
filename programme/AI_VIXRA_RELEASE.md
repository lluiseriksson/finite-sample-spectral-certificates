# ai.viXra release record

This file freezes the exact preprint and the checks used for submission.  Values
marked `TO_BE_FROZEN` are replaced only after the final PDF has passed visual QA.

## Artifact

- Title: *Irreducible Channel Mixing and Exponential Calibration Laws for
  Matrix-Polynomial and Linear-Phase MIMO FIR Filters*
- Author: Lluis Eriksson
- Submission venue: ai.viXra.org
- Category: Computational Science — Digital Signal Processing
- Version: `v1.0-ai-vixra`
- Scientific PDF commit: `e62d06bf702d09b41a9c51937adf7508c94c65f8`
- PDF: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Pages: 18
- PDF SHA-256: `ad711a211c3457ba0fef30c827b5fa543403023cc485ab16b34d0ead393c0a81`

## One-command reproduction

From the repository root in the frozen environment:

```powershell
python verification/run_ai_vixra_release_checks.py
```

Expected result: nine `PASS` lines for the environment, six numerical or
exact-arithmetic replays, the Lean/mathlib build, the Lean source boundary, and
the PDF hash, followed by `PASS  all ai.viXra release gates`.

## Frozen environment

- Python 3.12.6
- NumPy 2.5.1
- SciPy 1.18.0
- CVXPY 1.9.2
- Lean `leanprover/lean4:v4.29.0-rc6`
- mathlib commit `07642720480157414db592fa85b626dafb71355b`
- TeX engine: pdfTeX 3.141592653-2.6-1.40.28 (MiKTeX 25.12)

The Python dependency lock is `requirements.txt`; Lean dependencies are frozen
by `formal/lean-toolchain` and `formal/lake-manifest.json`.

## Verification boundary

The scripts replay exact arithmetic, frozen experiments, and numerical controls.
Lean/mathlib kernel-checks the scalar polynomial root-count implication, without
`sorry` or user axioms.  The full asymptotic theorem remains a conventional
mathematical proof; the repository does not claim a complete formalization.

## Frozen source hashes

- `paper_noncommuting/main.tex` — `a331f215750ce1d783ae612b0d9a2c0a9345f2fc505505d672113c1ec6aad8cb`
- `verification/run_ai_vixra_release_checks.py` — `77de678445eef90992a6b5dddcef38d352d2e44789749a88da23e7ed1ab83e5d`
- `verification/verify_constant_channel_scaling.py` — `971c46b515e4c3b8ab32305fab1c1f4e81d4ce9c81e80620dd27c7a7c8e95412`
- `verification/verify_quadratic_filter_witness.py` — `11facbdd582dbde8cca8c5845e8de451cec7c3efe0129156cd6147862d2d5474`
- `verification/verify_vbl_va001_witness.py` — `a83995e9a4f95a8c096a6aec0796fc919d5f185cd7c0cc41518cebe45ba074f3`
- `research/route_c_block_krylov_benchmark.py` — `ac3c49efb10f5c5f6e512066bf220df9d7fa6a4b46f1a023f61dc6b06353c81f`
- `research/route_c_graph_filter_pilot.py` — `f7ca69d6272f960b8ff4677fcfdb6ea7002bc077f727ccdebac93ce1523803ac`
- `research/route_c_vibration_calibration.py` — `11308acb9b968bbba90adf2c1d7185749297f7c0748b6118404f92104f4d6531`
- `formal/SpectralCertificates/RootCount.lean` — `bd31700645631f10ee82c849e8b38b01dc4feddcfe8a3a0ecab4d72f7644e5c5`
- `formal/lake-manifest.json` — `5be13515668905c3b2c53c1b1a26fa6b3c5fd400ebaad3556adceafab63d2881`
- `data/route_c/vbl_va001_calibration.json` — `857eb351ccde567cdec6aa889acb8cc19e0ccbd7a3d9825472a8904462d36bec`
- `requirements.txt` — `c922db4e2debf4352bd9423cfc9053ca3805ca8e1a6bd5e7210886295febf15d`
