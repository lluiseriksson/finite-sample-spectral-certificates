# Adversarial evaluator packet: irreducibility and calibration-law manuscript

This cover sheet identifies the exact artifact to inspect and asks for a blind,
claim-level assessment.  It is not a certificate of novelty and it is not a
request to reward repository size.

## Artifact identity

- Repository: <https://github.com/lluiseriksson/finite-sample-spectral-certificates>
- Research branch: `research/endogenous-visibility`
- Draft pull request: <https://github.com/lluiseriksson/finite-sample-spectral-certificates/pull/2>
- Scientific source commit containing the PDF bytes below:
  `cb1c8f4b9267cf55e586fcb17db8ae3a3115d4a2`.
- PDF path: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Title visible on page 1: *Irreducible Channel Mixing and Exponential
  Calibration Laws for Matrix-Polynomial and Linear-Phase MIMO FIR Filters*
- Page count: 16
- SHA-256 of the repository PDF bytes:
  `1de63a470acc5ee8e02268b9b73779e474c25b45c05139a6108c239946f9a6ee`
- Embedded manuscript date: `10 August 2026`
- CI status at that commit: pending.  The workflow now has separate numerical
  replay and Lean/mathlib kernel-check jobs.
- Independent review: pending for these exact bytes.  Reviews of the preceding
  15-page PDF are superseded and must not be reused.

Before scoring, independently hash the downloaded PDF and copy the result into
the response.  If the hash, title or page count differs, stop: the score belongs
to another artifact.  The release hash was checked both locally and from the
public commit URL.  If a later commit changes the PDF, this packet must be
regenerated rather than silently reused.

## Claims actually submitted for assessment

1. **Theorem 2.1:** for every sufficiently large degree `N`, already at fixed
   dimension `d=3` and `M=N+3`, there are full-spark tangential pass
   constraints for which an irreducible real-symmetric matrix polynomial
   attains stopband leakage `C exp(-cN)`, whereas every feasible reducible
   real-symmetric coefficient family has leakage at least one.  This comparator
   strictly includes noncommuting `1+2` block architectures; the commuting
   subcase is identically the identity.  Quantitatively, the signatures have
   `gamma_N >= 2^(-39N)/(48(N+4)^6)` and every approximately calibrated
   symmetric competitor has leakage at least
   `1-sqrt(3)(24e/5)^N(delta+beta_u)/gamma_N`, where `beta_u` is its sampled
   coupling away from a candidate invariant channel line.
2. **Corollary 4.1:** the same construction gives a fixed-latency reciprocal
   linear-phase MIMO FIR separation between irreducible taps and every
   reducible symmetric architecture, including the robust near-reducibility
   law.
3. **Proposition 2.7:** for arbitrary pass nodes and directions in the same
   separated bands, a scalar polynomial has both pass residual and stop leakage
   at most `exp(-2N/81)`.  Consequently no universal separation can survive a
   constant calibration tolerance; the exponential tolerance scale is
   necessary as a scale class.
4. **Propositions 5.1--5.2:** an exactly replayable rational five-tap witness
   has continuum leakage at most `25/32`, and every comparator satisfies the
   finite lower bound `1-60(delta+beta_u)`; the reducible-class gap persists
   under calibration residual `delta < 7/1920`.
5. **Lemma 6.1 and Proposition 6.2:** tangential calibration has an exact FDD
   modal-component interpretation, and a public VBL-VA001 triaxial split gives
   a separate exact rational `<24/25` versus `>=1` continuum certificate
   against every reducible symmetric coefficient family, plus a held-out modal
   replay.

The candidate priority claim is the conjunction of items 1--3: a fixed-channel,
fixed-order reducible/irreducible gap with a two-sided exponential joint
calibration/reducibility law.  Items 4--5 test exactness, robustness and physical meaning; they are not
substitutes for the asymptotic theorem.

## Three closest collision tests

The evaluator should compare theorem statements, not keywords.

1. G. Ljungars and M. Fu, *Design of Multi-Channel Linear Phase FIR Filters*
   (1998), Sec. 4.1 and Eqs. (4.5)--(4.6): matrix-valued real linear-phase FIR
   operator-norm minimax design on a frequency grid by SDP.  This kills any
   novelty claim for MIMO linear phase, operator-norm minimax design or SDP.
   Test whether its formulation nevertheless implies the fixed-order
   reducible-versus-irreducible separation.
2. J. Stefanovski and D. Georgijevic, *Interpolation with constraint on
   frequency region and systems & control application* (2016), Problem 1 and
   Theorem 1: real stable rational
   bitangential interpolation with arbitrarily small regional norm.  This
   kills novelty of small regional norm under tangential constraints.  Test
   whether its growing-order rational construction implies the fixed ordinary
   degree Hermitian-polynomial lower bound against the reducible subclass.
3. P. J. Kootsookos, *FIR(q) Filter Designs Using H-infinity Techniques*
   (1991), Problems 1.1--1.2 and Lemmas 3.4, 6.1: fixed-length MIMO FIR
   matrix-polynomial approximation in the uniform operator norm, with global
   lower bounds and linear-phase algorithms.  This kills novelty of that whole
   combination.  Test whether its target-approximation framework implies the
   directional interpolation obstruction for every fixed orthogonal scalar bank.

Ball--Kang (1990), Fuhrmann (2010), and Alpay--Lewkowicz (2014) are mandatory
adjacent checks: they establish tangential matrix-polynomial interpolation and
structured/Hermitian matrix-polynomial interpolation separately.

The bibliography and the clause-by-clause status of these and adjacent papers
are recorded in `programme/ROUTE_C_PRIORITY_MATRIX.md`.  Restricted full-text
access must be reported as uncertainty, never converted into positive evidence
of priority.

## Blind scoring request

Score **scientific contribution** and **manuscript quality** separately on a
0--10 scale.  Do not award points for repository size, the number of programs,
or exact replay of auxiliary facts.  Report:

- the exact theorem and proof step carrying the scientific score;
- the strongest primary-reference collision you found;
- whether the fixed-order comparator separation follows from that reference;
- whether the exponential construction is mathematically substantive or a
  padded interpolation example;
- whether full spark is a defensible generic calibration hypothesis or an
  application-destroying device;
- whether the FDD lemma gives meaningful task semantics without being
  overstated as a fault-diagnosis or deployment result;
- which claim, if any, would remain publishable after deleting all numerics;
- the strongest single reason the paper remains below 7, if either score is
  below 7.

## Required adversarial attacks

1. Check Lemma 2.2's invariant-line counting argument, including the `M=d+N`
   threshold, the exact role of full spark, and why a common invariant subspace
   of symmetric `3 x 3` coefficients supplies a common invariant line.
2. Recompute the norm of the symmetric interpolation right inverse used in
   Theorem 2.1, retaining the factorial in the Lagrange bound, and verify that
   the rational full-spark perturbation does not erase the exponential
   stopband estimate.  Run
   `verification/verify_constant_channel_scaling.py` as an arithmetic check,
   but inspect the proof rather than treating this script as a certificate.
3. Check the off-block decomposition leading to the `delta+beta_u` bound and
   verify that `beta_u=0` at more than `N` nodes forces coefficient-level
   reducibility.  Run the Lean root-count replay but do not mistake its scalar
   lemma for a formalization of the full theorem.
4. Try to construct a reducible counterexample, including a genuinely
   noncommuting `1+2` block, at the claimed number of pass signatures.
5. Run `verification/verify_quadratic_filter_witness.py` and
   `verification/verify_vbl_va001_witness.py` using only the standard library;
   distinguish their exact rational conclusions from CVXPY baselines.
6. Check that the cosine-to-FIR transformation really gives reciprocal
   palindromic taps, preserves coefficient commutation, and has the stated
   latency.
7. Inspect the deliberately failed graph-denoising pilot and decide whether it
   adequately limits the application claim.

## Explicit nonclaims

- Matrix-valued tangential interpolation is not new.
- MIMO FIR, linear phase, operator-norm minimax design and SDP are not new.
- Noncommuting matrix taps by themselves are not new.
- The VBL experiment is not a diagnostic-accuracy, damping-estimation,
  real-time, hardware or deployment study.
- The heuristic fixed-basis search is not a proof about all commuting filters;
  Lemma 2.2 supplies that proof.
- The manuscript does not claim to solve the Yang--Mills mass gap, the Riemann
  hypothesis or a Hamiltonian spectral-gap problem.

## Response block to fill

- PDF SHA-256 actually inspected:
- Title and page count:
- Commit or download URL:
- Scientific contribution score / 10:
- Manuscript quality score / 10:
- Headline theorem checked:
- Closest primary collision:
- Proof defect or counterexample found, if any:
- Strongest reason below 7, if applicable:
- Recommendation: reject / major revision / minor revision / submit:
