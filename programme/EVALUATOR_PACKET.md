# Adversarial evaluator packet: Route C research draft

This cover sheet identifies the exact artifact to inspect and asks for a blind,
claim-level assessment.  It is not a certificate of novelty and it is not a
request to reward repository size.

## Artifact identity

- Repository: <https://github.com/lluiseriksson/finite-sample-spectral-certificates>
- Research branch: `research/endogenous-visibility`
- Draft pull request: <https://github.com/lluiseriksson/finite-sample-spectral-certificates/pull/2>
- Source commit inspected when this packet was prepared:
  `7d31b2f0ceaf336dbe6e7f2b2511dddd01062530`
- PDF path: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Title visible on page 1: *Exponential Tangential Advantages of
  Noncommuting Matrix-Polynomial and Linear-Phase MIMO FIR Filters*
- Page count: 13
- SHA-256 of the repository PDF bytes:
  `dda469dc4c1714aab8bc519f2356f2f6dd42d8423cf2ed89775b5c8fbfafd2dc`
- Embedded manuscript status: `Research draft, 9 August 2026`
- CI status at the commit above: push run `31335423165` and pull-request run
  `31335426123` both passed.

Before scoring, independently hash the downloaded PDF and copy the result into
the response.  If the hash, title or page count differs, stop: the score belongs
to another artifact.  If a later commit changes the PDF, this packet must be
regenerated rather than silently reused.

## Claims actually submitted for assessment

1. **Theorem 2.1:** for every sufficiently large degree `N`, already at fixed
   dimension `d=3` and `M=N+3`, there are full-spark tangential pass
   constraints for which a real-symmetric noncommuting matrix polynomial
   attains stopband leakage `C exp(-cN)`, whereas every
   feasible pairwise-commuting real-symmetric polynomial is identically the
   identity and has leakage one.
2. **Corollary 4.1:** the same construction gives a fixed-latency reciprocal
   linear-phase MIMO FIR separation between fully coupled taps and every fixed
   orthogonal bank of scalar linear-phase FIR filters.
3. **Propositions 5.1--5.2:** an exactly replayable rational five-tap witness
   has continuum leakage at most `301/304`, and the commuting obstruction
   persists under calibration residual `delta < 3/32300`.
4. **Lemma 6.1 and Proposition 6.2:** tangential calibration has an exact FDD
   modal-component interpretation, and a public VBL-VA001 triaxial split gives
   a separate exact rational `<24/25` versus `1` continuum certificate plus a
   held-out modal replay.

Only item 1, together with its fixed-latency interpretation in item 2, is the
candidate priority claim.  Items 3--4 test exactness, robustness and physical
meaning; they are not substitutes for the asymptotic theorem.

## Three closest collision tests

The evaluator should compare theorem statements, not keywords.

1. G. Ljungars and M. Fu, *Design of Multi-Channel Linear Phase FIR Filters*
   (1998), Sec. 4.1 and Eqs. (4.5)--(4.6): matrix-valued real linear-phase FIR
   operator-norm minimax design on a frequency grid by SDP.  This kills any
   novelty claim for MIMO linear phase, operator-norm minimax design or SDP.
   Test whether its formulation nevertheless implies the fixed-order
   commuting-versus-fully-coupled separation.
2. J. Stefanovski and D. Georgijevic, *Real stable all-pass and minimum-phase
   solutions to the bitangential interpolation problem with a low norm on a
   frequency region* (2016), Problem 1 and Theorem 1: real stable rational
   bitangential interpolation with arbitrarily small regional norm.  This
   kills novelty of small regional norm under tangential constraints.  Test
   whether its growing-order rational construction implies the fixed ordinary
   degree Hermitian-polynomial lower bound against the commuting subclass.
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

1. Check Lemma 2.2's counting argument in a common eigenbasis, including the
   `M=d+N` threshold and the exact role of full spark.
2. Recompute the norm of the symmetric interpolation right inverse used in
   Theorem 2.1 and verify that the full-spark perturbation does not erase the
   exponential stopband estimate.  Run
   `verification/verify_constant_channel_scaling.py` as an arithmetic check,
   but inspect the proof rather than treating this script as a certificate.
3. Look for circularity: the noncommuting coefficients used to prove the upper
   bound must satisfy the same constraints as the commuting lower-bound class.
4. Try to construct a nonidentity commuting counterexample at the claimed
   number of pass signatures.
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
