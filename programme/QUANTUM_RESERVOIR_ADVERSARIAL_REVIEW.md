# Adversarial release review

## Decision

Internal provisional score: **6.15/10**, with a defensible interval of
**5.95--6.30**.  This clears the requested 6+ threshold on the frozen artifact
and is materially above the 5.72 assessment of the preceding MIMO paper.  It
is not presented as an independent referee score.

## Why the score is higher than the preceding paper

1. The central root-count separation now lives inside an explicit causal inner
   six-port completion instead of an unconstrained Hermitian polynomial class.
2. The bath congruence converts the signal gap into a squared Kossakowski-rate
   gap under visible spectral assumptions.
3. The dwell-time theorem and delay-budget corollary identify the resource
   paid for exponential suppression; the paper does not advertise a free
   coherence gain.
4. The exact quadratic no-go explains why the preceding zero-phase witness
   cannot be rescaled into a globally passive device.
5. A passive-preserving tolerance experiment supplies a negative result:
   stop suppression is stable in the frozen component model, but nominal pass
   calibration is fragile.
6. Priority collisions are named directly: quantum filter functions,
   passivity-preserving tangential interpolation, passive reservoir
   engineering, and Blaschke--Potapov network modes are not claimed as new.

## Main remaining attacks

- The comparator theorem is exact, not uniformly robust to calibration error.
- The Kossakowski interface is linear, stationary and weak-coupling; it is not
  a theorem for arbitrary interacting reservoirs.
- The six-port construction is a discrete-time/traveling-field network, not a
  fabricated device or a continuum QFT.
- The specific conjunction appears distinct in the literature searched, but
  exhaustive priority cannot be certified without specialist review.
- The largest numerical orders are ill-conditioned and are included as an
  asymptotic audit, not a fabrication claim.

## Release gates

- [x] analytic node count for every integer `S >= 5` at `alpha = 9/20`;
- [x] full-spark Vandermonde proof;
- [x] explicit inner completion and global Schur proof;
- [x] rational root-count comparator theorem;
- [x] squared bath-rate corollary with assumptions visible;
- [x] general delay lower bound and delay-budget frontier;
- [x] exact no-go and deterministic numerical replay;
- [x] passive-preserving tolerance stress test with frozen seed;
- [x] primary-source priority audit and explicit non-claims;
- [x] clean 13-page A4 PDF, rendered and inspected page by page;
- [x] immutable artifact hashes and one-command release gate.

Independent frontier-model review was attempted, but the available bridge
returned an organization-level 403 before any manuscript content was sent.
No billable fallback or account rotation was used.
