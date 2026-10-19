# Research record: exact memory of direct-sum spectral fan-out

## Research decision

The selected tenth-paper direction closes the principal geometric gap left by
the orthogonal fan-out theorem. It does not move to an unrelated conjecture.
The decisive observation is that orthogonality is a sufficient coordinate
choice, while direct-sum position is the exact generic mechanism forcing the
maximum degree.

## Proved statements

For distinct boundary nodes, a fixed input `k`-plane, and arbitrary output
`k`-planes with joint span dimension `r`, the minimum degree of a square
finite rational-inner interpolant satisfies

```text
r - k <= d_min <= k(L - 1).
```

The upper bound comes from a positive block Pick completion with a zero Schur
complement. The lower bound has two audits: a common-range argument from a
minimal realization and a Pick-Stein inertia argument. Therefore direct-sum
targets (`r = Lk`) have exact degree `k(L - 1)`. When `N >= Lk`, this condition
is open and dense.

The paper also proves:

- direct-sum targets can converge to one common plane while retaining maximum
  exact degree for every nonzero opening;
- the distinguishing projector-sum eigenvalue closes quadratically;
- a gauge-invariant, fail-closed noisy eigenvalue certificate;
- the exact three-line phase diagram for `L=3`, `k=1`.

## Novelty boundary

Boundary matrix Schur interpolation, lurking-isometry realizations, positive
semidefinite completion, and rank/degree analysis are established subjects.
The paper does not claim those tools as new. The claimed contribution is the
exact direct-sum law, its universal sandwich, the collision discontinuity,
the operational projector certificate, and the singular three-line
classification. A directed literature audit found no primary source stating
this conjunction.

## Adversarial checks incorporated

- Replaced ambiguous “full-span” terminology by “direct-sum targets”.
- Added the explicit factorization of the zero-Schur completion.
- Added and proved deletion of decoupled unit-circle realization modes.
- Replaced calibrated block-Gram noise by the gauge-invariant sum of target
  projectors.
- Added clustered-frequency fixtures and separated exact rank from floating
  rank and conditioning.
- Checked transfer unitarity on a 257-point boundary grid.
- Stated that double-precision colligation synthesis is not certified at the
  smallest node gaps.

## Reproduction

The producer, frozen JSON certificate, independent verifier, paper source,
figure, and hash manifest are all part of the release gate. The verifier does
not import the producer.
