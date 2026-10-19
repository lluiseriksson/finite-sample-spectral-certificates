# Research memo: exact action--memory regions

Date: 2026-08-10

## Selection decision

The eighth-paper search swept the public `lluiseriksson` programme, including
the passive-network series, finite-window spectral certificates, GKLS
maintenance bounds, transfer-matrix formalizations, Yang--Mills repositories,
zero-free-region work, rooted-tree expansions, and Riemann/resolvent projects.

Four serious cross-repository candidates survived the first pass:

1. a passive--active maintenance trilemma for filled SSH chains;
2. a flatness/blindness boundary for finite-window OS spectroscopy;
3. a zero-free-region/passive-delay dictionary;
4. an exact joint action--memory region for rational inner networks.

The fourth was selected. It had the lowest dependence on unproved physical
interfaces, a complete converse, closed-form sharp constants, explicit
saturators, and a direct falsification protocol. The passive--active SSH law
remains a promising later application but presently inherits more model
assumptions and is algebraically closer to a composition of existing bounds.

## Main result and proof gates

For endpoint principal angles `beta_j`, define `B = sum beta_j` and let `r`
count the positive angles. For proper-arc trace action `A` and degree budget
`d`, the selected theorem is

```text
B > 0:  feasible iff d >= r and 2B <= A <= 2*pi*d - 2B.
B = 0:  feasible iff A = 0 or d >= 1 and 0 < A < 2*pi*d.
```

The proof was not accepted until five independent gates closed:

- the forward trace-speed bound `A >= 2B`;
- the complementary-arc return bound `2*pi*n - A >= 2B`;
- the rank obstruction `n >= r` from rank-one Potapov factors;
- a scalar degree-one factor realizing every phase in `(0, 2*pi)` on a fixed
  proper arc;
- a rank-one unitary gate moving a line by `alpha` exactly when its phase lies
  in `[2*alpha, 2*pi - 2*alpha]`.

Subdividing principal geodesics and taking the Minkowski sum of the gate
intervals proves sufficiency, including both closed faces.

## Precedence audit

The audit found classical precedence for every ingredient:

- Qiu--Zhang--Li for symmetric-gauge Grassmann metrics;
- Antezana--Larotonda--Varela for optimal unitary/Grassmann paths;
- Potapov and Alpay--Jorgensen--Lewkowicz for rank-one para-unitary
  factorization;
- Bolotnikov--Dym for matrix tangential boundary interpolation;
- Bharath et al. for matrix all-pass design with nodewise group-delay data.

No source found in the targeted audit states the joint attainable diamond,
its complementary-return upper face, or the exact minimum-degree formula at
prescribed arc action. Accordingly, the manuscript claims novelty only for
that conjunction. It explicitly does not claim a new Grassmann metric,
Potapov factorization, Pick theorem, or generic quantum speed limit.

A proposed multi-node minimum-rank theorem was demoted after the audit showed
that it is a gauge-quotient corollary of Bolotnikov--Dym. The paper retains the
classical corollary transparently because it yields a useful exact example:
three pairwise degree-one routing tasks require degree two jointly, detected
by a non-real cycle product.

## Verification contract

The producer and verifier are separate programs. The verifier does not import
producer functions. Frozen checks cover:

- exact lower and upper faces;
- orthogonal collapsed diamonds;
- zero principal angles and the discontinuity at zero motion;
- the open upper endpoint for coincident subspaces;
- analytic derivatives, adaptive phase integration and full winding;
- noncommuting Potapov factor order;
- a unitary, controllable and observable degree-two colligation;
- a fail-closed noisy cycle-holonomy witness.

The artifact manifest hashes the manuscript, PDF, code, certificate, figure,
workflow and this memo. The release gate regenerates the certificate in a
temporary directory, compares it numerically with the frozen artifact, runs
the independent verifier, and checks hashes before and after replay.

## Scope and falsifiers

The result assumes square finite rational inner transfer matrices regular on
the unit circle. Loss, active gain, singular inner factors, pure propagation
delays, and nonrational networks are outside scope. Trace action is a
Wigner--Smith phase/dwell resource, not thermodynamic work.

Within scope, any example with `A < 2B`, `A > 2*pi*d - 2B`, or `d < r` would
refute necessity. Any point inside the region not realizable by the explicit
compiler would refute sufficiency.
