# Week 3 lab note: noncommuting tangential filters

Date: 2026-08-09

## Result found

For tangential pass constraints `P(lambda_j)v_j=v_j` and stopband `[-1,0]`, an
exact real-symmetric affine `2 by 2` filter has leakage

```text
(1+sqrt(61))/10 < 0.882,
```

whereas every feasible commuting real-symmetric affine filter has leakage at
least one.  A separate asymptotic construction now supports the candidate gap

```text
noncommuting symmetric leakage <= C exp(-cN),
commuting symmetric leakage       = 1,
channel dimension                 = N^2.
```

The construction uses explicit Lagrange--Chebyshev diagonal filters, a
symmetric tangential-interpolation right inverse on at most four nodes per
matrix entry, and an explicit moment-curve perturbation with
`epsilon=exp(-N^3)` to enforce full spark.

## Numerical stress test

The script `research/route_c_scaling_pilot.py` constrains every coefficient to
be real symmetric.  With degree 3, dimension 9, 12 targets, seed 7 and an
81-point optimization grid, it returned:

| Perturbation | SDP objective | 4001-point validation max | Max commutator norm | Full-spark SVD margin |
|---:|---:|---:|---:|---:|
| 0.001 | 0.2983490 | 0.2983682 | 0.466925 | 1.50e-5 |
| 0.010 | 0.3199793 | 0.3200281 | 0.391179 | 5.95e-5 |
| 0.030 | 0.3977121 | 0.3977588 | 0.708604 | 8.22e-5 |

All solves were reported `optimal_inaccurate` by Clarabel, so these numbers are
exploratory.  Interpolation residuals were below `5.1e-12`, and dense-grid
overshoot stayed below `4.9e-5`.  The commuting value one is an analytic lower
bound conditional on exact full spark; the reported SVD margins are not exact
certificates.

## Exact non-affine certificate

A degree-2, dimension-2 witness with four rational pass constraints is now
certified on the full stopband.  Its targets have slopes `-2,-1,1,2`, hence are
full spark.  Positive rational Bernstein coefficients for the principal minors
of `I+P` and `I-P` prove

```text
sup_{[-1,0]} ||P(x)||_op <= 301/304 < 1,
commuting symmetric optimum = 1.
```

The replay script uses only Python's exact `Fraction` arithmetic and is included
in CI as `verification/verify_quadratic_filter_witness.py`.

## Controlled Hermitian block benchmark

`research/route_c_block_krylov_benchmark.py` embeds the exact quadratic witness
in a deterministic `256 by 256` Hermitian problem with four pass eigenvalues
and 252 stopband eigenvalues.  After 100 degree-2 applications, the
noncommuting filter keeps pass geometry to `6e-16`, reduces the stop/pass ratio
from `5` to `0.0443`, and reduces the desired-subspace largest-angle sine from
`0.9892` to `0.0525`.  The exact-feasible commuting comparator is the identity.

The adversarial scalar degree-2 Chebyshev filter has stopband maximum `1/7` and
suppresses stop components much faster, but its pass multipliers are
`(1,17/7,31/7,7)`.  Even after optimal global rescaling, its pass-geometry error
is `0.8018` at 100 iterations and the desired-subspace sine is `0.9982`.  This
supports the theorem's controlled error interpretation but also exposes the
limitation: a natural application must genuinely require multipoint pass-row
geometry.  The construction is not evidence of a runtime advantage.

## Collision outcome so far

Tangential polynomial interpolation, matrix-valued Nevanlinna--Pick theory,
matrix Chebyshev optimization, and the block-Krylov/matrix-polynomial
correspondence are all established.  Searches and theorem abstracts examined
so far have not located the precise fixed-ordinary-degree, real-interval,
full-spark, commuting-versus-noncommuting exponential separation.  Absence from
this search is not evidence of priority.

## Gate G1 blockers

1. inspect the full statements of the closest norm-constrained tangential
   interpolation papers, not only abstracts;
2. find a naturally arising multichannel problem that mandates pass-row
   geometry, and test adaptive/nonstationary scalar baselines; the controlled
   Hermitian benchmark alone is not sufficient scientific payoff.

The quantified stability proof and the former rational non-affine-witness
blocker are closed.  Until the remaining two close, Route C is the lead but
Gate G1 remains open.

## FIR reformulation and failed natural graph pilot

The change of variables `x=(9 cos(omega)+5)/4` converts the asymptotic
polynomial into a reciprocal, real-symmetric, palindromic `2N+1`-tap MIMO FIR
response with common delay `N`.  The pass interval maps to low/mid frequencies
and `[-1,0]` maps to the high-frequency band
`[arccos(-5/9),pi]`.  Because both the affine power-basis change and the cosine
basis change are invertible, pairwise commutation is preserved in both
directions.  The commuting obstruction therefore becomes a lower bound for
every fixed orthogonal bank of scalar linear-phase filters at the same latency.

For the exact quadratic witness, `x=(3 cos(omega)+1)/2` produces the rational
palindromic taps now checked in
`verification/verify_quadratic_filter_witness.py`; the high-frequency bound is
still `301/304` versus one for the commuting class.

The 414-node grid-graph pilot failed honestly.  Symmetric taps collapsed to
identity, while a scalar Chebyshev baseline reached stop norm `1/17` with only
`0.00524` rescaled pass error at two iterations.  Consequently ordinary graph
denoising is rejected as the motivating application.  FIR frequency-direction
calibration is structurally better aligned, but a real calibration data set
and the FIR-specific collision audit are now the remaining application gates.
