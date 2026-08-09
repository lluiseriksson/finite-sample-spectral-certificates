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

## Collision outcome so far

Tangential polynomial interpolation, matrix-valued Nevanlinna--Pick theory,
matrix Chebyshev optimization, and the block-Krylov/matrix-polynomial
correspondence are all established.  Searches and theorem abstracts examined
so far have not located the precise fixed-ordinary-degree, real-interval,
full-spark, commuting-versus-noncommuting exponential separation.  Absence from
this search is not evidence of priority.

## Gate G1 blockers

1. typeset the full proof with fixed norms and explicit `N^O(1)` stability
   bounds;
2. inspect the full statements of the closest norm-constrained tangential
   interpolation papers, not only abstracts;
3. benchmark the proved block-Krylov pass/stop proposition against independent
   scalar Chebyshev filtering on realistic eigenproblems.

The former fourth blocker, a rational or interval-certified non-affine finite
witness, is closed.  Until the remaining three close, Route C is the lead but
Gate G1 remains open.
