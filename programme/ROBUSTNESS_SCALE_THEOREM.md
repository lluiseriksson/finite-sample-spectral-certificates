# Robustness-scale upgrade

Date: 2026-08-10

## Trigger

The 5.47/10 evaluation identified the superexponential perturbation
`epsilon_N = exp(-N^3)` as the main scientific weakness.  It asked whether the
class separation survives a polynomial full-spark margin or a
degree-independent calibration floor.

The response is not to hide the small margin.  The revised theorem proves two
new facts:

1. the construction can be made fully rational before normalization with
   `epsilon_N = 2^(-13N)`, quantitative full-spark margin
   `2^(-39N)/(48(N+4)^6)`, and an explicit exponentially small robust
   calibration radius; and
2. a degree-independent calibration floor is impossible for *any* signatures,
   because an explicit scalar commuting polynomial has both pass error and
   stopband leakage at most `exp(-2N/81)`.

Thus `-log(delta) = Theta(N)` is both sufficient for the constructed
noncommuting advantage and necessary in scale for any asymptotic class
separation.  The exponent constants are not claimed to be sharp.

## Proof chain

### Retaining the factorial

For a scalar right-inverse interpolation problem on `t` nodes, put `n=t-1`.
The exact cardinal-function estimate is

```text
(27S/5)^n 2^n / n!.
```

The old draft discarded `n!` and obtained `exp(O(N log N))`.  Since every
entry problem has `S-1 <= t <= 2S`, Stirling's elementary lower bound gives

```text
(27S/5)^n 2^n / n! <= (108e/5)^n,
```

and hence an `exp(O(N))` right inverse.  This is the decisive improvement.

### Quantitative full spark

For any three unnormalized signatures, expand their determinant as a cubic in
`epsilon`.  Its cubic coefficient is a nonzero Vandermonde determinant.  All
coefficients have denominator dividing `(M+1)^6`, every nonzero coefficient is
therefore at least `(M+1)^(-6)` in magnitude, and every coefficient is at most
18.  At `epsilon=2^(-13N)` the first nonzero term dominates the tail.  After
normalization and conversion from determinant to smallest singular value,

```text
gamma_N >= 2^(-39N) / (48 (N+4)^6).
```

This replaces the former transcendence-only argument with a quantitative
rational certificate.

### Robust near-reducibility lower bound

For any unit direction `u`, quantitative full spark implies that at most two
signatures have overlap below `gamma_N/sqrt(3)`.  Write `beta_u` for the largest
sampled coupling from `u` into its orthogonal complement.  The remaining
`N+1` nodes control the scalar degree-`N` compression `u^T Q(x) u`.  The
Lagrange sum at zero is at most `(24e/5)^N`, giving

```text
rho(Q) >= 1 - sqrt(3) (24e/5)^N (delta + beta_u) / gamma_N.
```

The lower bound is at least `3/4` for the explicit combined-error radius
recorded in the manuscript.  If `beta_u=0` at all `M>N` pass nodes, polynomial
root counting forces every coefficient to preserve the line spanned by `u`.
In three symmetric channels, any nontrivial common invariant subspace yields
such a line.  The exact comparator therefore strictly contains the commuting
class, including arbitrary noncommuting `2 x 2` blocks.

### Universal scalar barrier

Map the union of stop and pass intervals into `[0,2/9] union [4/9,1]` by
`r(x)=2(x+1)/9`.  The binomial-tail polynomial with threshold `ceil(N/3)` is
at most `exp(-2N/81)` on the first interval and differs from one by at most the
same amount on the second, by Hoeffding.  Multiplying it by the identity gives
a commuting comparator for arbitrary signatures.  Therefore constant-error
robustness is not a scientifically coherent target in this fixed-gap setting.

## Independent checks and failed routes

- Exact determinant polynomials were enumerated with rational arithmetic for
  every target triple through degree 24.  The smallest observed determinant to
  proved-envelope ratio was 5324.
- Neumann, correction, power-of-two, Lagrange, and Hoeffding envelopes were
  replayed through degree 512 by
  `verification/verify_constant_channel_scaling.py`.
- A pilot with random, nonclustered signatures found stopband optima below one
  already at degrees 3--6, suggesting that noncommutative advantage is not
  confined to the coordinate perturbation.  Higher-degree SDP solves became
  numerically ill-conditioned and were stopped; they are not evidence for the
  theorem.
- An attempted independent Fable High audit did not run: the authenticated
  organization returned HTTP 403 because subscription access to Claude Code
  is disabled.  No review credit is assigned to that attempt.

## Stronger finite certificate

A discrete search over small rational nodes and primitive integer directions
also found a much stronger auditable quadratic witness.  With pass nodes
`(1/2,1,3/2,7/2)` and directions
`(1,2),(0,1),(-1,2),(-1,1)`, its simple rational coefficients satisfy

```text
sup[-1,0] ||P(x)|| = 25/32,
commuting leakage >= 1 - 60 delta.
```

Strict separation therefore survives every `delta < 7/1920`, about 0.00365.
This is over 39 times the old certified radius `3/32300`, while the leakage
improves from about 0.9901 to 0.78125.  Nonnegative rational Bernstein
coefficients certify `(25/32)I +/- P(x)` on the whole interval; no grid or SDP
output enters the result.  The same object has an exact palindromic five-tap
realization and is replayed by
`verification/verify_quadratic_filter_witness.py`.

## Remaining limitations

- The sufficient and impossible exponential rates have very different
  constants; closing that exponent gap is open.
- The full-spark margin is exponential rather than polynomial.  The scalar
  barrier proves that constant calibration tolerance is impossible, but does
  not rule out a construction with polynomial geometric margin and exponential
  calibration radius.
- The proof is exact mathematics plus standard-library rational replay, not a
  kernel-checked formalization.
- Priority remains qualified where the two closest publisher full texts are
  inaccessible.
