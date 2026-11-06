# Optional sign test for raw elliptical sketches

`monomial_sign_test.py` is a separate, explicitly enabled, standard-library
entry point. It tests **one fixed monomial and channel direction** in the
12-dimensional raw sketch model (degree two, two channels). It does not
replace the Wishart band, optimize a separating direction on test data, or
reinterpret the existing Student misspecification campaigns.

## Model and null

Write each raw row as six two-channel blocks `z_0,...,z_5`. Clean rows are
iid, elliptical about a **known zero centre**, with finite covariance `K`.
Radial Student laws with degrees of freedom greater than two are included;
independent Student coordinates and arbitrary symmetric laws are not covered.
The eighteen quadratic features are not raw elliptical rows.

Let `B_k = K[k,k]` denote the diagonal time blocks and
`L[i,j] = theta B_(i+j) - B_(i+j+1)`, for `i,j=0,1,2`. The null is `L >= 0`.
Choose one `i` and nonzero two-channel `v`, independently of the clean test
sample, and put `X = v' z_(2i)`, `Y = v' z_(2i+1)`. For every covariance,

```
E[theta X^2 - Y^2] = v' L[i,i] v >= 0.
```

For `theta > 0`, the expression is the product of two linear forms
`(sqrt(theta) X-Y)(sqrt(theta) X+Y)`. The elliptical angular identity gives
`P(product < 0) <= 1/2` whenever its expectation is nonnegative. For
nondegenerate forms their correlation is nonnegative and this probability
is `P(R>0) arccos(rho)/pi`; zero projections and radial mass at zero only
reduce the negative probability. Hence the upper tail of `Bin(n,1/2)` at
the observed negative count is a conservative p-value. The code compares
`theta X^2-Y^2` directly, without approximating a square root.

This is an application of the classical quadrant/arcsine identity and
binomial sign test, not a new general method. For the elliptical identity,
see Langworthy et al., *Canonical correlation analysis for elliptical
copulas*, J. Multivar. Anal. 183:104715, sections 2 and 3.1,
[doi:10.1016/j.jmva.2020.104715](https://doi.org/10.1016/j.jmva.2020.104715).

## Selection, zeros, and measurement errors

- Fix `n`, `theta`, the degree, direction and test level before viewing the
  test sample. An independent training sample may select the direction if
  the null implies the constraint for every permitted selection. This does
  not permit scanning support thresholds and claiming simultaneous coverage.
- A centre estimated from training is not a known zero centre. There is no
  centring-error correction here.
- Keep ties and uncertain signs in `n`. Do not drop rows, stop early, choose
  the smallest p-value, or run multiple unadjusted tests.
- Supply rational enclosing intervals for the true raw coordinates. Only
  a product whose entire interval is negative counts. The count can be
  smaller than the true count; upper-tail testing remains conservative.
  Declaring that a binary float is exact does not establish measurement
  accuracy. The program does not construct certified measurement bounds.
- If a genuinely known bound `b` limits arbitrary row replacements in the
  clean sample, use `max(0, count-b)`. This follows from deterministic count
  dominance, even for adversarial replacements. It does not justify a
  direction chosen through a corrupted training process that reveals the
  clean holdout.

The required declarations are an explicit contract, **not model validation**.
Nonrejection is inconclusive; it certifies neither restricted support nor a
positive spectral gap. The test may miss violations visible to multigrade
localizers or the SDP. No power comparison or new simulation is supplied.

## Why only one monomial

For a general polynomial direction `c=(c_0,c_1,c_2)`, replacing `c' L c` by
`theta E(sum c_i' z_(2i))^2 - E(sum c_i' z_(2i+1))^2` requires additional
half-moment identities between off-diagonal time blocks. The full-covariance
Wishart null does not impose those identities.

For example, let `theta=9/16`, `B_k=theta^k I_2`, and let all off-diagonal
entries of `K` vanish except the symmetric pair between coordinates `(1,0)`
and `(3,0)`, set to `theta^2/2`. Then `K` is positive definite and `L=0`, but
`c_0=c_1=(1,0)`, `c_2=0` makes that proposed two-form mean `-theta^2`.
Thus it would reject a valid full-covariance null. The regression test
retains this counterexample. Numerical eigenvalue clipping is not a proof
that exact half-moment identities remain true.

## Input and execution

The JSON object contains `degree`, two-entry `direction`, `theta`, `n`,
optional `alpha` (default `"1/20"`), optional `max_replaced` (default zero),
`assumptions` and `rows`. Every row has twelve `[lower, upper]` pairs in
time-major, then channel-major order. Use integers or quoted rational strings
such as `"-2/3"`. Unquoted decimals, booleans and nonfinite numbers are rejected.

`assumptions` must contain exactly these five declarations:

```json
[
  "clean_iid_elliptical_finite_covariance",
  "known_zero_centre",
  "pair_theta_fixed_without_clean_holdout",
  "fixed_n_single_test",
  "certified_intervals_for_raw_sketches"
]
```

```sh
python monomial_sign_test.py contract-and-rows.json --enable
python -m unittest discover -s verification -p test_monomial_sign.py -v
```

The entry point refuses to test without `--enable`. It prints the exact
rational p-value, decision, counts and input hash. The hash identifies the
input; it does not prove that the choices were made before data collection.
The tests use exact algebra and fixed fixtures, not a coverage simulation.
