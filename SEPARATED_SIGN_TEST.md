# Optional test for one separated-time pair

`separated_sign_test.py` accepts exactly one fixed pair `[0,3]`, `[2,5]`,
or `[0,5]`. It extends the available population inequalities through exact
SOS identities. The original `monomial_sign_test.py` degree API is unchanged.
This is an application of polynomial factorization and the classical sign
test, not a new general statistical method or a new power benchmark.

## Full-covariance justification

Use six two-channel raw blocks `z_0,...,z_5`, their diagonal covariance
blocks `B_k = E[z_k z_k']`, and the null block localizer
`L[i,j] = theta B_(i+j) - B_(i+j+1) >= 0`, `i,j=0,1,2`.
The full 12-dimensional covariance may have arbitrary cross-time blocks.
No half-Hankel equalities are assumed. For `theta > 0`, put

```
q2(t) = t^2 + theta*t + theta^2
      = (t + theta/2)^2 + (3/4)*theta^2

q4(t) = t^4 + theta*t^3 + theta^2*t^2 + theta^3*t + theta^4
      = (theta^2 + theta*t/2)^2
        + (3/4)*(theta*t + (2/3)*t^2)^2 + (2/3)*t^4
```

The identities `(theta-t)q2 = theta^3-t^3`,
`(theta-t)t^2 q2 = theta^3 t^2-t^5` and
`(theta-t)q4 = theta^5-t^5` show respectively

```
theta^3 B_0 - B_3 >= 0
theta^3 B_2 - B_5 >= 0
theta^5 B_0 - B_5 >= 0.
```

Indeed, apply the linear matrix map `T(t^k)=B_k`. For every polynomial
`p(t)=sum c_i t^i` of degree at most two and channel vector `v`,
`v' T((theta-t)p^2) v = sum c_i c_j v' L[i,j] v >= 0`.
This map is purely algebraic and does not assume a representing measure.
The positive weights in the decompositions prove the matrix inequalities.
The rational certificate and replay are in
`verification/separated_time_certificate.json` and
`verification/check_separated_sos.py`.

For the selected pair `(a,b)`, take `X=v'z_a`, `Y=v'z_b`. The tested
quantity is `theta^(b-a) X^2-Y^2`, the product of the two linear forms
`theta^((b-a)/2)X-Y` and `theta^((b-a)/2)X+Y`.
Under the [original elliptical contract](MONOMIAL_SIGN_TEST.md), its
nonnegative mean implies negative probability at most one half. The same
binomial upper tail therefore applies. Exact rational comparisons use
`theta^(b-a)` directly, without an approximate square root.

The classical quadrant/arcsine relation is attributed to Langworthy et al.,
*Canonical correlation analysis for elliptical copulas*, J. Multivar. Anal.
183:104715, sections 2 and 3.1,
[doi:10.1016/j.jmva.2020.104715](https://doi.org/10.1016/j.jmva.2020.104715).

## What this adds, and what it cannot certify

For `theta=1/2`, let `K` be block diagonal with
`B=(1,1/2,1,1/2,1,1/2) I_2`. The three original diagonal inequalities
all vanish, but `theta^3 B_0-B_3=-(3/8)I_2`. Thus the additional pairs
can detect a population violation missed by the old inequalities. This is
an algebraic example, not an observed or uniform finite-sample power gain.

The resulting six inequalities are still incomplete. With
`B=(1,1/2,1/10,1/20,1/100,1/200) I_2`, all six are nonnegative, yet
`c'L c=-3/10` for `c=((1,0),(-1,0),(0,0))`. Both examples have positive
definite block-diagonal `K`. Passing these tests certifies neither the
full null, restricted support nor a positive spectral gap.

## One test, with an explicit contract

Fix the pair, channel direction, `theta`, sample size and level before the
clean holdout. Independent training may select them only when the required
independence and null premise hold. Selecting a support threshold does not
produce simultaneous confidence over thresholds. The entrypoint accepts
one pair; it does not scan pairs, choose a minimum p-value or adjust multiple
calls. Running it repeatedly on the same holdout needs its own valid
multiplicity procedure, which this interface does not supply.

Clean rows must be iid elliptical about a known zero centre, with finite
covariance. Quadratic features, estimated centring and arbitrary symmetric
laws are not covered. The declarations are not verified by the software.
All twelve raw coordinate intervals are validated. Ties and ambiguous signs
stay in `n`, and only certified negatives count. The existing genuinely
known row-replacement bound `max_replaced` remains available. See the original
documentation for the independence, contamination and measurement limits.

## Input

Use the original JSON contract, replacing `degree` with `time_pair`, e.g.
`"time_pair": [0,3]`. Keep `direction`, `theta`, `n`, `assumptions`, `rows`
and optional `alpha`/`max_replaced`. The required five assumptions are exactly
those in `MONOMIAL_SIGN_TEST.md`. Values must be integers or quoted rational
strings; the raw intervals must enclose the true coordinates. Floats,
unsupported pairs, lists of pairs, `degree` and scan options are rejected.

```
python separated_sign_test.py contract-and-rows.json --enable
python -m unittest discover -s verification -p "test_*sign.py" -v
```

The result records `time_pair`, the original `theta`, `quadratic_weight`,
counts, exact p-value and decision. It contains no monomial `degree` field.
The input hash identifies bytes; it is not proof of prior registration.
The existing conservative counting kernel is reused after certified linear
interval projection. No sampling, optimization or model fitting is performed.
