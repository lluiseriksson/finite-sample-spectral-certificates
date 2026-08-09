# Theorem draft: adaptive covariance-exact localizer certificates

These are research statements under audit, not yet publication claims.

## Setup

Let `y in R^p` collect all independent measured entries of a finite scalar or
block-matrix moment window.  Fix a proposed support edge `theta` and degree
`N`.  Its localizer is an affine symmetric matrix

    L(y) = L_0 + sum_k y_k L_k  in Sym(r).

Let the joint confidence region be the ellipsoid

    E = { y_hat + R z : ||z||_2 <= q }.

For `X >= 0`, define `ell(X)_k = tr(L_k X)` and

    U_E(X) = tr(L(y_hat) X) + q ||R^T ell(X)||_2.

The normalization `tr X = 1` is a compact slice through every nonzero PSD ray.

## Theorem A: exact robust truncated infeasibility

Assume `E` is nonempty and compact.  Then

    E intersect {y : L(y) >= 0} is empty

if and only if

    Phi(E) := min_{X >= 0, tr X = 1} U_E(X) < 0.

Thus the SDP/SOCP is necessary and sufficient for ruling out every moment
vector in the confidence ellipsoid that satisfies this degree-N localizer
condition.  It is not merely a sufficient polynomial witness.

### Proof route

The image `K=L(E)` is compact and convex.  If `K` is disjoint from the closed
PSD cone, strong separation supplies a symmetric separator.  Finiteness on the
unbounded PSD cone forces the separator to lie in the dual cone, which is the
PSD cone itself; rescale it to trace one.  The support function of `E` is
exactly `U_E(X)`.  Conversely, `U_E(X)<0` implies `tr(L(y)X)<0` for every
`y in E`, impossible when both `L(y)` and `X` are PSD.

Audit point: write the separation orientation carefully and handle the trivial
`r=0`/zero-separator cases explicitly.

## Corollary B: noiseless reduction and rank

For a singleton `E={y}`, the optimum is

    Phi({y}) = lambda_min(L(y)),

attained by a rank-one projector onto a minimum eigenvector.  With uncertainty,
an optimal `X` can have higher rank: it is an ensemble/separating certificate,
not necessarily one polynomial direction.  Any claim of rank-one exactness in
the noisy problem would be false without additional hypotheses.

## Theorem C: adaptive coverage without data splitting

Suppose the reported random set `E(Data)` covers the true moment vector with
probability at least `1-alpha`.  Under the null `L(y_true)>=0`, the adaptive
test

    reject iff Phi(E(Data)) < 0

has type-I error at most `alpha`, even though `X` is optimized after observing
the same data.  Indeed, rejection is impossible on the coverage event.

This is the key distinction from choosing a filter on the data and then using
an unadjusted pointwise error bar.

## Corollary D: exact unknown-covariance Gaussian test

Let `Y_1,...,Y_n` be iid `N_p(y,Sigma)`, with `n>p`, sample mean `y_bar`, and
unbiased sample covariance `S`.  Put

    q_H^2 = p(n-1)/(n-p) F_{p,n-p}(1-alpha)

and

    E_H = { y : n(y_bar-y)^T S^{-1}(y_bar-y) <= q_H^2 }.

Then `E_H` has exact coverage `1-alpha`, so Theorem C gives a finite-sample
level-alpha adaptive localizer test with estimated covariance and no data
splitting.

Audit points: specify nonsingularity almost surely for positive-definite Sigma;
state what changes under singular covariance, batching, Markov-chain samples,
or non-Gaussian observations.

## Extension E: full truncated Hausdorff feasibility

The physical null is an intersection of several affine conic constraints, for
example `H_N(y)>=0`, `G_N(y)>=0`, and `theta H_N(y)-G_N(y)>=0`.  The exact
ellipsoidal feasibility problem is itself an SDP/SOCP.  Its conic dual uses one
PSD multiplier per constraint.  A production paper should compare:

1. the single-localizer separator above;
2. the full Hausdorff-feasibility primal;
3. its dual certificate and strong-duality conditions.

The full formulation may certify incompatibility when no individual localizer
has a common negative separator.

## Theorem F: non-Gaussian Gaussian-sketch covariance

Let `u_{k,a}=T^(k/2) psi_a`, draw `g~N(0,I)`, and observe

    Z_{k,ab}=<g,u_{k,a}><g,u_{k,b}>.

Then `E Z_{k,ab}=B_k[a,b]`.  If

    K[(k,a),(ell,c)] = <u_{k,a},u_{ell,c}>,

Isserlis' identity gives the exact feature covariance

    Cov(Z_{k,ab},Z_{ell,cd})
      = K[(k,a),(ell,c)] K[(k,b),(ell,d)]
      + K[(k,a),(ell,d)] K[(k,b),(ell,c)].

For the mean of `n` independent sketches this covariance is `C/n`.  Without
assuming Gaussian feature vectors, multivariate Chebyshev gives

    P((m_bar-m)^T (C/n)^dagger (m_bar-m) >= r/alpha) <= alpha,

provided the error is supported on the rank-`r` range of `C`.  Combining this
ellipsoid with Theorems A and C yields a distribution-free level-alpha block-
localizer test for the sketch model.  Singular covariance requires explicit
range constraints; silently adding a ridge changes the theorem.

Audit points: prove the range statement, distinguish mathematical rank from a
floating-point cutoff, and provide an interval/rational certificate for the
production covariance factorization.

## Theorem G: covariance-free finite-sample Wishart certificate

Stack all sketch coordinates into `z in R^d`.  In the Gaussian-sketch model,
`z_1,...,z_n` are iid `N(0,K)` for an unknown positive-semidefinite Gram
matrix `K`, and the empirical Gram matrix is

    S = (1/n) sum_r z_r z_r^T.

Let

    eta = (sqrt(d) + sqrt(2 log(2/alpha))) / sqrt(n) < 1,
    a = (1+eta)^(-2),  b = (1-eta)^(-2).

The standard extreme-singular-value tail bound for a rectangular Gaussian
matrix implies, simultaneously with probability at least `1-alpha`,

    a S <= K <= b S.

This remains valid when `K` is singular by restricting the Gaussian matrix to
the range of a factor of `K`; using the ambient `d` in `eta` is conservative.

More explicitly, write the `n x d` sketch matrix as `Z=G K^(1/2)` with `G`
standard Gaussian.  On the event

    sqrt(n)-sqrt(d)-t <= s_min(G) <= s_max(G)
                         <= sqrt(n)+sqrt(d)+t,

whose failure probability is at most `2 exp(-t^2/2)`, congruence by `K^(1/2)`
gives `(1-eta)^2 K <= S <= (1+eta)^2 K`.  Taking
`t=sqrt(2 log(2/alpha))` and rearranging yields the displayed band.  The final
paper must cite the exact theorem used and check whether its lower-tail formula
uses `sqrt(d)` or `sqrt(d-1)`; the more conservative valid version will be used.

Let `P(K)` extract the diagonal time blocks `K[(k,:), (k,:)]`, which are the
block moments `B_k`, and let `L_theta(P(K))` be the proposed upper-support
localizer.  The test

    reject iff there is no K such that
        a S <= K <= b S  and  L_theta(P(K)) >= 0

is therefore a finite-sample level-`alpha` test with *unknown covariance*.
It is a single semidefinite feasibility problem.  Unlike Hotelling, it does
not assume that the quadratic feature vector is Gaussian; unlike Theorem F's
Chebyshev ellipsoid, it does not require oracle knowledge of its covariance.

Audit points: cite/prove the Gaussian singular-value inequality with constants;
derive and numerically check the conic dual; distinguish centered and known-zero
mean sketches; study whether the Loewner band is too conservative in practice.

### Semidefinite alternative

After a fixed diagonal congruence scaling, write the confidence interval as
`A <= Q <= B` and the localizer as `L_D(Q)>=0`.  A dual certificate consists of
`Y_A,Y_B,Y_L>=0` satisfying

    Y_A - Y_B + L_D^*(Y_L) = 0

and

    -<Y_A,A> + <Y_B,B> < 0.

Indeed, pairing the three primal inequalities with these multipliers makes all
terms containing `Q` cancel, while the remaining constant would have to be
nonnegative under primal feasibility.  The implementation fixes the homogeneous
scale by `tr(Y_A)+tr(Y_B)+tr(Y_L)=1` and minimizes the constant.  When the
empirical Gram matrix is positive definite, the interval has a strict midpoint
and is compact; standard conic duality excludes weak infeasibility.  Singular
sample Grams require an explicit reduction to their range before making that
claim.

For approximate multipliers with stationarity residual `R`, feasibility would
imply `0 <= beta + tr(R Q)`.  Since the scaled interval obeys `0<=Q<=B`,

    |tr(R Q)| <= ||R||_* ||B||_op.

Thus `beta + ||R||_* ||B||_op < 0` is an a-posteriori incompatibility
certificate even before exact rational reconstruction.  The production code
first projects multipliers PSD, repairs `R=R_+-R_-` by adding `R_-` to the
lower-band multiplier and `R_+` to the upper-band multiplier, adds an explicit
roundoff-scale identity shift, and finally applies this nuclear/operator-norm
bound.

## Theorem H: uncertainty can change the optimal witness

For a population localizer `M` and mean covariance `C/n`, define

    kappa_PSD = sup_{X>=0,tr X=1} [-tr(MX)]_+ /
                sqrt(ell(X)^T C ell(X)).

The known-covariance Chebyshev certificate becomes negative exactly when
`sqrt(n) kappa_PSD > sqrt(p/alpha)` (with the zero-denominator cases treated
separately).  Thus `kappa_PSD` is the certificate signal-to-noise ratio and
the corresponding population sample threshold is explicit.

The noiseless minimum-eigenvector need not maximize this ratio.  Moreover the
PSD relaxation can genuinely prefer a mixed witness: in the synthetic family
`M=-epsilon I_d` with isotropic Frobenius uncertainty, the optimum is `I_d/d`
and its uncertainty norm is `1/sqrt(d)`, whereas every rank-one projector has
norm one.  Hence mixed certificates can reduce required sample size by a
factor `d`.  In the current ANNNI pilot, however, many optimized witnesses near
the detection boundary are rank one; the empirical advantage there must be
described as covariance-aware versus the plug-in eigenvector, not falsely as
mixed versus pure.

## Model-specific work still missing

These results solve the statistical/adaptive layer but are not yet a 7+ physics
paper.  The remaining required producer is an interacting-model campaign that
shows what new gap statement becomes possible, how power scales with volume
and operator support, and where covariance estimation or visibility breaks.
