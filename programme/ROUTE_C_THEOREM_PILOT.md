# Route C theorem pilot: tangential matrix filters

Status: a strict noncommutative separation has been found.  It is a pilot, not
yet a priority claim; tangential Nevanlinna--Pick and matrix interpolation are
serious collision risks.

## 1. Extremal problem

For a matrix polynomial `P(x)` and stopband `[-1,0]`, impose tangential pass
constraints

```text
P(lambda_j) v_j = v_j.
```

Minimize `sup_{x in [-1,0]} ||P(x)||_op`.  The constraints model several
spectral atoms with different channel vectors.  They differ essentially from
the full normalization `P(lambda)=I`: under full normalization, compression to
any vector proves that a scalar Chebyshev polynomial times the identity is
already optimal.

## 2. Exact 2 by 2 affine separation

Take degree one, `P(x)=A+xB`, with

```text
lambda_1 = 1/2,   v_1 = (1,0)^T,
lambda_2 = 1,     v_2 = (1,1)^T.
```

The rational matrices

```text
A = [[4/5,  1/5],       B = [[ 2/5, -2/5],
     [1/5, 3/10]]            [-2/5,  9/10]]
```

satisfy both interpolation constraints exactly.  Since the operator norm of an
affine matrix function is convex in the scalar argument, its maximum on
`[-1,0]` occurs at an endpoint.  Direct diagonalization gives

```text
||A||_op     = (11 + sqrt(41))/20,
||A-B||_op   = (1 + sqrt(61))/10,
sup ||P(x)|| = (1 + sqrt(61))/10 < 0.882.
```

Moreover

```text
[A,B] = [[0,-1/10],[1/10,0]],
```

so the coefficients do not commute.

### Commuting lower bound

Suppose instead that `A` and `B` are commuting real symmetric matrices.  In a
common orthonormal eigenbasis, `P` is diagonal with scalar affine entries
`p_r`.  Because `v_1` and `v_2` are neither orthogonal nor parallel, at least
one common eigenvector has nonzero overlap with both targets.  The corresponding
scalar affine entry obeys

```text
p_r(1/2)=p_r(1)=1,
```

hence `p_r` is identically one.  Therefore every commuting filter has stopband
norm at least one.  The rational construction gives the strict separation

```text
noncommuting optimum <= (1+sqrt(61))/10 < 1 <= commuting optimum.
```

This closes the logical possibility, left open in the earlier paper, that every
matrix-filter problem necessarily reduces to scalar Chebyshev filtering.

## 3. Numerical minimax check

The reproducible SDP in `research/route_c_affine_filter.py` discretizes the
stopband.  Convexity makes the endpoint constraints exact for degree one; the
interior grid is retained as an implementation check.  At a 45-degree second
target it obtains approximately

```text
optimal stopband norm  0.870840672
commutator norm        0.098647696.
```

The rational construction is deliberately slightly suboptimal so that every
claim can be replayed without trusting a solver.

## 4. Exact full-spark quadratic separation

A solver-independent non-affine witness now matches the mechanism of the
asymptotic theorem.  Take degree `N=2`, dimension `d=2`, four pass points and
targets

```text
lambda = (1/2, 1, 3/2, 2),
v      = ((1,-2), (1,-1), (1,1), (1,2)).
```

Every pair of targets is linearly independent, so the tuple is full spark.  Set
`P(x)=C_0+x C_1+x^2 C_2`, where

```text
C_0 = [[ 2/5,  -3/16],    C_1 = [[15/16, 3/20],
       [-3/16, 29/32]]           [ 3/20, 3/32]],

C_2 = [[-3/8,     0],
       [   0, -3/80]].
```

Direct rational substitution verifies all four tangential constraints.  The
coefficients are noncommuting; for example,

```text
[C_0,C_1] = [[0, 1053/12800],[-1053/12800,0]].
```

For a pairwise-commuting symmetric quadratic filter, every common eigenvector
overlaps at least `M-d+1=3=N+1` targets.  Its scalar entry is therefore one at
three distinct nodes and hence identically one.  The only commuting feasible
filter is `I`.

The noncommuting leakage is controlled on the full interval, not a grid.
After `x=t-1`, the leading principal minor and determinant of each of
`I-P(x)` and `I+P(x)` have strictly positive rational Bernstein coefficients on
`t in [0,1]`.  Exact determinant coefficients are

```text
det(I-P): 81/256, 1701/10240, 219/2560, 441/10240, 27/1280,
det(I+P): 53/1280, 8389/10240, 783/512, 21913/10240, 3371/1280.
```

Combining their determinant lower bounds with Bernstein upper bounds on the
traces gives

```text
I-P(x) >= (3/304) I,
I+P(x) >= (3/304) I,
sup ||P(x)||_op <= 301/304 < 1.
```

`verification/verify_quadratic_filter_witness.py` replays the complete
rational calculation.  This closes the earlier requirement for a certified
non-affine finite witness.

## 5. Scaling theorem candidate

For every sufficiently large integer `N`, set `d=3` and `M=N+3=d+N`.  There exist `M`
distinct pass points `lambda_j>0` and full-spark unit vectors `v_j in R^d` for
which the following two statements hold:

1. every degree-`N` matrix polynomial with pairwise commuting real symmetric
   coefficients and `P(lambda_j)v_j=v_j` is identically the identity; but
2. a degree-`N` polynomial with real symmetric, generally noncommuting
   coefficients obeys all the same constraints and has

```text
sup_{x in [-1,0]} ||P(x)||_op <= C exp(-cN).
```

Here `C,c>0` are independent of `N`.  Consequently the noncommuting Hermitian
class has an exponential advantage at a fixed three-channel dimension over
its commuting subclass.

This is now supported by a complete existence-proof skeleton below.  It remains
a **candidate theorem**, rather than a priority claim, until the proof is fully
written and the collision audit is complete.

### 5.1 Exact commuting obstruction

Let the targets be full spark, meaning that every `d` of them span `R^d`.  Any
nonzero vector can be orthogonal to at most `d-1` targets.  If the symmetric
coefficient matrices commute, choose their common orthonormal eigenbasis
`u_r`, and let `p_r` be the scalar degree-`N` polynomial on the `r`-th common
eigenspace.  Projecting `P(lambda_j)v_j=v_j` onto `u_r` gives

```text
(p_r(lambda_j)-1) <u_r,v_j> = 0.
```

Since `M=d+N`, at least `N+1` of these overlaps are nonzero.  The pass points
are distinct, so `p_r-1` has at least `N+1` roots and is zero identically.
This holds for every `r`; hence `P(x)=I` and its stopband norm is exactly one.

### 5.2 Exponentially small diagonal base point

Put `S=ceil((N+3)/3)` and split the `M=N+3` nodes among three coordinate
directions with group sizes `s_g in {S-1,S}`.  Use the rational interlaced grids

```text
lambda_(g,r) = 1 + 5(3r+g)/(6S),
g=0,1,2,  0<=r<s_g.
```

For the Lagrange cardinals `L_(g,r)` on one group and
`m_g=N-s_g+1`, set

```text
q_g(x) = sum_r L_(g,r)(x)
         T_(m_g)(2x+1)/T_(m_g)(2 lambda_(g,r)+1).
```

This has degree `N` and equals one on all nodes assigned to channel `g`.
Equispacing and the factorial denominators give

```text
sum_r |L_(g,r)(x)| <= (36e/5)^(s_g-1),  x in [-1,0].
```

Combining this with Chebyshev decay yields

```text
max_g sup[-1,0] |q_g| <= C0 exp(-cN),
c = (2 arcosh(3)-log(36e/5))/3 > 0.
```

The fixed diagonal polynomial `P_0=diag(q_0,q_1,q_2)` is the required base.

### 5.3 Symmetric interpolation survives a full-spark perturbation

Let `L_v` map the `N+1` real symmetric coefficient matrices to the stacked
vectors `(P(lambda_j)v_j)_j`.  At the coordinate target tuple, an off-diagonal
entry polynomial `p_rs=p_sr` is evaluated at the union of the nodes assigned to
coordinates `r` and `s`.  This union has at most
`2S<=N+1` nodes for `N>=7`, so degree-`N` evaluation is onto.  Diagonal entries
see at most `S` nodes.  Entry by entry, `L_v` is therefore surjective.

Surjectivity is open in finite dimension and supplies a continuous local right
inverse.  The perturbation can be made explicit.  Choose distinct rational
numbers `t_j in (0,1)`, put

```text
w_j = (1,t_j,t_j^2)^T,
v_j(epsilon) = normalize(e_(g(j)) + epsilon w_j),
epsilon = exp(-N^3),
```

where `g(j)` is the coordinate group at the base point.  For any three-element
subset, the unnormalized target determinant is a polynomial in `epsilon` with
rational coefficients.  Its leading coefficient is the nonzero Vandermonde
determinant of the corresponding `w_j`.  Since `exp(-N^3)` is transcendental,
none of these finitely many nonzero rational polynomials vanishes.  The
normalized target tuple is therefore full spark.

Choose the rational pass nodes as in the manuscript.  Their minimum separation
is at least `5/(6S)`.  With the maximum stop/pass operator norm on polynomials
and maximum Euclidean block norm on the data, entrywise Lagrange interpolation
gives

```text
||R_N|| <= 6S (27S/5)^(2S-1) = exp(O(N log N)).
```

The perturbation operator has norm at most `sqrt(3)`.  The base interpolant at
the pass nodes is at most

```text
2 (4e)^(S-1) exp(N DeltaEta).
```

A Neumann-series correction through the right inverse is valid because

```text
exp(-N^3) ||R_N|| = exp(-N^3+O(N log N)).
```

The manuscript gives explicit envelopes: the Neumann product is below
`2^(-288)` at `N=7`, and the correction divided by `exp(-cN)` is below
`2^(-217)` there; both envelopes decrease.  The corrected coefficients remain
real symmetric, and their correction is at most `exp(-cN)`.  Once
`(C0+1)exp(-cN)<1`, the commuting obstruction forces the family to be
noncommuting.

The perturbation is explicit but extremely small.  The remaining gates are:

1. verify the phenomenon is not already an immediate corollary of tangential
   interpolation or `H-infinity` control theory; and
2. move beyond the controlled Hermitian block benchmark to a naturally arising
   multichannel task that genuinely mandates multipoint pass-row geometry, and
   test adaptive/nonstationary scalar baselines without overstating physical
   consequences.

The earlier algebraic obstacles are now removed: the multi-node interpolant is
explicit and rationally located, the construction lies in the
real-symmetric/Hermitian coefficient class, and the full-spark perturbation is
deterministic.

### 5.4 Block-Krylov interpretation

For a Hermitian matrix `H` and starting block `B`, a block Krylov output has the
form

```text
sum_(k=0)^N H^k B C_k.
```

In an eigenbasis `H q_j=lambda_j q_j`, the row signature
`b_j^*=q_j^*B` is transformed to `b_j^* P(lambda_j)`, up to the harmless choice
of left-versus-right convention.  Tangential pass constraints preserve selected
spectral signatures, while `sup ||P(x)||` controls every unwanted signature in
the stopband.  Pairwise commuting symmetric `C_k` reduce, after one fixed
orthogonal channel rotation, to independent scalar filters.  The theorem
candidate therefore separates genuinely coupled block postprocessing from all
fixed-basis scalar filtering.  Turning this interpretation into an algorithmic
gain requires convergence and cost experiments; the algebraic separation alone
does not establish a faster eigensolver.

## 6. Scaling numerical stress test

`research/route_c_scaling_pilot.py` perturbs the coordinate tuple, checks every
`d`-target minor numerically, and solves the real-symmetric coefficient SDP.
For the first dimension in the surjective regime, `N=3`, `d=9`, `M=12`, seed 7
and perturbations `0.001`, `0.01`, and `0.03`, an 81-point optimization grid
gives respective objective values approximately `0.29835`, `0.31998`, and
`0.39771`, versus the exact commuting lower bound one whenever the numerical
full-spark check is trusted.  The coefficient commutator norms are respectively
about `0.467`, `0.391`, and `0.709`.  The script also evaluates the returned
polynomial on an independent 4001-point validation grid to expose between-node
overshoot.

For `N=2`, where symmetric surjectivity is not claimed, the strict advantage
persists but shrinks to approximately `0.99635`, `0.97016`, and `0.94834` over
the same perturbations.  This is exploratory evidence only: neither a
floating-point minor test nor a discretized stopband is a proof.

## 7. Literature collision warning

Bitangential matrix Nevanlinna--Pick interpolation supplies general matrix
Schur-class interpolation criteria, and matrix-valued polynomial interpolation
is classical.  Our interval-polynomial, fixed-degree, commuting-vs-general
separation is not identical to those results, but terminology-based searches
are insufficient.  Gate G1 requires comparison with the exact tangential
interpolation data after mapping the interval complement conformally to the
disk.

## 8. Exact reciprocal linear-phase FIR corollary

Set `x=(9 cos(omega)+5)/4`.  A real-symmetric degree-`N` polynomial becomes a
real-symmetric cosine polynomial and hence the zero-phase response of a causal
palindromic `2N+1`-tap MIMO FIR filter with common group delay `N`.  The theorem
maps `[-1,0]` to the high-frequency band
`[arccos(-5/9),pi]` and maps all pass points into cosine values `[-1/9,1]`.
The change of polynomial/cosine bases is invertible, so commuting taps are
equivalent to commuting original coefficients.  They are therefore forced to
the pure delay, whereas the coupled response has norm `C exp(-cN)`.

The rational quadratic witness uses `x=(3 cos(omega)+1)/2` and exactly the taps

```text
B0=B4=[[-27/128, 0], [0, -27/1280]]
B1=B3=[[27/64, 9/80], [9/80, 27/640]]
B2=[[113/320, -9/80], [-9/80, 577/640]].
```

Its pass cosine values are `0,1/3,2/3,1`, its stop band is
`[arccos(-1/3),pi]`, and its exact norm bound is `301/304`.  This closes a
natural fixed-latency interpretation, not the priority or deployed-application
gate.

## 9. Rational approximate-calibration margin

Normalize the four quadratic target directions.  A two-by-two Gram calculation
shows that no unit vector can have overlap below `4/25` with two different
targets.  Projecting an approximate tangential residual `delta` onto either
common eigenvector of a commuting comparator therefore controls the scalar
quadratic at three pass nodes by `(25/4)delta`.  Lagrange interpolation at zero
has worst absolute weight sum `17`, giving the exact lower bound

```text
rho(Q) >= 1 - (425/4) delta.
```

Thus the `301/304` noncommuting certificate remains strictly separated for
every `delta < 3/32300`.  The same statement transfers to the zero-phase
five-tap FIR response.  This is an exact theorem, not a numerical condition
estimate.
