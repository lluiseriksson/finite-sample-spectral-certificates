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

## 4. Scaling theorem under investigation

Let `M` distinct pass points and target vectors in `R^d` satisfy the full-spark
condition that every `d` targets span `R^d`.  If `M >= d+N`, then every common
eigenvector of commuting normal coefficient matrices overlaps at least
`M-d+1 >= N+1` target vectors.  A scalar degree-`N` entry satisfying the
corresponding tangential constraints is therefore identically one.  Thus every
simultaneously unitarily diagonalizable degree-`N` filter has stopband norm at
least one.

The constructive side starts at a rank-deficient but exactly solvable data set.
Take `d=N^2` coordinate target directions and `M=d+N` distinct pass points.
Assign two pass points to each of the first `N` coordinate directions and one
pass point to every remaining direction.  For a direction with one pass point,
a normalized degree-`N` Chebyshev polynomial has exponentially small stopband
norm.  For a direction with two pass points `y_1,y_2>1` after mapping the
stopband to `[-1,1]`, solve

```text
c_N T_N(y_i) + c_(N-1) T_(N-1)(y_i) = 1,  i=1,2.
```

The determinant is nonzero for distinct positive pass points.  Hyperbolic
cosine asymptotics give

```text
|c_N|+|c_(N-1)| <= poly(N) exp(-N eta_min),
eta_min = min_i arcosh(y_i)>0,
```

when the pass-point separation is only polynomially small.  Placing these
scalar polynomials on the diagonal produces a degree-`N` matrix polynomial with
stopband leakage `poly(N) exp(-N eta_min)`.

At the coordinate data, the tangential interpolation map is surjective: after
splitting by input coordinate it is a direct sum of Vandermonde evaluation maps
with at most two distinct nodes, well below the `N+1` available coefficients.
Surjectivity is open.  Full-spark target tuples are dense, so the vectors can be
perturbed by an arbitrarily small amount to become full spark while retaining a
nearby interpolant and essentially the same stopband norm.  The commuting lower
bound then becomes one.

If every stability estimate closes, this yields an existence separation of the
form

```text
general matrix optimum <= poly(N) exp(-c N),
simultaneously diagonalizable optimum = 1,
dimension d = N^2.
```

Thus the advantage is exponential in `sqrt(d)`.  This argument is presently a
proof sketch, not a theorem.  The unresolved points are:

1. an explicit right-inverse bound for the perturbed interpolation operator;
2. an explicit full-spark perturbation small enough to preserve the exponential
   leakage, rather than a purely existential density argument;
3. extension from arbitrary real coefficient matrices to a physically natural
   Hermitian or adjoint-paired filter class; and
4. whether the separation is already a corollary of tangential interpolation or
   `H-infinity` control theory.

Without these four items, the 2 by 2 example is a useful lemma but not a 7+
paper.

## 5. Literature collision warning

Bitangential matrix Nevanlinna--Pick interpolation supplies general matrix
Schur-class interpolation criteria, and matrix-valued polynomial interpolation
is classical.  Our interval-polynomial, fixed-degree, commuting-vs-general
separation is not identical to those results, but terminology-based searches
are insufficient.  Gate G1 requires comparison with the exact tangential
interpolation data after mapping the interval complement conformally to the
disk.
