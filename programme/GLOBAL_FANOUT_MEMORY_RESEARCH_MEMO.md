# Global fan-out memory: research memo

## Target theorem

Let `L >= 2`, let `zeta_1,...,zeta_L` be distinct boundary nodes, let `X` be an
isometric `N x k` input frame, and let `Y_1,...,Y_L` span mutually
orthogonal `k`-planes (`N >= L k`).  Among square finite rational-inner
matrices regular at the nodes and satisfying

```text
ran S(zeta_i) X = ran Y_i,
```

the exact minimum McMillan degree is

```text
k (L - 1).
```

Every two-node restriction has minimum degree `k`, so the ratio between
the global cost and the largest pairwise cost is `L-1`.

## Proof skeleton

For calibrated target frames `W_i`, every positive boundary-Pick completion
`P` has off-diagonal blocks

```text
P_ij = (I - W_i^* W_j) / (1 - conj(zeta_i) zeta_j).
```

With `Z = diag(zeta_i I_k)` and the block Gram matrix
`G = [W_i^* W_j]`, this gives the Stein identity

```text
P - Z^* P Z = (11^*) tensor I_k - G.
```

If `P` is positive semidefinite, the negative inertia of the right-hand
side is at most `rank P`.  Therefore every interpolant obeys the general
calibrated-frame lower bound

```text
mcdeg S >= nu_-((11^*) tensor I_k - G).
```

It also obeys the gauge-independent span bound

```text
mcdeg S >= rank G - k
          = dim span(Y_1,...,Y_L) - k.
```

For mutually orthogonal target planes, `G = I_(Lk)`, so the Stein right-hand
side has `k(L-1)` negative eigenvalues.

For sufficiency, set

```text
C_ij = 1 / (1 - conj(zeta_i) zeta_j), i != j,
C_ii = 0,
P_0 = (C - lambda_min(C) I_L) tensor I_k.
```

Then `P_0` is a feasible positive completion and has rank at most
`k(L-1)`.  The lower bound forces equality.  The boundary lurking-isometry
construction produces a rational-inner interpolant of that degree.

For regular-polygon nodes, the eigenvalues of `C` are

```text
-(L-1)/2, -(L-3)/2, ..., (L-3)/2, (L-1)/2,
```

so the canonical completion has spectrum `0,1,...,L-1`, each repeated
`k` times.

## Robust theorem

For any measured Hermitian approximation `H_hat` to the Stein right-hand
side with `||H-H_hat|| <= eta`, the fail-closed certificate is

```text
mcdeg S >= number of eigenvalues of H_hat below -eta.
```

If `G` is the target-frame Gram matrix and `||G-I|| < 1`, the full
orthogonal lower bound `k(L-1)` survives.  A simpler sufficient condition is

```text
(L-1) max_{i != j} ||Y_i^* Y_j|| < 1.
```

This norm condition is uniform over all frame gauges.  The noisy inertia
count itself uses coherently calibrated frames; subspace-only data require
minimization over endpoint gauges.

## Action comparison

Order the nodes cyclically.  Orthogonality makes all `k` principal angles on
each adjacent arc equal to `pi/2`, so endpoint-only positive-delay geometry gives

```text
total trace action >= L k pi,
mcdeg S >= ceil(L k / 2).
```

The exact global degree `k(L-1)` is strictly stronger for `L >= 3` and is
asymptotically twice this endpoint-bound count.  Measuring the full
Wigner--Smith action still recovers the exact degree.  At minimum degree, the total
winding action exceeds the sum of the local geodesic minima by
`pi k (L-2)`.

## Falsifiers and kill criteria

The programme fails if any of the following occurs:

1. a primary source already states the exact orthogonal fan-out law in the
   same boundary rational-inner/McMillan setting;
2. the boundary Pick rank theorem does not apply to subspace interpolation
   after quotienting endpoint frames;
3. a numerical lurking-isometry construction at rank `k(L-1)` fails to be
   regular, inner, interpolating, controllable, or observable;
4. the robust inertia certificate produces a false positive under its stated
   operator-norm error model.

## Required artifact

- analytic proof of the inertia and construction theorems;
- independent numerical implementations of the Pick completion and inertia
  certificate;
- arbitrary-node and regular-polygon campaigns over multiple `L` and `k`;
- explicit colligations for representative multi-channel fan-out instances;
- noisy fail-closed tests and negative controls;
- a theorem-level comparison with the eight preceding papers and the closest
  interpolation literature.
