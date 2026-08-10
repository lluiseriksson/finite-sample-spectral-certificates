# Route B theorem pilot: calibration-robust bulk-gap exclusion

Status: proof architecture under Gate G1.  None of the statements below is yet
labelled a theorem of the project.  The purpose of this file is to expose the
logical load-bearing points before code or manuscript prose is built around
them.

## 1. Setting

Let `theta` lie in a compact parameter set `Theta` in `R^p`.  For each `theta`,
let `Phi_theta` be a translation-invariant finite-range interaction whose local
terms depend polynomially (affinely in the first implementation) on `theta`.
The interaction range and local Hilbert-space dimension are uniform over
`Theta`.

For a fixed proposed gap `gamma`, write `P_gamma(theta, omega)` when `omega` is
a KMS ground state for `Phi_theta` and satisfies the locally nondegenerate bulk
gap inequality

```text
omega(a* [H_theta,a]) >= gamma (omega(a* a) - |omega(a)|^2)
```

for every local observable `a`.  Define

```text
G_gamma = {theta in Theta : there exists omega with P_gamma(theta,omega)}.
```

The Xu et al. hierarchy supplies, for a fixed `theta`, nested finite SDP
feasibility tests `F_k(theta)` with

```text
theta in G_gamma  iff  F_k(theta) is feasible for every k.
```

Here `k` abbreviates the support and state-polynomial degrees.

## 2. Lemma pilot R1: finite rejection is open in Hamiltonian space

**Candidate statement.**  Fix a hierarchy level `k`.  Suppose all normalized
truncated state-polynomial variables at level `k` lie in a common compact set
`X_k`, and the constraint map is continuous in `(theta,x)`.  Then

```text
K_k = {theta in Theta : F_k(theta) is feasible}
```

is compact.  Consequently, if `F_k(theta_star)` is infeasible, there exists
`rho_k(theta_star)>0` such that `F_k(theta)` is infeasible whenever
`||theta-theta_star|| < rho_k(theta_star)`.

**Proof skeleton.**  The full finite-level feasible set is closed in the
compact product `Theta x X_k`.  Its projection onto `Theta` is compact and
hence closed.  The complement of `K_k` is open.

**Load-bearing check.**  Normalization and the archimedean state-polynomial
constraints must actually bound every auxiliary moment variable used in the
implementation.  Closedness fails if solver-only free variables or an
unbounded homogenization variable are introduced.

## 3. Lemma pilot R2: computable rejection radius

Topological openness alone gives no sample complexity.  Introduce a normalized
minimum-violation problem

```text
v_k(theta) = min over x in X_k of the smallest s >= 0 such that
             every equality has residual at most s and
             every PSD constraint becomes feasible after adding s I.
```

The exact normalization will use Frobenius/Euclidean residuals represented by
second-order or semidefinite cones.  It must satisfy

```text
v_k(theta)=0  iff  F_k(theta) is feasible.
```

If the coefficient maps of the level-`k` constraints are Lipschitz in `theta`
with a uniform bound over `X_k`, then

```text
|v_k(theta)-v_k(theta')| <= L_k ||theta-theta'||.
```

Therefore a certified lower bound `v_k(theta_star) >= m_k > 0` gives the
explicit robust radius

```text
rho_k >= m_k / L_k.
```

**Required proof work.**

1. choose a violation gauge for which zero is exactly feasibility, including
   complex Hermitian equality constraints;
2. prove a uniform bound on all truncated variables from the archimedean
   relations;
3. derive `L_k` from the interaction basis and that variable bound;
4. show how a rational primal/dual certificate proves `m_k>0` without trusting
   floating-point solver status.

This quantitative stability statement is a candidate decisive lemma.  A mere
appeal to continuity is insufficient for the paper.

## 4. The joint uncertainty--state hierarchy

Let `C` be a compact basic semialgebraic confidence set for `theta`.  The exact
robust question is

```text
does there exist theta in C and omega such that P_gamma(theta,omega)?
```

Adjoin commuting parameter variables to the state-polynomial algebra and add
the polynomial inequalities defining `C`.  The KMS and gap localizers contain
products of parameter monomials and state moments.  A Lasserre/state-polynomial
moment hierarchy then gives convex SDP outer relaxations of the joint problem.

### Candidate completeness statement R3

Assume the parameter quadratic module and the state-polynomial quadratic module
are archimedean.  The joint hierarchy is feasible at every degree if and only
if there exists a pair `(theta,omega)` in the exact robust set.

**Forward direction.**  A feasible pair induces every truncated moment
functional.

**Reverse direction under audit.**  Extract a limiting positive functional on
the tensor product of the commutative parameter algebra and the state-
polynomial algebra.  A representation theorem should produce a probability
measure supported on feasible pairs.  Nonempty support yields one pair.

**Mixture trap.**  A representation may instead produce a direct integral in
which state expectations and parameters are correlated.  That is acceptable
only if almost every fibre is a genuine state for the same fibre parameter and
satisfies the fibrewise KMS/gap constraints.  Constraints imposed merely after
averaging could allow cancellations and would not prove existence of a
fibrewise feasible pair.  The localizers must be multiplied by arbitrary
parameter test polynomials to force fibrewise validity.  This point is a
mandatory proof or a stop condition.

## 5. Statistical semidecision theorem pilot R4

Let `(C_t)` be a confidence sequence satisfying

```text
Pr_theta_star(theta_star in C_t for every t) >= 1-alpha
```

and `d_H(C_t,{theta_star}) -> 0` almost surely.  At arbitrary data-dependent
times, solve arbitrary levels of the joint hierarchy over `C_t`; stop only on
certified infeasibility.

**Validity candidate.**  With probability at least `1-alpha`, every reported
conclusion `Delta_bulk(theta_star) < gamma` is correct.  This needs no union
bound over times, hierarchy levels, gap thresholds selected predictably from a
countable grid, or solver witnesses, because on the simultaneous-coverage event
the true parameter belongs to every set searched.

**Almost-sure termination candidate.**  If
`Delta_bulk(theta_star) < gamma`, fixed-parameter completeness yields some
finite `k` for which `F_k(theta_star)` is infeasible.  Lemma R1 gives an open
rejection neighbourhood.  Eventually `C_t` lies inside it.  Exact emptiness of
the resulting compact joint semialgebraic feasibility set must then be detected
at a finite joint hierarchy degree.  A dovetailed schedule over `(t,k,degree)`
therefore halts almost surely.

**Quantitative corollary candidate.**  If R2 supplies radius `rho` and the
calibration confidence sequence satisfies a deterministic radius bound
`r_t(alpha)`, then the statistical part stops once `r_t(alpha)<rho`.  For a
sub-Gaussian affine calibration model this gives, up to time-uniform logarithms,

```text
t = O(sigma^2 (p + log(1/alpha)) / rho^2).
```

The computational degree required to certify the joint empty set remains a
separate quantity and must not be hidden inside the sample bound.

## 6. Falsification checklist

- Construct a finite-dimensional counterexample where averaged joint moments
  satisfy KMS/gap constraints although no fibre `(theta,omega)` does.  If the
  proposed localizers do not eliminate it, R3 fails.
- Add every parameter-polynomial multiplier required for fibrewise constraints
  and prove archimedeanity of the tensor-product module.
- Test whether the minimum-violation gauge can be zero at weakly infeasible SDP
  instances.  Compactness should rule this out only when the gauge and variable
  set are correctly normalized.
- Check whether Xu et al.'s state-polynomial framework already permits external
  commuting parameters with exactly this completeness statement.
- Check calibration relevance: the showcase must estimate genuinely uncertain
  local interaction coefficients, not pretend that a precisely programmed
  coupling is experimental shot noise.
- Separate exact statistical validity from numerical SDP reliability.

## 7. Gate-G1 criterion for Route B

Route B passes Gate G1 only if all of the following hold:

1. R2 is proved with an explicit and exactly certifiable radius;
2. R3 survives the mixture trap;
3. the theorem-to-theorem literature audit finds no prior statistical
   semidecision result for uncertain local interactions; and
4. a small nonintegrable-chain pilot exhibits a nonzero robust radius that is
   not destroyed by the first realistic calibration confidence set.

