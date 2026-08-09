# Primary-literature collision audit

Status: live audit for Gate G1.  A citation in this file is evidence against an
overbroad novelty claim, not evidence that the proposed replacement claim is
new.  Claims remain forbidden until the theorem statements and proofs of the
closest works have been compared line by line.

## Route A: endogenous one-particle visibility

### Established results that collide with the naive claim

| Work | Hypotheses | Established conclusion | Collision |
|---|---|---|---|
| D. A. Yarotsky, [*Quasi-particles in weak perturbations of non-interacting quantum lattice systems*](https://arxiv.org/abs/math-ph/0411042) (2004) | Translation-invariant weak finite-range perturbation of a decoupled lattice Hamiltonian; unique gapped product ground state; isolated nondegenerate on-site excitation | The isolated band persists as a one-particle subspace with real-analytic dispersion.  More sharply, the spectral projection of a local bare excitation has the expansion in Eqs. (19)--(20), and the orthonormal one-particle basis has the exponentially localized expansion in Eqs. (23)--(24). | **Fatal to the claim that weak interaction alone newly implies a visible quasiparticle.**  The paper already contains a qualitative uniform overlap/localization mechanism, although its small constants are not made model-explicit. |
| S. Bachmann, W. Dybalski, P. Naaijkens, [*Lieb--Robinson bounds, Arveson spectrum and Haag--Ruelle scattering theory for gapped quantum spin systems*](https://arxiv.org/abs/1412.2970) (2016) | Translation-invariant gapped system with an isolated regular or pseudo-relativistic mass shell | Almost-local observables create single-particle vectors; Haag--Ruelle scattering states and an S-matrix are constructed.  The strong-field Ising regime is included among the rigorous examples. | **Fatal to presenting almost-local creation operators or mass-shell filtering as new.** |
| J. Haegeman et al., [*Elementary Excitations in Gapped Quantum Spin Systems*](https://arxiv.org/abs/1305.2176) (2013) | Unique gapped ground state, local Hamiltonian | A local-operator/tangent-space ansatz approximates isolated elementary excitations with controlled error. | Collides with any claim that local probes approximate isolated excitations merely from locality and a gap. |
| A. Ukai, [*Volume-Independent Spectral Stability of Energy-Truncated Effective Hamiltonians in Quantum Spin Systems*](https://arxiv.org/abs/2605.07410) (2026) | Bounded finite-range interactions, fixed target region | Volume-uniform spectral-overlap bounds in finite volume and the GNS thermodynamic limit. | Collides with generic claims of first volume-uniform spectral-overlap control; a new theorem must exploit the specified probe/form factor, not just low-energy projection stability. |

### Surviving candidate claim

For a fixed, experimentally meaningful translation/parity-resolved observable
in a concrete nonintegrable chain, derive a **numerical, computable lower bound**
on its one-particle spectral weight which:

1. is uniform in volume and survives the thermodynamic limit;
2. has an explicit nonempty coupling domain;
3. uses the bare observable, not a spectrally dressed operator whose
   construction presupposes the target band;
4. separates the reference-phase isolation estimate from the unknown numerical
   mass subsequently inferred from data; and
5. gives constants strong enough to close a finite-sample certificate.

This is not yet a novelty claim.  Yarotsky's Eqs. (19)--(24) make the candidate
plausible but also set the minimum bar: replacing an unspecified `epsilon` by an
explicit model-dependent constant must require a genuinely new quantitative
argument, not bookkeeping.

### Circularity test

If the proof invokes a spectral flow whose construction assumes a gap at least
as strong as the conclusion to be certified, Route A fails.  A permissible
theorem may use a coarse, independently proved reference-phase isolation window
and use data only to sharpen the excitation energy inside that window.  The
manuscript must then say that it calibrates the mass within a proven phase; it
must not claim to prove gappedness from the same data.

## Route B: finite statistics and thermodynamic hierarchies

### Established results that collide with the naive claim

| Work | Hypotheses | Established conclusion | Collision |
|---|---|---|---|
| X. Xu et al., [*The bulk spectral gap is semi-decidable: a convergent family of certified upper bounds*](https://arxiv.org/abs/2606.03836) (2026) | Known finite-range interaction and the quasi-local operator algebra | A complete state-polynomial SDP hierarchy: feasibility at all levels is equivalent to existence of a KMS ground state with locally nondegenerate bulk gap at least the proposed value. | **Fatal to any claim of first convergent or first certified thermodynamic upper-gap hierarchy.** |
| L. Mortimer et al., [*Bounding many-body properties under partial information and finite measurement statistics*](https://arxiv.org/abs/2601.10408) (2026) | Partial finite-shot measurements plus moment/RDM relaxations and optional ground-state, symmetry or steady-state constraints | Probabilistic many-body bounds from measurement confidence regions embedded in scalable moment SDPs. | **Fatal to the generic claim that finite-shot confidence regions can newly be inserted into many-body moment relaxations.** |
| S. Rao, [*A convergent hierarchy of spectral gap certificates for qubit Hamiltonians*](https://arxiv.org/abs/2510.08427) (2025) | Finite qubit Hamiltonian, universal-enveloping-algebra/NPA constraints and a ground-energy upper bound | Convergent lower-gap certificate hierarchy. | Prevents conflating the upper-bulk problem with general lower-gap certification. |
| K. S. Rai et al., [*A Hierarchy of Spectral Gap Certificates for Frustration-Free Spin Systems*](https://doi.org/10.22331/q-2026-04-13-2065) (2026) | Translation-invariant frustration-free spin systems | Thermodynamic lower-gap SDP hierarchy subsuming finite-size criteria. | Collides with generic thermodynamic lower-bound claims in the frustration-free setting. |
| S. R. Howard et al., [*Time-uniform, nonparametric, nonasymptotic confidence sequences*](https://arxiv.org/abs/1810.08240) (2021) | Sequential sub-Gaussian/Bernstein/self-normalized or matrix-martingale observations | Confidence sequences valid under optional stopping, including covariance-matrix applications. | Makes `anytime-safe by a confidence sequence' a standard statistical component, not the headline theorem. |

### Surviving candidate claim

The only currently credible Route-B upgrade is a **joint uncertainty--state
polynomial hierarchy**.  Its input is a shrinking, simultaneous confidence set
for either local interaction coefficients or a growing collection of local
state moments.  The target theorem would prove all three statements:

1. finite-level infeasibility is valid uniformly over the entire confidence
   set, even after adaptive selection of hierarchy level and witness;
2. feasibility at every algebraic level is equivalent to the existence of a
   coefficient/state pair in the limiting uncertainty set satisfying the KMS
   and bulk-gap inequalities; and
3. under a separating gap margin and shrinking confidence sets, a diagonal
   schedule of shots and hierarchy depth terminates almost surely, with an
   explicit finite-sample bound when a finite-level dual margin exists.

Items 2--3 are the required nonclassical work.  A union bound around the Xu
hierarchy or the Mortimer finite-shot constraints fails Gate G1.

## Route C: matrix-valued minimax filters

### Established adjacent theory

| Work/family | Established conclusion | Collision risk |
|---|---|---|
| Matrix orthogonal-polynomial theory; e.g. S. Delvaux and H. Dette, [*Zeros and ratio asymptotics for matrix orthogonal polynomials*](https://arxiv.org/abs/1108.5155) | Matrix recurrences, matrix Chebyshev measures, zero and ratio asymptotics | `Matrix Chebyshev polynomial' is established terminology and cannot itself be a novelty claim. |
| Block Lanczos and block Krylov literature | Matrix-valued recurrences and multiple starting vectors are standard | A multichannel implementation or faster convergence example is insufficient. |
| Classical constrained minimax approximation | Alternation and extremal scalar filters are classical | A commuting matrix coefficient solution that diagonalizes into scalar problems is insufficient. |

The search has not yet found a primary source proving or disproving the precise
noncommutative minimax advantage proposed here.  This absence is not evidence of
novelty; Route C remains high risk until a sharper search by problem statement,
not terminology, is complete.

### Required pilot before a novelty claim

Construct the smallest PSD matrix-valued spectral measure for which an SDP over
matrix-polynomial coefficients has a unique optimizer with noncommuting
coefficients.  Prove a lower bound for every commuting/scalar filter and an
upper construction with a dimension-dependent strict gap.  Then search the
operator approximation, robust control and block-Krylov literatures using the
exact extremal formulation.  If the optimizer always jointly diagonalizes in
the tested classes, Route C stops.

## Provisional Gate-G1 ranking

1. **Route B, reformulated as joint uncertainty--state polynomial
   completeness:** strongest surviving theorem architecture, but it needs a
   real compactness/robust-infeasibility proof and must beat Mortimer's scope.
2. **Route A, explicit bare-probe residue in nonintegrable ANNNI:** strongest
   physical payoff, but Yarotsky already supplies the qualitative mechanism
   and circularity is a major risk.
3. **Route C:** least collided by the current search, but its physical payoff
   and existing approximation-theory coverage are both uncertain.

No route has passed Gate G1 yet.
