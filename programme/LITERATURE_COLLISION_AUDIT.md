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
| P. A. Fuhrmann, [*On tangential matrix interpolation*](https://doi.org/10.1016/j.laa.2010.07.018) (2010) | Algebraic polynomial and rational matrix interpolation with tangential Lagrange/Hermite/Newton data | **Direct collision:** existence and parametrization of matrix polynomials satisfying `P(lambda_j)v_j=w_j` are classical.  Novelty can only concern the fixed-degree stopband minimax separation and its spectral-certificate role. |
| J. A. Ball and J. Kang, [*Matrix polynomial solutions of tangential Lagrange--Sylvester interpolation conditions of low McMillan degree*](https://doi.org/10.1016/0024-3795(90)90145-3) (1990) | Explicit realization-form matrix-polynomial solutions for tangential Lagrange--Sylvester data with low McMillan degree | **Direct collision:** low-complexity tangential polynomial construction predates the proposed work.  The candidate theorem must be a norm separation under a fixed ordinary degree, not an existence or low-degree interpolation claim. |
| J. A. Ball and V. Bolotnikov, [*The bitangential matrix Nevanlinna--Pick interpolation problem revisited*](https://arxiv.org/abs/1611.07097) (2016) | Schur-class matrix interpolants, Pick-matrix criteria and linear-fractional parametrization for bitangential data | A minimum-supremum-norm tangential interpolant without a fixed polynomial-degree restriction belongs to established Nevanlinna--Pick theory. |
| A. Blomqvist, A. Lindquist and R. Nagamune, [*Matrix-valued Nevanlinna--Pick interpolation with complexity constraint: an optimization approach*](https://doi.org/10.1109/TAC.2003.820227) (2003) | Matrix-valued analytic interpolation with bounded complexity, convex parametrization and robust-control applications | Collides with any broad claim that combining a matrix norm, interpolation and a complexity constraint is new.  Its Schur-class/rational setting is not the same as the fixed ordinary-degree interval-polynomial commuting obstruction. |
| Y. Kuroiwa and A. Lindquist, [*Bi-Tangential Nevanlinna--Pick Interpolation with a Complexity Constraint*](https://mathweb.ucsd.edu/~helton/MTNSHISTORY/CONTENTS/2006KYOTO/CONFERENCEWEBSITE/papers/0289.pdf) (2006) | Convex construction of bitangential rational matrix interpolants under McMillan-degree bounds; lower-degree sensitivity-shaping controllers | **Direct degree-constraint collision:** bounded-complexity bitangential matrix interpolation is established.  The remaining distinction is ordinary polynomial degree on a real stopband plus the exact commuting-subclass obstruction. |
| J. Stefanovski and D. Georgijević, [*Interpolation with constraint on frequency region and systems & control application*](https://doi.org/10.1016/j.sysconle.2016.08.009) (2016) | Bitangential interpolation by real stable rational matrices with arbitrarily small spectral norm on a proper frequency region; MIMO filtering and control applications | **Closest norm-region collision found so far.**  Small stopband norm under tangential matrix constraints is established.  The surviving delta is the fixed ordinary polynomial degree, Hermitian coefficients, explicit `d=N^2` tradeoff, and exact obstruction for the commuting subclass. |
| K.-C. Toh and L. N. Trefethen, [*The Chebyshev polynomials of a matrix*](https://doi.org/10.1137/S0895479896303739) (1998) | Monic scalar polynomials minimizing `||p(A)||_2`, computed by SDP, with Krylov motivation | Minimax polynomial filtering and SDP formulations are classical.  This work does not use matrix coefficients or tangential pass constraints. |
| M. Rinelli and R. Vandebril, [*Block Krylov subspaces and orthogonal matrix polynomials: a structural correspondence with applications to unitary matrices*](https://arxiv.org/abs/2605.16954) (2026) | Isometric correspondence between block Krylov spaces and matrix-polynomial spaces; spectral matrix measures and short recurrences | Supplies the natural application language and prevents claiming the block-Krylov/matrix-polynomial correspondence.  Its stated theorems concern orthogonality and recurrences, not a commuting-vs-noncommuting stopband minimax gap. |
| M. Hartz, [*On von Neumann's inequality on the polydisc*](https://doi.org/10.1007/s00208-024-03040-2) (2025), Sec. 2 | Sharp and near-sharp norm inequalities for one-variable polynomials with operator coefficients, including distinct commuting and noncommuting regimes | Confirms that separations caused by operator coefficients are established phenomena.  Our exact interpolation-constrained separation must not be advertised as the first benefit of noncommuting coefficients. |

The search has not yet found a primary source proving or disproving the precise
fixed-degree, real-interval, tangential interpolation separation in
`ROUTE_C_THEOREM_PILOT.md`.  Stefanovski--Georgijević comes especially close at
the problem level, but its headline theorem permits stable rational matrices of
unrestricted growing degree and does not state a commuting-coefficient lower
bound.  The surrounding ingredients are emphatically classical, so this
remaining difference is not yet evidence of novelty.

### Required pilot before a novelty claim

The first part is met in dimension two: an exact rational affine filter has
stopband norm below `0.882`, whereas every commuting symmetric affine filter
has norm at least one.  A candidate asymptotic theorem now gives leakage
`C exp(-cN)` for noncommuting symmetric coefficients versus exactly one for the
commuting class at dimension `d=N^2`.  The proof skeleton includes an explicit
full-spark perturbation and a block-Krylov interpretation.  The remaining G1
bar is a line-by-line norm-constrained interpolation collision check and a
fully quantified proof; the two-dimensional example alone is not a paper-scale
contribution.

## Provisional Gate-G1 ranking

1. **Route C:** now has an exact finite witness, an asymptotic exponential
   separation proof skeleton in the symmetric/Hermitian class, reproducible
   SDP stress tests, and a block-Krylov interpretation.  Novelty collision and
   quantitative proof closure remain mandatory before G1 passes.
2. **Route B, reformulated as joint uncertainty--state polynomial
   completeness:** strongest thermodynamic architecture, but it needs a real
   compactness/robust-infeasibility proof and must beat Mortimer's scope.
3. **Route A, explicit bare-probe residue in nonintegrable ANNNI:** strongest
   physical payoff, but Yarotsky already supplies the qualitative mechanism
   and circularity is a major risk.

No route has passed Gate G1 yet, but Route C is the current lead.
