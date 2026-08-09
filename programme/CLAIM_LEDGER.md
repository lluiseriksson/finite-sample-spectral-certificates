# Claim ledger

Every headline sentence starts as `candidate`.  It becomes `supported` only
after a proof, collision audit and linked artifact exist.  A failed claim stays
in this ledger.

| ID | Candidate claim | Current status | Evidence needed / reason forbidden |
|---|---|---|---|
| A0 | Weak perturbations of a product Hamiltonian produce a visible one-particle band. | **Forbidden as novel** | Yarotsky (2004) already proves a one-particle subspace and exponentially localized projected bare excitations. |
| A1 | A specified bare ANNNI momentum/parity probe has an explicit positive one-particle residue uniform in volume on a nonzero high-field domain. | Candidate, high risk | Model-explicit constants, thermodynamic passage, noncircular proof and collision check against strong-field Ising/cluster expansions. |
| A2 | Finite correlator data alone prove gappedness in the high-field phase. | **Forbidden** | The phase-isolation input would already be an independent gap theorem; data may sharpen a mass inside a proven window but cannot be credited with the phase proof. |
| B0 | Finite-shot confidence sets can be inserted into many-body moment SDPs. | **Forbidden as novel** | Mortimer et al. (2026) already do this. |
| B1 | A complete SDP hierarchy exists for upper bounds on the thermodynamic bulk gap. | **Forbidden as novel** | Xu et al. (2026) already prove completeness. |
| B2 | A joint compact-parameter/state-polynomial hierarchy is complete for the existence of a Hamiltonian inside a semialgebraic uncertainty set with a KMS ground state of bulk gap at least `gamma`. | Candidate | Formal hierarchy, archimedean compactness proof, both directions, comparison to Xu's fixed-H theorem. |
| B3 | Shrinking anytime-valid confidence sets plus the joint hierarchy give a statistically valid and almost-surely terminating semidecision whenever `gamma` is strictly above the true bulk gap. | Candidate, current lead | Need the finite-level robust-neighbourhood lemma, stopping-time coverage proof, and an explicit calibration model/sample bound. |
| B4 | The method proves a positive Yang--Mills mass gap. | **Forbidden** | The proposed result concerns lattice Hamiltonian bulk-gap exclusion/upper bounds, not constructive four-dimensional Yang--Mills. |
| C0 | Matrix Chebyshev polynomials or block Krylov methods are new. | **Forbidden** | Both are established bodies of literature. |
| C1 | Tangential degree-`N` matrix-polynomial filtering admits an exponential-in-`N=sqrt(d)` separation between real-symmetric noncommuting coefficients and pairwise-commuting real-symmetric coefficients. | Candidate, quantified theorem proof and exact finite certificate found | Full-spark quadratic witness is exactly replayed with continuum bound `301/304<1`; the asymptotic proof now fixes the norms and gives the correction bound `64152 N^9 exp(-N^3+N DeltaEta)`.  Still needs the completed tangential-interpolation collision audit and an application with genuine scientific payoff. |

## Current headline under test

> The absence of a proposed thermodynamic bulk gap is statistically
> semi-decidable even when the local interaction is known only through
> sequential calibration data.

This wording remains provisional.  It must be weakened if the hierarchy only
handles a restricted parameter class or if termination needs an additional
dual-margin assumption.

## Explicit nonclaims

- No solution of the Yang--Mills mass-gap problem.
- No lower bound proving gappedness for a generic interacting Hamiltonian.
- No first SDP hierarchy for spectral gaps.
- No claim that optional stopping or confidence sequences are new.
- No claim that numerical SDP solver output is a proof without rational or
  interval replay.
