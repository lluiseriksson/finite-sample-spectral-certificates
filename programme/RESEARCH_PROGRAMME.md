# Ten-week frontier programme

## Decision target

The programme will not extend the existing paper by accumulation.  It will
select one of three mutually competing headline contributions:

### Route A — endogenous visibility in a concrete phase

Target a translation- and parity-resolved momentum probe in the high-field
paramagnetic phase of the one-dimensional quantum ANNNI chain.  At the
decoupled point the probe creates the one-particle excitation exactly.  The
research question is whether spectral perturbation or quasi-local spectral
flow yields an **explicit quasiparticle-residue lower bound uniform in system
size**, on a nonzero parameter domain.  Combined with correlator confidence
sets, this would turn a visible-edge exclusion into a Hamiltonian-gap statement
without inserting an empirical overlap floor.

Main danger: a gap-stability theorem used to construct the spectral flow may
already assume enough gap information to make the inference circular.  The
proof must separate a reference-phase stability assumption from the unknown
numerical value being certified.

### Route B — calibration-robust bulk-gap semidecision (current priority)

Let the finite-range interaction itself be known only through a shrinking
simultaneous confidence set for its local coefficients.  Extend the
state-polynomial hierarchy to optimize jointly over the interaction parameters
and the KMS state.  The desired theorem is a probably-valid, anytime-safe
semidecision procedure whose coverage survives adaptive hierarchy depth and
data-dependent stopping and whose termination is almost sure whenever the
proposed gap is strictly above the true bulk gap.

The nontrivial step is a robust completeness theorem: if the singleton
Hamiltonian is rejected at some finite hierarchy level, compactness must imply
that a full neighbourhood of Hamiltonian coefficients is rejected at that same
level.  A shrinking confidence sequence then eventually enters that
neighbourhood.  This would convert algebraic semidecidability into statistical
semidecidability for a calibrated-but-uncertain many-body interaction.

Main danger: the parameter extension may follow too directly from existing
state-polynomial Positivstellensatz machinery, or the physical calibration
model may be too artificial.  Gate G1 requires a precise theorem and a
theorem-to-theorem delta against both Xu et al. and Mortimer et al.; optional
stopping alone is not novel.

### Route C — noncommutative minimax filters (high-risk alternative)

Replace scalar Chebyshev filters with matrix-valued polynomials for
multichannel spectral measures.  Seek a family in which the minimax optimum is
provably noncommuting and beats every scalar or rank-one construction by a
dimension-dependent factor.

Main danger: operator-valued Chebyshev and matrix approximation theory may
already contain the desired result, or the physical interpretation may remain
too methodological for a 7+ paper.

## Weekly outputs and gates

| Week | Output | Gate |
|---:|---|---|
| 1 | Artifact manifest, rubric matrix, clean evaluation packet, initial claim ledger | Evaluation must name the correct SHA |
| 2 | Primary-literature map and theorem collision table for A/B/C | At least three closest works per route |
| 3 | Proof/computation pilots and counterexample search | **G1:** one route has a nonclassical lemma and a precise novelty delta |
| 4 | Full proof skeleton and constants for the selected route | No hidden finite-volume or visibility assumption |
| 5 | Uniform/thermodynamic closure | **G2:** uniform constants close; otherwise pivot to Route B |
| 6 | ED/DMRG/series or hierarchy campaign with modern baselines | Evidence tests the theorem, not merely solver performance |
| 7 | Interval/rational certificate and Lean formalization of the decisive lemma | **G3:** headline evidence survives independent replay |
| 8 | New manuscript and complete claim-to-artifact links | No inherited text without a current role |
| 9 | Blind adversarial review, counterexample audit, clean reproduction | **G4:** score 7+ or concrete blocking report |
| 10 | Repair remaining blockers, freeze public hashes, CI, data, proof and paper | Publish only after all gates pass |

## Stop rules

- Route A stops if the only overlap lower bound scales to zero for every
  physically local or experimentally aggregate probe considered, or if its
  proof is circular in the target gap.
- Route B stops if finite-sample validity follows by a one-line union bound and
  no new robustness/convergence theorem is needed.
- Route C stops if the optimum always reduces to a commuting scalar problem in
  the tested operator-valued classes, or if the closest literature already
  supplies the claimed advantage.
- A failed route is recorded with code/counterexample and is not repackaged as
  a positive headline.
