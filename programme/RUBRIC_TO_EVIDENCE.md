# Rubric-to-evidence matrix: 5.31 to 7+

The critique of the superseded identifiability paper is converted here into
exit conditions.  More pages, more random seeds, and more finite-volume exact
diagonalization do not close any row by themselves.

| Blocking criticism | What the later Wishart paper already changes | Evidence required for the next paper | Status |
|---|---|---|---|
| Individual ingredients are classical | Adds an exact Loewner band, SDP alternative, power theorem, and rational replay | A theorem whose conclusion is not obtained by composing standard moment/localizer and concentration results | Open |
| Visibility is external | Still labels every conclusion as a visible gap | Derive a quantitative overlap/form-factor lower bound from the Hamiltonian, symmetry, and probe construction, uniformly in volume; or prove an equivalent noisy-data-to-bulk theorem that removes the external datum | Open, headline blocker |
| No thermodynamic result | Extends finite ANNNI numerics only to `L=16` | A theorem uniform in `L` with a controlled `L -> infinity` passage, or a certified bulk hierarchy with finite-sample coverage | Open, headline blocker |
| Chebyshev layer is not formalized | Replaces that layer in a different paper but does not formalize the old positive theorem | Machine-check or exact-certify the decisive new lemma, not merely an auxiliary identity | Open |
| Componentwise noise ignores correlations | Closed by the simultaneous Wishart--Loewner set | Retain correlated uncertainty and compare to modern covariance-aware baselines | Closed as infrastructure |
| Sign condition is only sufficient | Exact SDP compatibility is necessary and sufficient for the chosen finite relaxation | State clearly what remains relaxation-dependent; prove hierarchy convergence or a sharp alternative if this is part of the headline | Partly closed |
| Optimality only for scalar filters | Degree ablation does not address matrix-filter optimality | Either solve a genuinely matrix-valued minimax problem with a noncommutative advantage, or explicitly discard this route after collision/pilot audit | Open route, not mandatory if another route wins |
| Finite ANNNI example illustrates but does not discover physics | Larger statistical campaign still uses ANNNI as a benchmark | Produce a model-specific analytic statement: explicit parameter domain, uniform constant, or new scaling law verified beyond ED | Open |
| Numerical evidence could be solver-dependent | Representative exact-rational dual replay exists | Exact/interval proof must cover the new headline constant or lemma | Partly closed |
| “Minimal visibility datum” was overstated | Later paper avoids this wording | Characterize necessity only in an ordered class of assumptions, or make no minimality claim | Editorially closed |

## Non-negotiable 7+ evidence

The next manuscript must contain all of the following before Gate G4:

1. **Novelty delta:** a theorem-to-theorem comparison with at least three closest
   primary references, including hypotheses and conclusions.
2. **Hard theorem:** a proof with a nontrivial uniform or noncommutative step;
   an SDP formulation alone is insufficient.
3. **Physical closure:** the observable/probe must be linked to the Hamiltonian
   gap without an unexplained external visibility floor.
4. **Beyond finite ED:** either a thermodynamic proof or a certified hierarchy
   whose validity is uniform in volume.
5. **Adversarial evidence:** attempted counterexamples, a negative-results
   ledger, baseline failures, and reproduction from a clean environment.
6. **Artifact identity:** exact public hash included in every evaluation prompt.

