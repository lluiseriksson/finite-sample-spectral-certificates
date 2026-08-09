# Long-horizon programme: model-derived visibility and certified spectroscopy

Status: exploratory; no paper claim has been selected yet.

## Why the 5.31 paper is not the target

The current paper establishes a clean information boundary, but its positive
branch imports visibility and its interacting experiment stops at finite-size
illustration.  More examples of the same theorem will not cross the target.

## Exit criteria for a 7+ manuscript

A candidate advances to manuscript only if it supplies all of the following:

1. A model-specific producer of visibility, or a theorem showing precisely
   how visibility scales with subsystem size, volume, symmetry sector, and
   proximity to criticality.
2. A result not reducible to a convex mixture, standard scalar Chebyshev
   extremality, or a known spectrum inserted into the method.
3. At least two independent computational routes, with one withheld as ground
   truth and one using correlator data only.
4. A scaling campaign beyond exact-diagonalization demonstration sizes,
   including finite-size, noise, covariance, and operator-family ablations.
5. Comparison against GEVP/Ritz, exponential fitting, scalar filtering, and a
   correlated-noise baseline.
6. A negative-results ledger recording failed conjectures and parameter
   regimes where the proposed certificate loses power.
7. Fixed-commit artifacts, deterministic small checks, statistical replay,
   and formal or interval certification of the mathematical core where viable.
8. Explicit separation from the two earlier spectral papers and an honest
   thermodynamic/continuum claim boundary.

## Candidate A: local visibility length in interacting spin chains

Define the transition operator

    X_A = Tr_{A^c} |psi_0><psi_1|

and the optimal normalized A-local transition strength.  An orthonormal local
operator basis converts the sum of all squared matrix elements into the
Hilbert--Schmidt norm of X_A.  The research question is whether this quantity
has a universal scaling law in gapped, critical, symmetry-broken, edge-mode,
and frustrated phases, and whether it gives a model-derived gamma for finite
Euclidean gap certificates.

Pilot models: TFIM (analytic control), ANNNI (nonintegrable/frustrated), cluster
or SSH-type edge phase (designed local invisibility), and disordered controls.
Scale route: exact diagonalization for audit, then MPS/DMRG contractions.

Kill criterion: discard as a main paper if the result remains only the elementary
operator-basis identity plus numerical phase plots.

## Candidate B: covariance-optimal multichannel moment certificates

Replace componentwise error bars and scalar filters by a joint confidence
ellipsoid for all matrix moments.  Derive the dual robust localizer problem,
matrix-polynomial filter, exact type-I error control, and an operator-selection
criterion.  Apply it to synthetic spin data and genuine lattice-gauge Monte
Carlo correlator matrices.

Kill criterion: discard if the optimization is a routine SDP reformulation
without a sharp theorem, demonstrable sample-complexity gain, or new physical
conclusion.

## Candidate C: gauge-theory visibility from exact 2D to Monte Carlo 4D

Use the exact SU(2) heat-kernel/character transfer spectrum to derive explicit
Wilson-channel weights and a complete calibration case.  Then test which
features survive for smeared Wilson-loop/glueball operator bases in 2+1D or 4D
Wilson-action ensembles, with correlated bootstrap and volume/coupling scans.

Kill criterion: the exact 2D part cannot be the headline if it merely rewrites
known character expansion.  The higher-dimensional campaign must produce a
new quantitative visibility or certification result, not a smoke ensemble.

## Selection protocol

Run one discriminating pilot per candidate.  Score novelty, theorem depth,
physical relevance, achievable scale, independent verification, and overlap
with existing Eriksson papers.  Select only after literature collision checks
and after at least one candidate has been falsified or downgraded.

## Initial literature collision notes

- Matrix moments/Ritz and finite-spectrum recovery are already established;
  they are baselines, not claims.
- TFIM form factors and DMRG excitation extraction are established; TFIM is a
  calibration model unless a new visibility-length theorem emerges.
- SU(N) glueball spectra and continuum extrapolations are established; a new
  gauge paper needs a certification/visibility result rather than another mass
  fit.
- Transition density matrices occur in quantum chemistry and integrability,
  but the specific optimal-local-probe visibility programme requires a deeper
  collision search before any novelty claim.

