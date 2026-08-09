# Literature positioning (live collision audit)

This is a working novelty ledger, not a claim of priority.

## Closest recent papers

1. Xu et al., *The bulk spectral gap is semi-decidable: a convergent family of
   certified upper bounds*, arXiv:2606.03836.  This is now the essential
   comparator for any claim about certified gap upper bounds.  It starts from a
   local Hamiltonian/operator algebra and builds a convergent SDP hierarchy in
   the thermodynamic bulk.  The present project must not claim the first SDP
   upper bounds.  Its distinct target is finite-sample validity when the input is
   a noisy random Gram/correlator sketch with unknown covariance.

2. Cho, Gabai, Lin, Yeh and Zheng, *Bootstrapping Euclidean Two-point
   Correlators*, arXiv:2511.08560.  This supplies the closest correlator-bootstrap
   context: reflection positivity, equations of motion and KMS/ground-state
   positivity constrain exact continuous-time correlators.  Our differentiator
   must be statistical: a confidence region, adaptive witness selection, and an
   auditable finite-sample type-I guarantee.

3. Rao, *A convergent hierarchy of spectral gap certificates for qubit
   Hamiltonians*, arXiv:2510.08427.  This gives lower-gap SDP certificates from
   noncommutative polynomial constraints.  It reinforces an important boundary:
   moment-localizer incompatibility from a two-point function primarily rejects
   an overly large *visible* gap; it is not a general lower bound on the
   Hamiltonian mass gap.

4. Rudelson and Vershynin, *Non-asymptotic theory of random matrices: extreme
   singular values*, arXiv:1003.2990, and Vershynin, *Introduction to the
   non-asymptotic analysis of random matrices*, arXiv:1011.3027.  These are the
   source family for the rectangular-Gaussian singular-value tail that produces
   the Loewner confidence band.  The manuscript must quote a theorem with exact
   constants and hypotheses rather than cite a generic covariance-concentration
   slogan.

5. de Klerk and Laurent, *A survey of semidefinite programming approaches to
   the generalized problem of moments and their error analysis*,
   arXiv:1811.05439.  This is background for generalized moment SDPs and error
   analysis.  The affine-ellipsoid separation theorem must be framed as a
   specialized exact robust-feasibility dual, not as inventing robust moment
   optimization.

## Defensible novelty target

The strongest currently defensible package is the conjunction:

- a simultaneous, finite-sample Loewner confidence set for the *unknown* Gram
  matrix underlying Gaussian correlator sketches;
- an exact SDP compatibility test between that set and a truncated block
  Hausdorff localizer;
- an explicit semidefinite alternative whose negative margin can be replayed
  and numerically certified;
- model-specific power/volume scaling in an interacting ANNNI chain, including
  honest comparison to point GEVP/Ritz estimates and covariance-known controls.

Each ingredient has classical ancestors.  Priority, if any, can only concern
their precise combination and the demonstrated statistical/physical regime.

## Claims forbidden by the current evidence

- “First certified upper bound on a many-body spectral gap.”
- “Solution or substantial progress on the Yang--Mills mass gap.”
- “Best rank-one witness” for the plug-in minimum-eigenvector baseline.
- “Exact 95% Hotelling coverage” for the non-Gaussian quadratic features.
- “Proof of a positive mass gap” from non-rejection of a finite localizer.
- “Experimentally direct protocol” until the acquisition cost of the Gaussian
  sketch coordinates is specified for a concrete platform.
