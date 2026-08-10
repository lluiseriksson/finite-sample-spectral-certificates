# Pilot and negative-results ledger

This file records results before a manuscript claim is selected.  Numerical
values here are exploratory and may be superseded by fixed production runs.

## A. Reduced-transition visibility

### Exact pilot

- Models: TFIM paramagnetic/critical/ferromagnetic; the prior ANNNI point; a
  frustrated ANNNI point.
- Sizes: L = 6, 8, 10, 12; centered blocks through six sites.
- Ground truth: sparse diagonalization of the lowest three states.
- Quantity: `X_A = Tr_{A^c}|psi0><psi1|`; Frobenius norm for aggregate
  Hilbert--Schmidt visibility and nuclear norm for the optimal contraction.

Observed pattern:

- ordered TFIM/ANNNI: one-site visibility approaches a nonzero value while the
  even--odd finite-size gap closes exponentially;
- critical TFIM: one-site visibility decreases with the expected order-field
  finite-size power;
- paramagnetic TFIM/ANNNI: fixed-support visibility decreases with L, while a
  growing block recovers the delocalized excitation.

### Collision and decision

The operator-duality identities are elementary, TFIM form factors and critical
scaling are established, and a May-2026 paper studies low rank of reduced
transition matrices in another setting.  Candidate A is therefore downgraded
from standalone headline to a possible model-derived input/diagnostic for the
covariance-certificate paper.  It can be promoted only if the MPS campaign
finds and explains a new universal visibility-length law beyond known form
factors.

## B. Covariance-optimal PSD localizers

### Mathematical pilot

For a joint moment confidence ellipsoid, the support function gives the exact
simultaneous margin for every lifted PSD polynomial mixture `X`.  Minimizing

    tr(A(m_hat) X) + q ||C^(1/2) l(X)||_2

over `X >= 0`, `tr X = 1`, is a compact SDP/SOCP.  The test remains valid after
data-dependent optimization because the single ellipsoid controls all `X`
simultaneously.  In exact data, negative trace against a PSD `X` implies that
the localizer is not PSD.  Under noise, the aggregate `X` is itself the robust
witness; no rank-one extraction is assumed.

### Conditioning incident

The first implementation used monomial coordinates and point normalization.
At the boundary null it caused solver failures and could become unbounded in
directions vanishing at the normalization point.  This was not hidden:

1. coordinates were changed to shifted Chebyshev polynomials;
2. point normalization was replaced by the homogeneous trace-one slice;
3. every solver output is projected to PSD, renormalized, and the reported
   certificate is recomputed from that explicit feasible matrix.

### Boundary-null campaign

Parameters: claimed edge 0.8, visible alternative atom at 0.9 with weight
0.01, degree 5, AR(1) correlation 0.92, moment standard-deviation scale 1e-6,
95% chi-square ellipsoid, 1,000 null and 1,000 alternative trials.

| Method | Boundary-null detection | Alternative detection |
|---|---:|---:|
| Fixed Chebyshev + ellipsoid | 0.0% | 0.0% |
| Fixed Chebyshev + Bonferroni box | 0.0% | 0.0% |
| Optimized PSD ellipsoid | 1.1% | 100.0% |

The empirical ellipsoid miss rate was 4.5% under the null; there were zero
optimized detections while the true moment vector lay inside the ellipsoid.
Thus the observed 1.1% false-positive rate is compatible with and below the 5%
simultaneous guarantee.  The alternative run had a 5.2% ellipsoid miss rate
and 100% detection.

This passes the pilot power/coverage kill criterion.  It does not yet establish
novelty, estimated-covariance coverage, robustness to non-Gaussian MCMC data,
or a physical conclusion.

### Unknown covariance: Hotelling pilot

The same adaptive optimization was repeated with covariance estimated from
40 independent multivariate-normal replicates.  The simultaneous radius used
the exact one-sample Hotelling law

    T^2 = n (m_bar-mu)' S^{-1} (m_bar-mu),
    (n-p) T^2 / (p(n-1)) ~ F_{p,n-p}.

At degree 5, per-replicate noise scale 5e-6, AR(1) correlation 0.92, and 1,000
null plus 1,000 alternative trials:

| Method | Boundary-null detection | Alternative detection |
|---|---:|---:|
| Fixed Chebyshev + Hotelling | 0.0% | 0.0% |
| Fixed Chebyshev + Bonferroni-t | 0.0% | 0.0% |
| Optimized PSD + Hotelling | 1.7% | 100.0% |

The observed Hotelling miss rate was 6.1% under the boundary null (sampling
variation around the nominal 5%); no optimized null detection occurred while
the true vector lay inside the Hotelling ellipsoid.  This supports the exact
finite-sample argument and shows that the gain does not rely on known
covariance.  Production must enlarge the trial count, vary n/p, and include
non-Gaussian and autocorrelated batches.

### Primal--dual equivalence replay

The ellipsoidal primal maximizes the smallest attainable localizer eigenvalue
over all moment vectors in the confidence ellipsoid.  Its conic dual is the
trace-one PSD separator used above.  Across 50 boundary-null and 50 alternative
draws there were zero sign disagreements.  The largest absolute primal--dual
objective discrepancy was 1.37e-9 for the null and 8.38e-10 for the
alternative.  This numerically supports the strong-separation/minimax theorem;
it is not a substitute for its proof.

### Interacting ANNNI with non-Gaussian sketches

For L=12, two parity-resolved probes `(X,Z)`, and degree 2, exact sparse
propagation constructs the joint Gaussian process
`z_{k,a}=<g,T^(k/2) psi_a>`.  Each replicate reports the rank-one products
`z_k z_k^T`, an unbiased estimator of the block moment `B_k`.  The resulting
18-dimensional samples are Wishart-like and strongly correlated across time
and channel; they are not multivariate normal.

The first result was negative and important: with 120 samples, a nominal 95%
Hotelling ellipsoid covered the true ANNNI moment vector only 83.4% of the time
over 500 ensembles.  Its 97.4% optimized detection rate therefore was not a
valid 95% certificate.  Coverage approached the nominal level only gradually
with increasing sample count.  The manuscript must not use Hotelling outside
its Gaussian hypothesis without calibration.

For this Gaussian-sketch model the feature covariance is available exactly by
Isserlis' theorem.  Multivariate Chebyshev therefore gives the distribution-
free ellipsoid `q^2=p/alpha`.  Though conservative, the covariance-optimized
PSD separator showed:

| Samples | PSD detection | Adaptive rank-one detection | observed coverage |
|---:|---:|---:|---:|
| 120 | 0% | 0% | 100% |
| 500 | 0% | 0% | 100% |
| 1,000 | 97% | 0% | 100% |
| 2,000 | 100% | 0% | 100% |
| 5,000 | 100% | 100% | 100% |

Thus the covariance-aware SDP reduces the observed certified sample threshold
by roughly a factor of five relative to the *plug-in minimum-eigenvector*
direction in this instance.  This pilot did not optimize over all rank-one
directions; calling the comparator the best rank-one certificate would be
incorrect.  In fact, many SDP optima above threshold were numerically rank one,
showing that covariance-aware direction selection, rather than mixing alone,
drives much of the observed gain.  At the true ANNNI transfer edge (gap factor 1.0),
1,000 independent ensembles of 1,000 sketches produced zero rejections for
both methods; distribution-free coverage was 100% and nominal Hotelling
coverage was 94.7%.

These are still pilot rather than production numbers.  The known-covariance
Chebyshev ellipsoid is exact but conservative; unknown-covariance,
autocorrelated, and model-misspecified variants remain open.

### Unknown covariance: Wishart--Loewner pilot

The Gaussian-sketch coordinates themselves are jointly Gaussian even though
their quadratic moment features are not.  This permits a different exact
finite-sample route.  Extreme-singular-value concentration gives a simultaneous
Loewner confidence band `a S <= K <= b S` for the unknown joint covariance
`K`, based only on its empirical Gram matrix `S`.  Compatibility of any `K` in
that band with the proposed spectral edge is one SDP feasibility problem.

In a first L=8 smoke (10 ensembles per cell, degree 2), the true edge produced
0/10 rejections at n=500 and n=1,000.  A gap exaggerated by 35% produced 8/10
rejections at n=500 and 10/10 at n=1,000.  This is promising, not a power
estimate: production needs many more ensembles, several volumes, explicit dual
infeasibility certificates, and solver cross-checks.  Null solves are also
ill-conditioned and often report `optimal_inaccurate`; no scientific claim may
rest on that status without a residual audit.

The explicit semidefinite alternative was then implemented.  In a 10-ensemble
dual smoke at L=8 and n=500, its null values lay between `1.8e-10` and
`6.5e-10`, while the 35%-exaggerated-gap values lay between `-8.2e-4` and
`-4.1e-4`.  Maximum stationarity residual was `2.3e-9`.  A separate primal--
dual replay exposed two false `optimal_inaccurate` primal classifications among
five alternatives; the dual certificates remained negative by more than
`4e-4`.  Production therefore uses the dual margin and records residuals rather
than interpreting primal solver status as scientific evidence.

## C. Gauge Monte Carlo route

The current `ym-lattice-numerics` engine is a correct reference/smoke code for
plaquettes and rectangular Wilson loops.  It has no production glueball
operator basis, smearing, temporal correlation matrices, autocorrelation-time
analysis, or continuum-grade ensembles.  Candidate C is downgraded as a near-
term headline: presenting its smoke data would weaken the project.  It remains
a later external validation target after the statistical theorem and spin-chain
production pipeline are mature.

## Current provisional selection

Candidate B is the lead.  Candidate A supplies the concrete interacting-model
visibility study if it survives MPS scaling; Candidate C is deferred.  Before
paper drafting, the programme still requires:

- a collision search specific to adaptive robust localizer SDPs;
- a proof and dual audit of the unknown-covariance Wishart--Loewner theorem;
- baselines including data-trained rank-one filters and GEVP;
- ANNNI correlator ensembles with realistic estimated covariance;
- covariance misspecification and non-Gaussian/autocorrelated stress tests;
- interval or rational replay of selected negative certificates.

## Frozen production outcome (supersedes provisional checklist)

The selected contribution is now the unknown-covariance Wishart--Loewner
certificate with a finite-sample power theorem and exact-rational
representative replay.  The frozen campaigns contain 16,800 SDPs and 32.4
million primitive sketch vectors:

- main degree-two grid: 0/2,400 null rejections, 2,400/2,400 rejections at 35%
  gap inflation, and a resolved 15% sample-size transition through `L=16`;
- degree ablation: 0/3,600 null rejections; degree three strongly dominates at
  large volume despite its wider confidence band;
- primitive-law stress: Gaussian coverage 800/800, Student-t5 coverage
  607/800, and Rademacher coverage 800/800.  Only the Gaussian row is licensed
  by the theorem;
- exact replay: all five rational matrices pass positive LDL pivots and the
  conservative contradiction upper bound is `-6.47697914332e-4`.

The Hotelling route remains a documented negative result.  The Yang--Mills
route remains deferred: no claim in this paper converts finite ANNNI evidence
into a thermodynamic or gauge-theory mass-gap theorem.

## Route C natural-application gate: negative graph result

`research/route_c_graph_filter_pilot.py` tested a standard degree-two MIMO
graph filter on two smooth signals over an 18-by-23 grid graph.  The four
largest distinct graph frequencies supplied full-spark two-channel pass rows
(minimum absolute two-by-two minor `0.02718`); 207 nonpositive eigenvalues
formed the stopband.

The symmetric-tap SDP returned validation leakage `0.99999992` and commutators
of order `1e-14`, numerically collapsing to the identity.  A general
nonsymmetric MIMO filter reached leakage `0.60291`.  The scalar Chebyshev
baseline reached `1/17 = 0.05882`; after two iterations it had best globally
rescaled pass error `0.00524` and desired-subspace sine `0.0292`, while the
symmetric MIMO result remained at `0.997`.

Decision: **ordinary graph denoising fails the application gate**.  Exact
directional preservation is task-designed rather than scientifically compelled
here.  The numbers remain public in
`results/route_c/graph_filter_pilot.json` and the pilot runs in CI.  They must
not be reframed as positive evidence.

The surviving application interpretation is reciprocal linear-phase MIMO FIR
filtering.  The theorem's affine cosine transform fixes order, latency,
reciprocity and directional frequency calibrations exactly, and the quadratic
certificate becomes a rational five-tap filter.  The next gate is the measured
calibration reported below.

## Route C measured-calibration gate: positive VBL-VA001 result

The public CC-BY-4.0 VBL-VA001 archive (Zenodo DOI
`10.5281/zenodo.7006575`) supplies triaxial pump-vibration recordings at 20 kHz.
The acquisition manifest fixes twelve normal-condition CSV members, their
CRC-32 values and SHA-256 hashes.  Six `normal_000` records form the training
split and six `normal_001` records are held out.  After polyphase decimation to
2.5 kHz, the five largest-trace training cospectral peaks between 20 and 1000
Hz give five three-channel signatures.  No held-out record enters peak
selection or fitting.

The six-decimal rationalization is full spark in exact arithmetic (minimum
normalized 3-by-3 minor `0.009489456049`).  The exact verifier constructs a
noncommuting symmetric quadratic meeting all five rational calibrations and
proves continuum leakage below `24/25` on `[-1,-0.929]`.  It checks
`(191/200)^2 I-P(x_i)^2` by rational Sylvester tests on 101 grid points and
closes the cells with the exact derivative bound, producing the sharper global
upper bound `0.957037877842`.  Lemma 2.2 forces every exactly calibrated
commuting symmetric quadratic to be the identity, so its leakage is exactly
one.

The rational witness has measured-training residual `3.50e-7` from rounding
and maximum held-out residual `0.0038491`; train-to-held-out direction angles
are at most `0.4599` degrees.  A general nonsymmetric exact-fit SDP reaches
leakage `0.14514`, but its held-out residual is `0.10447` and maximum
coefficient Frobenius norm `36.06`, versus `2.623` for the symmetric witness.
A scalar quadratic allowed 1% calibration error reaches `0.97520`; the best of
64 seeded fixed orthogonal bases reaches `0.96966`.  The latter is explicitly a
stress test, not a global certificate over bases.

For a task-level check, the complete `3`-by-`3` cospectral matrices at all five
peaks were frozen and replayed as an FDD preprocessing problem.  The constraint
`P(x_j)v_j=v_j` exactly preserves an ideal rank-one modal component
`lambda_j v_j v_j^T`, up to the common FIR delay.  On the held-out cospectra,
the rational witness rotates the recomputed leading FDD direction by at most
`0.2217` degrees and changes its leading spectral ordinate by at most `2.40e-5`
relatively.

Decision: **the measured-data and modal-task grounding gates pass**.  The
deployment gate does not: this experiment proves neither a fault-classification
or damping-estimation improvement nor a latency, hardware or runtime advantage.
Replays are
`verification/verify_vbl_va001_witness.py`,
`research/route_c_vibration_calibration.py`, and
`results/route_c/vbl_va001_calibration.json`.
