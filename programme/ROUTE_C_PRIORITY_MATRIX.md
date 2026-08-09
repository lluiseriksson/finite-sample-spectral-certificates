# Route C priority matrix

Status: bounded Gate G1 audit with documented external access blocker,
2026-08-09.  “Not found” below is not a proof of novelty.  The candidate claim
remains provisional until the restricted-access rows can be checked from full
text or an author preprint.  Independent adversarial review of the exact
revised PDF returned 7/10 scientific and 7/10 manuscript, with no proof defect
identified; that assessment does not eliminate this priority uncertainty.

## Candidate theorem delta

At ordinary polynomial degree `N` and fixed channel dimension `d=3`,
full-spark tangential calibrations can force every pairwise-commuting real
symmetric coefficient polynomial to equal the identity, while an explicitly
constructed noncommuting symmetric polynomial has real-stopband operator norm
`O(exp(-cN))`.  After an affine cosine substitution, the same statement is a
fixed-latency separation between reciprocal linear-phase MIMO FIR filters and
a fixed orthogonal bank of scalar linear-phase FIR filters.

## Clause-by-clause comparison

| Source | Interpolant/domain | Complexity notion | Regional norm conclusion | Symmetry / comparator | Collision decision | Access actually checked |
|---|---|---|---|---|---|---|
| Ball--Kang (1990) | Matrix polynomials; tangential Lagrange--Sylvester data | Low McMillan degree | No regional minimax bound in the theorem-equivalent dissertation chapter | No Hermitian commuting lower bound | Kills novelty of low-degree tangential polynomial construction, not the stated separation | Kang's full 1990 dissertation, Ch. III, read; it states its main results also appear in Ball--Kang.  The journal pagination itself remains inaccessible |
| Fuhrmann (2010) | Polynomial and rational matrix tangential interpolation | Algebraic degree/realization structure | No regional minimax theorem in the abstract/portal record | No commuting comparison located | Kills novelty of the interpolation formalism | Institutional record read; OpenAlex exposes a publisher PDF URL, but ScienceDirect returns a Cloudflare error and no repository manuscript was found |
| Alpay--Lewkowicz (2014), Props. 3.1, 6.1--6.2 and Sec. 9(v) | Matrix-polynomial interpolation with Hermitian, positive-semidefinite and related symmetries | Minimal polynomial degree, sometimes reduced through structured vanishing corrections | Positivity/Hermiticity on the imaginary axis, not stopband norm minimization | Structured full matrix-value data; tangential structured interpolation is explicitly left as future work | Kills novelty of structured/Hermitian matrix-polynomial interpolation, but not the tangential commuting-class separation | Complete 20-page arXiv text inspected, including construction, minimal-degree propositions and future-work clause |
| Blomqvist--Lindquist--Nagamune (2003) | Matrix Nevanlinna--Pick interpolation on the disk | Rational complexity constraint | Schur/positive-real norm control | No fixed ordinary degree or commuting lower bound located | Kills broad “norm + interpolation + complexity is new” language | Abstract and bibliographic record |
| Kuroiwa--Lindquist (2006), Thms. 1 and 3, Sec. 3.4 | Bitangential positive-real analytic interpolation on the disk | McMillan degree at most `2(n_l+n_r)`; lower degree through restricted `Psi` | Feasibility via a Pick/positive sequence condition and a unique convex optimizer for positive `Psi` | No Hermitian real-interval polynomial or commuting comparison | Closest degree-constrained analytic interpolation result; does not imply the fixed-degree separation found here | Full conference PDF read |
| Stefanovski--Georgijević (2016), Problem 1 and Thm. 1 | Real stable rational matrix with bitangential constraints | Rational order may grow | For every `gamma>0`, norm below `gamma` on a closed proper frequency region | No Hermitian coefficient restriction or commuting comparator | Fatal to claiming that arbitrarily small regional norm under tangential constraints is new; surviving delta is fixed order/latency plus the comparator lower bound | Primary publisher HTML sections read; publisher PDF blocked by CAPTCHA, so proof-level audit remains open |
| Toh--Trefethen (1998) | Scalar polynomial `p(A)` | Ordinary polynomial degree | Matrix-operator minimax/Chebyshev problem | Same scalar polynomial, not matrix taps | Kills novelty of SDP/minimax polynomial filtering | Full theorem metadata and established paper scope |
| Gama--Marques--Ribeiro--Leus (2018), Eqs. (5)--(7) | MIMO graph filter `Y=sum S^k X H_k^T` | Number of graph taps | Learnt/parametric filtering, not the present interval minimax theorem | General MIMO taps and parsimonious variants; no commuting obstruction | Establishes the natural multichannel polynomial-filter class, not the separation | Full arXiv PDF read |
| Wu--Boyd--Vandenberghe (1998) | Scalar FIR frequency response | FIR order | Convex/minimax spectral-mask design | Scalar filters; no tangential channel signatures | Kills novelty of convex FIR design and scalar stopband optimization | Author-hosted bibliographic page and standard chapter scope |
| Kootsookos (1991), Problems 1.1--1.2, Lemmas 3.4 and 6.1, Chs. 4--6 | Matrix transfer functions approximated by `q`-coefficient MIMO FIR matrix polynomials | Fixed FIR length `q`; rational target McMillan degree | Uniform-frequency/operator `H-infinity` approximation, global lower bounds and constructive algorithms | Covers linear phase and MIMO coefficients; no tangential constraints, per-coefficient Hermiticity or commuting comparator | Fatal to claiming fixed-length MIMO uniform-norm approximation, its lower bounds or linear-phase algorithms as new; does not imply the class separation | Complete 119-page ANU doctoral thesis downloaded from the primary repository and searched/read through the problem statements, main bounds, algorithms and conclusion |
| Ljungars--Fu (1998), Sec. 4.1, Eqs. (4.5)--(4.6) | Matrix-valued multi-channel real linear-phase FIR response | Per-entry tap counts (examples use unequal diagonal/off-diagonal lengths) | Minimizes the maximum spectral/operator norm of the matrix error over a frequency grid by SDP | Fully coupled matrix response; no reciprocal-symmetry restriction, tangential pass data, or commuting comparator | Fatal to claiming that operator-norm minimax design of multi-channel linear-phase FIR filters is new; does not imply the class-separation theorem | Complete author-hosted 20-page proceedings paper inspected, including SDP formulation, experiments, limitations and conclusion |
| Hurley--Hurley (2012) | Paraunitary polynomial matrices for filter banks | FIR/matrix-polynomial constructions | Exact unitary response, not regional attenuation with directional pass data | Explicit symmetric idempotent constructions and noncommuting products | Kills any claim that noncommuting FIR matrix taps or rational examples are new; it does not supply the minimax separation | Full arXiv text read |
| Gallivan--Vandendorpe--Van Dooren (2004), Problem 1.1 and Thm. 4.2 | MIMO reduced transfer functions with left, right and two-sided tangential interpolation | Minimal McMillan degree, certified by Loewner-matrix rank | No regional or minimax objective | General strictly proper rational MIMO systems; no Hermitian or commuting comparator | Establishes generic uniqueness and minimal McMillan degree of the rational interpolant, not the stated fixed-ordinary-degree separation | Author-deposited full text inspected through the complete theorem and minimality proof |
| Gosea--G\"uttel (2021), Secs. 3--4 | Matrix-valued rational approximation and block-AAA | Rational support count / barycentric representation | Discrete-set approximation accuracy, including noisy-data comparisons | Matrix weights and tangential interpolation noted; no symmetric-polynomial or commuting lower bound | Kills novelty of matrix-weighted rational approximation and noisy block interpolation | CC-BY SIAM full text, abstract and relevant tangential-interpolation section inspected |
| Jonas--Bamieh (2026), v2, Thms. 3.1 and 4.1--4.2, Sec. V-A | Iterative frequency-domain MIMO tangential interpolation | Low-rank interpolation data and reduced-system state dimension | Uses spectral-norm maximum error to choose the next frequency; proves monotonicity only for a weighted `H2` objective and leaves monotone `H-infinity` bounds as future work | General rational reduced models; no Hermitian coefficients or simultaneously diagonalizable lower bound | Closest modern algorithmic link to maximum-error tangential MIMO approximation; theorem-level audit finds no collision | Full arXiv v2 text, theorem statements, algorithms and conclusion inspected |

## Restricted-access audit log

The remaining uncertainty is access-specific, not positive evidence of
priority.  On 2026-08-09 the following independent routes were checked:

- For Fuhrmann (2010), OpenAlex and Semantic Scholar both label the article
  bronze/open-archive but expose only the same ScienceDirect PDF URL.
  Elsevier's metadata endpoint confirms `openaccessArticle=true`, while its
  full-view text-mining endpoint returns `401 AUTHENTICATION_ERROR`; the direct
  PDF returns a robot challenge/HTTP 403 both locally and from a fresh Google
  Compute Engine egress.  The Ben-Gurion institutional record contains
  metadata but no deposited file, and ResearchGate explicitly reports no full
  text.  No CAPTCHA or authentication control was bypassed.
- For Stefanovski--Georgijević (2016), the publisher HTML exposes the abstract,
  Problem 1 setup and the introduction's statement that Theorem 1 solves it
  for every positive `gamma`.  OpenAlex and Semantic Scholar both report the
  work closed with no repository PDF.  Elsevier's full-view endpoint again
  requires authentication, and the publisher PDF is robot-challenged.

Thus the audit can already rule out broad claims of first tangential polynomial
interpolation or first arbitrarily small regional norm.  It cannot yet certify
that neither proof contains an implicit fixed-ordinary-degree,
Hermitian-coefficient or commuting-comparator specialization.

## What survives the audit so far

The following conjunction has not been found in the checked sources:

1. fixed **ordinary** degree, equivalently fixed FIR latency and tap count;
2. real symmetric/Hermitian matrix coefficients and exact directional pass
   calibrations;
3. a full pairwise-commuting comparator, equivalently one fixed orthogonal bank
   of scalar filters;
4. an explicit exponential upper/lower separation with only three channels;
   and
5. an exact rational five-tap certificate whose finite gap survives commuting
   calibration residuals `delta < 3/32300`.

Each ingredient separately is close to classical work.  Priority, if it
survives, belongs only to their theorem-level conjunction.  Gate G1's priority
clause therefore remains externally blocked by the inaccessible Fuhrmann and
Stefanovski--Georgijević proof texts.  The independent adversarial review gate
has passed at 7/10 scientific and 7/10 manuscript.  A targeted search of
reciprocal/passive MIMO FIR approximation and simultaneous diagonalization found
no theorem with the surviving conjunction, but this is not positive evidence of
priority.  The broad fixed-length MIMO `H-infinity` and linear-phase
operator-norm branches are now represented by the full Kootsookos and
Ljungars--Fu audits.

## Application audit

The graph-filter pilot is deliberately negative.  On an 18-by-23 grid graph,
the symmetric MIMO optimum collapses numerically to the identity, while a
scalar Chebyshev filter has stopband maximum `1/17` and only `0.00524` best
rescaled pass error after two iterations.  Exact directional preservation is
therefore not scientifically mandated by ordinary graph denoising.

The FIR corollary is structurally natural because frequency-dependent channel
calibrations, reciprocity, linear phase, fixed delay and filter order are
native specifications.  More specifically, frequency-domain decomposition
(FDD) estimates modal shapes from the leading singular vectors of cross-spectral
matrices.  If a modal contribution is `lambda v v^T`, the constraint `P(x)v=v`
preserves that component exactly under the palindromic FIR response, up to its
common delay.  The measured-data gate is passed by a frozen
VBL-VA001 triaxial calibration.  Five leading cospectral directions extracted
from six training records give an exact rational `<24/25` versus `1` continuum
certificate; six held-out records give maximum directional residual
`0.0038491`.  An unconstrained nonsymmetric MIMO baseline reaches much lower
training leakage but degrades to `0.10447` held-out residual, while scalar and
64-basis approximate-commuting controls are explicitly reported.

The full held-out cospectra are also frozen.  Recomputing their FDD leading
eigenpairs after filtering gives at most `0.2217` degrees of angular drift and
`2.40e-5` relative change in the leading spectral ordinate.  This closes the
constructed-data and task-meaning gates at the modal-preprocessing level, not
the deployment gate.  No classification, damping-estimation, hardware, latency
or runtime claim is licensed by this application.
