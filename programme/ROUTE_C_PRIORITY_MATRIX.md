# Route C priority matrix

Status: open Gate G1 audit, 2026-08-09.  “Not found” below is not a proof of
novelty.  The candidate claim remains provisional until the restricted-access
rows have been checked from full text or an author preprint.

## Candidate theorem delta

At ordinary polynomial degree `N` and channel dimension `d=N^2`, fixed
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
| Blomqvist--Lindquist--Nagamune (2003) | Matrix Nevanlinna--Pick interpolation on the disk | Rational complexity constraint | Schur/positive-real norm control | No fixed ordinary degree or commuting lower bound located | Kills broad “norm + interpolation + complexity is new” language | Abstract and bibliographic record |
| Kuroiwa--Lindquist (2006), Thms. 1 and 3, Sec. 3.4 | Bitangential positive-real analytic interpolation on the disk | McMillan degree at most `2(n_l+n_r)`; lower degree through restricted `Psi` | Feasibility via a Pick/positive sequence condition and a unique convex optimizer for positive `Psi` | No Hermitian real-interval polynomial or commuting comparison | Closest degree-constrained analytic interpolation result; does not imply the fixed-degree separation found here | Full conference PDF read |
| Stefanovski--Georgijević (2016), Problem 1 and Thm. 1 | Real stable rational matrix with bitangential constraints | Rational order may grow | For every `gamma>0`, norm below `gamma` on a closed proper frequency region | No Hermitian coefficient restriction or commuting comparator | Fatal to claiming that arbitrarily small regional norm under tangential constraints is new; surviving delta is fixed order/latency plus the comparator lower bound | Primary publisher HTML sections read; publisher PDF blocked by CAPTCHA, so proof-level audit remains open |
| Toh--Trefethen (1998) | Scalar polynomial `p(A)` | Ordinary polynomial degree | Matrix-operator minimax/Chebyshev problem | Same scalar polynomial, not matrix taps | Kills novelty of SDP/minimax polynomial filtering | Full theorem metadata and established paper scope |
| Gama--Marques--Ribeiro--Leus (2018), Eqs. (5)--(7) | MIMO graph filter `Y=sum S^k X H_k^T` | Number of graph taps | Learnt/parametric filtering, not the present interval minimax theorem | General MIMO taps and parsimonious variants; no commuting obstruction | Establishes the natural multichannel polynomial-filter class, not the separation | Full arXiv PDF read |
| Wu--Boyd--Vandenberghe (1998) | Scalar FIR frequency response | FIR order | Convex/minimax spectral-mask design | Scalar filters; no tangential channel signatures | Kills novelty of convex FIR design and scalar stopband optimization | Author-hosted bibliographic page and standard chapter scope |
| Hurley--Hurley (2012) | Paraunitary polynomial matrices for filter banks | FIR/matrix-polynomial constructions | Exact unitary response, not regional attenuation with directional pass data | Explicit symmetric idempotent constructions and noncommuting products | Kills any claim that noncommuting FIR matrix taps or rational examples are new; it does not supply the minimax separation | Full arXiv text read |
| Gallivan--Vandendorpe--Van Dooren (2004) | MIMO reduced transfer functions with tangential interpolation | Minimal McMillan degree | Model-reduction interpolation, not a real-interval stopband minimax theorem | General rational MIMO systems; no Hermitian commuting comparator | Establishes the classical MIMO/Krylov meaning of tangential interpolation, not the stated fixed-degree separation | Primary SIAM abstract, metadata and reference list checked; author-shared full text located but theorem-level audit remains open |
| Gosea--G\"uttel (2021), Secs. 3--4 | Matrix-valued rational approximation and block-AAA | Rational support count / barycentric representation | Discrete-set approximation accuracy, including noisy-data comparisons | Matrix weights and tangential interpolation noted; no symmetric-polynomial or commuting lower bound | Kills novelty of matrix-weighted rational approximation and noisy block interpolation | CC-BY SIAM full text, abstract and relevant tangential-interpolation section inspected |
| Jonas--Bamieh (2026), v2 | Iterative frequency-domain MIMO tangential interpolation | Low-rank interpolation data and reduced-system order | Point selection includes maximum-error criteria; proves weighted-H2 monotonicity | General rational reduced models; no Hermitian coefficients or simultaneously diagonalizable lower bound | Closest modern algorithmic link to maximum-error tangential MIMO approximation; no theorem collision located | arXiv v2 abstract and method statements inspected; submitted to IEEE TAC |

## What survives the audit so far

The following conjunction has not been found in the checked sources:

1. fixed **ordinary** degree, equivalently fixed FIR latency and tap count;
2. real symmetric/Hermitian matrix coefficients and exact directional pass
   calibrations;
3. a full pairwise-commuting comparator, equivalently one fixed orthogonal bank
   of scalar filters;
4. an explicit exponential upper/lower separation with `d=N^2`; and
5. an exact rational five-tap certificate whose finite gap survives commuting
   calibration residuals `delta < 3/32300`.

Each ingredient separately is close to classical work.  Priority, if it
survives, belongs only to their theorem-level conjunction.  Gate G1 therefore
remains open pending the inaccessible Fuhrmann and Stefanovski--Georgijević
proof texts, plus a search targeted specifically at paraunitary,
reciprocal, and linear-phase MIMO FIR approximation.

## Application audit

The graph-filter pilot is deliberately negative.  On an 18-by-23 grid graph,
the symmetric MIMO optimum collapses numerically to the identity, while a
scalar Chebyshev filter has stopband maximum `1/17` and only `0.00524` best
rescaled pass error after two iterations.  Exact directional preservation is
therefore not scientifically mandated by ordinary graph denoising.

The FIR corollary is structurally natural because frequency-dependent channel
calibrations, reciprocity, linear phase, fixed delay and filter order are
native specifications.  The measured-data gate is now passed by a frozen
VBL-VA001 triaxial calibration.  Five leading cospectral directions extracted
from six training records give an exact rational `<24/25` versus `1` continuum
certificate; six held-out records give maximum directional residual
`0.0038491`.  An unconstrained nonsymmetric MIMO baseline reaches much lower
training leakage but degrades to `0.10447` held-out residual, while scalar and
64-basis approximate-commuting controls are explicitly reported.

This closes the constructed-data objection, not the deployment gate.  The
normal-condition spectral summary does not establish that a particular pump
diagnostic mandates exact direction preservation, reciprocity or the selected
fixed delay.  No classification, hardware, latency or runtime claim is
licensed by this application.
