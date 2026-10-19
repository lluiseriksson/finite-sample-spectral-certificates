# Quantum reservoir boundary programme

## Working title

*Architecture-Dependent Decoherence Suppression in Passive Quantum Networks:
Irreducible Channel Mixing, Squared Rate Gaps, and the Price in Dwell Time*

## Target contribution

The paper will not claim that quantum noise filtering, passive linear quantum
systems, reservoir engineering, or filter functions are new.  Its target is the
following conjunction:

1. an explicit passive lossless six-port network whose three-channel signal
   block is a rational Schur function;
2. exact delayed tangential calibration at a full-spark family of boundary
   frequencies;
3. exponentially small stopband signal leakage for the irreducible network and
   a unit lower bound for every constant-channel reducing comparator of the
   same rational bidegree;
4. a squared exponential separation after transfer through a uniformly
   nondegenerate bath spectral density;
5. an explicit delay-cost law showing that the construction pays exponential
   peak group delay for exponential suppression.

The fifth item is essential.  Without it the paper would hide the resource that
has been moved into high-Q passive memory and would overstate the operational
consequence for the Heisenberg cut.

## Construction

Let

`b_r(z)=(z-r)/(1-rz)` and `B_S(z)=eta_S b_r(z)^S`,

with `eta_S` chosen so that `B_S(-1)=-1` and
`r_S=1-exp(-alpha S)`.  For three phases
`phi_g in {-1,0,1}(1-r_S)`, define

`p_g(z)=(1+exp(i phi_g) B_S(z))/2`.

Each `p_g` is a scalar Schur function and the upper-left signal entry of the
lossless Mach-Zehnder block obtained by conjugating
`diag(1,exp(i phi_g)B_S)` with a balanced beamsplitter.

Let `F` be the normalized three-point DFT and

`U(z)=F diag(1,z,z^2) F*`,
`U_rev(z)=F diag(z^2,z,1) F*`.

Both are paraunitary polynomials and
`U_rev(z)=z^2 U(z)*` on the unit circle.  The signal transfer is

`G_S(z)=U(z) diag(p_0(z),p_1(z),p_2(z)) U_rev(z)`.

It is therefore the signal block of an explicit passive six-port inner
transfer.  At every boundary point where
`exp(i phi_g)B_S(z)=1`, the vector `v=U(z)e_g` obeys
`G_S(z)v=z^2v`.  After the invertible DFT, these vectors are Vandermonde
columns with parameters in three disjoint arcs, hence every three are linearly
independent.

## Comparator theorem skeleton

Let `Q=A/d` have a common scalar denominator of degree at most `S`, numerator
degree at most `S+4`, and a constant reducing channel line `Cu`.  Projection of
the delayed calibration onto `u` gives

`(q(z_j)-z_j^2) <u,v_j> = 0`.

Full spark removes at most two nodes.  The remaining nodes are roots of the
numerator of `q-z^2`, whose degree is at most `S+4`.  Since the construction
provides at least `3(S-1)` calibrations, the scalar branch is identically
`z^2` for `S>=5`.  Thus `||Q(z)||>=1` throughout the stopband.  The passive
irreducible signal block has norm `O(S exp(-alpha S))` there.

## Davies/resource interface

For a stationary Gaussian bath with
`j_- I <= J(omega) <= j_+ I` on the protected band, the filtered Kossakowski
matrix is

`Gamma_G(omega)=G(omega)J(omega)G(omega)*`.

The constructed network obeys

`||Gamma_G|| <= j_+ O(S^2 exp(-2 alpha S))`,

whereas any reducing comparator retains a channel with rate at least `j_-`.
The claim will first be stated for an exact multichannel pure-dephasing model;
extension to general Davies generators will be isolated as a conditional
corollary unless all domain and secular assumptions are proved.

## Exact finite no-go inherited from the previous witness

For the four degree-two calibrations of the previous paper, every symmetric
quadratic solution is `P_tau(x)=I+tau R(x)`, with

`R(0)` positive definite and `det R(5/2)=-2/3`.

If `tau>0`, contractivity fails at zero; if `tau<0`, it fails along the negative
eigendirection at `x=5/2`.  Thus global Hermitian contractivity forces
`tau=0`, the identity.  This explains why physical passivity requires the
causal complex boundary construction above rather than a rescaling of the old
zero-phase certificate.

## Release gates

- analytic proof of pass-node count and stopband bound;
- exact proof of the Vandermonde full-spark property;
- explicit six-port inner completion;
- rational root-count comparator proof;
- exact pure-dephasing/Davies interface with assumptions visible;
- numerical scaling replay and dense-grid unitarity checks;
- adversarial priority matrix;
- clean compiled and visually inspected PDF;
- independent evaluation must exceed the preceding MIMO paper on the exact PDF
  hash before the objective can close.
