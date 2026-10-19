# Routing-word memory research record

## Result

For distinct boundary nodes, one rank-`k` input plane, and mutually orthogonal rank-`k` targets, a routing word of length `L` with occupancies `n_a` has exact minimum McMillan degree

`d_min = k (L - min_a n_a)`.

The lower bound uses the total zero budget of compressed exterior-power minors. The upper bound uses node polynomials, scalar Fejér–Riesz factorization, degree-preserving rectangular-to-square lossless completion, tensor amplification, and constant unitary alignment.

## Strong separation

At `k=1`, `AAAB` has two cyclic switches and costs three states; `ABAB` has four cyclic switches and costs two. Switch count and target span therefore do not order passive memory.

## Audit boundary

The theorem does not cover repeated nodes, nonorthogonal targets, calibrated output frames, heterogeneous ranks, active devices, or approximate routing. A directed literature search located the component tools but not the closed rarest-occupancy law. Novelty claims are restricted to this exact theorem and compiler.

## Evidence

- 1,024 seeded random words.
- All 1,155 word classes through length seven modulo alphabet renaming.
- 256 permutation checks.
- Curated and clustered-node fixtures.
- Independent verifier and release-wide hash gate.
