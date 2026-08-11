# Exact Memory of Finite Spectral Routing Tables

This paper proves the direct-sum occupancy law

`d_min = k * (L - min_a n_a)`

for exact passive routing of one rank-`k` input subspace among repeated,
not necessarily orthogonal, direct-sum target subspaces.

For arbitrary equal-rank target arrangements it also proves the universal
detector-rank hierarchy

`d_min >= max_D sum_a n_a * (k - rank(D|Y_a))`.

For lines this is the maximum weighted occupancy contained in a proper
hyperplane.  It resolves the full three-line phase diagram, including an
explicit one-state `CP^1` compiler for three distinct coplanar lines.

Reproduce the frozen evidence with:

```text
python research/unified_routing_table_certificate.py --output results/unified_routing_table_memory/certificate.json --figure paper_unified_routing_table_memory/figures/unified_routing_table_memory.pdf
python research/detector_rank_hierarchy_certificate.py --output results/unified_routing_table_memory/detector_rank_certificate.json
python verification/verify_unified_routing_table_memory.py
python verification/verify_detector_rank_hierarchy.py
```

The floating-point artifact audits the direct-sum compiler.  The
detector-rank artifact uses exact rational arithmetic.  Both fail closed and
remain subordinate to the analytic proofs.
