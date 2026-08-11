# Exact Memory of Finite Spectral Routing Tables

This paper proves the direct-sum occupancy law

`d_min = k * (L - min_a n_a)`

for exact passive routing of one rank-`k` input subspace among repeated,
not necessarily orthogonal, direct-sum target subspaces.

Reproduce the frozen evidence with:

```text
python research/unified_routing_table_certificate.py --output results/unified_routing_table_memory/certificate.json --figure paper_unified_routing_table_memory/figures/unified_routing_table_memory.pdf
python verification/verify_unified_routing_table_memory.py
```

The numerical artifact audits the proof construction and fails closed; it
does not replace the analytic proof.
