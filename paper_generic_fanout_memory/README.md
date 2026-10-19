# Exact memory of direct-sum spectral fan-out

This directory contains the new paper **Exact Memory of Direct-Sum Spectral
Fan-Out: Generic Maximality, Colliding Targets, and the Three-Line Phase
Diagram**.

The main theorem replaces orthogonal destinations by the open dense condition
that the target planes are in direct-sum position. The exact degree remains
`k(L-1)`, even for planes converging to coincidence.

Reproduce from the repository root:

```text
python research/generic_fanout_memory_certificate.py
python verification/verify_generic_fanout_memory.py
python verification/run_generic_fanout_memory_release_checks.py
```

Build the manuscript from this directory with `pdflatex`, `bibtex`, and two
further `pdflatex` passes. Frozen outputs and their hashes are recorded in
`programme/GENERIC_FANOUT_MEMORY_ARTIFACT.json`.
