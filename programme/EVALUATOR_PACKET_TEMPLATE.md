# Evaluator packet (mandatory cover sheet)

An evaluation is valid only if the evaluator fills every field below from the
actual bytes inspected.

## Artifact identity

- Repository URL:
- Release/tag:
- Commit:
- PDF path in repository:
- Downloaded filename:
- Title visible on page 1:
- Page count:
- SHA-256 of downloaded bytes:
- Embedded manuscript version:

If any field disagrees with `ARTIFACT_MANIFEST.md`, stop: the score belongs to a
different artifact.

## Blind scoring request

Score scientific contribution and manuscript quality separately on a 0--10
scale.  Do not award points for repository size, number of experiments, or
formalization of auxiliary facts.  Evaluate:

1. theorem novelty against the three closest primary references;
2. mathematical difficulty of the decisive step;
3. whether the conclusion is thermodynamic or only finite-volume;
4. whether measurement/calibration assumptions are physically justified;
5. whether counterexamples and negative results are disclosed;
6. whether headline numerical claims have exact/interval replay;
7. whether the paper distinguishes upper-gap exclusion from positive-gap proof;
8. whether every headline claim links to a reproducible artifact.

## Adversarial questions

- Can the main theorem be obtained by directly composing a known confidence
  sequence with a known SDP hierarchy?
- Does the proof of visibility or stability assume the gap it later claims to
  infer?
- Is there a fixed finite hierarchy level with a genuine robust rejection
  margin, or only floating-point infeasibility?
- Does the claimed thermodynamic statement depend on a hidden boundary
  condition?
- Which result would remain publishable if every numerical figure were
  removed?

The evaluator must quote the exact theorem number supporting the score and list
the strongest reason the manuscript still falls below 7, if applicable.
