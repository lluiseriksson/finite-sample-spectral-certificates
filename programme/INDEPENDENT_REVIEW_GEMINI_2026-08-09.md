# Independent adversarial review: Gemini Advanced

Date: 2026-08-09

Reviewer surface: Google Gemini, Advanced reasoning mode, independent model
session. The reviewer received only the public 13-page PDF and the complete
text of EVALUATOR_PACKET.md; it was explicitly instructed not to reward
repository size and to reduce, rather than inflate, its scores when uncertain.

## Artifact submitted

- Scientific source commit:
  31bee2e8bc7d5b920b0b129940e63e30c4fae714
- PDF SHA-256 verified locally before upload:
  0fadc8a1d453f76897aa4375844d68ac1b8e058d9092342893fb52e307d57d12
- Title observed independently by the reviewer: *Exponential Tangential
  Advantages of Noncommuting Matrix-Polynomial and Linear-Phase MIMO FIR
  Filters*
- Page count observed independently by the reviewer: 13

Gemini stated that its attachment interface could not compute SHA-256 from raw
bytes. It therefore did **not** claim an independent hash check. The upload
path and local hash were checked by the submitting workflow; the reviewer
independently matched title, page count, date and scientific commit.

## First response block

- Scientific contribution score: **7 / 10**
- Manuscript quality score: **8 / 10**
- Headline theorem checked: Theorem 2.1, supported by Lemmas 2.2 and 2.5.
- Closest collision: Stefanovski--Georgijević (2016) for low regional norm
  under bitangential interpolation, and Ljungars--Fu (1998) for
  operator-norm multi-channel linear-phase FIR design.
- Proof defect or counterexample: **none found**.
- Recommendation: **Submit**, with minor revisions.

The reviewer identified the scientific score as coming from the fixed
three-channel exponential separation, not from the numerical volume of the
repository. It checked the common-eigenbasis root count, the full-spark
threshold M >= d+N, the Lagrange right-inverse construction and perturbation
domination, equality of constraints for the two compared classes, failure of a
commuting counterexample, the quadratic rational witness, and the
cosine-to-FIR transformation. It judged the full-spark condition generic and
the FDD statement meaningful at the preprocessing level without licensing a
deployment claim. It stated that Theorem 2.1 and Corollary 4.1 remain
publishable after deleting all numerics.

## Integrity challenge and corrected assessment

The reviewer was challenged on three overstatements in its first response:

1. its text extraction reversed the displayed sign of the commutator;
2. it described the right-inverse growth too loosely as exponential; and
3. it initially used “verified” without distinguishing PDF algebra from
   repository execution.

Its erratum states:

- reversing the commutator sign corresponds to reversing multiplication order
  and does not affect the required non-vanishing conclusion;
- the correct envelope is R_N = exp(O(N log N)), which is still dominated by
  exp(-N^3);
- the root count, Chebyshev/Lagrange estimates, matrix witness and FIR algebra
  were checked symbolically from the PDF;
- repository scripts, raw VBL binaries, hashes and solver outputs were **not**
  executed by Gemini and are not part of its independent verification; and
- restricted full texts remain a bibliographical uncertainty.

After those corrections the reviewer returned:

- Scientific contribution score: **7 / 10, unchanged**
- Manuscript quality score: **7 / 10, reduced from 8**
- Recommendation: **Minor revision**

The requested minor revision is to state R_N = exp(O(N log N)) explicitly in
Section 2.3. The second condition is the already disclosed proof-level access
audit of the closest restricted papers.

## Evidence limits

This is an independent model review, not peer review by a named human and not a
formal proof certificate. The scores apply only to the hash above. The
reviewer's script-name transcription contained one non-existent filename and
is not used as evidence of execution. Local exact verifiers and CI remain
separate evidence. Any PDF edit after this report requires a new artifact hash
and a fresh reviewer confirmation.
