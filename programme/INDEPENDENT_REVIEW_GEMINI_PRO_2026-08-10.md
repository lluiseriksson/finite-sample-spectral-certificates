# Independent Gemini 3.1 Pro review — 10 August 2026

## Artifact and method

- Scientific commit: `46ed63d12d1639b949625f67b1dfbc175ad4580f`
- PDF: `paper_noncommuting/noncommuting_filters_draft.pdf`
- Visible title reported by the reviewer: *Exponential Calibration Laws for
  Noncommuting Matrix-Polynomial and Linear-Phase MIMO FIR Filters*
- Visible page count reported by the reviewer: 15
- Local and independently downloaded public SHA-256:
  `d62ed80dba4a4c9cd42016aa23931277154a6bdfc070686ff19b91e8f3e7e323`
- Review surface: a new Gemini 3.1 Pro conversation, with the exact repository
  PDF attached after the model correctly refused to score the raw URL when its
  URL retriever could not parse the binary.

This is an independent model audit, not named-human peer review.  The reviewer
did not compute the SHA-256 and did not run repository scripts.  Those checks
were performed separately locally, from the public bytes, and in GitHub CI.

## Blind verdict

- Scientific contribution: **9.0 / 10**.
- Manuscript quality: **8.5 / 10**.
- Recommendation: **submit for human peer review**.
- Exact scientific carrier identified: Theorem 2.1, the robust exponential
  tangential separation.

The reviewer marked all eight requested mathematical targets `PASS`:

1. the full-spark/root-count obstruction at `M=d+N`;
2. the factorial Lagrange identity and exponential Chebyshev construction,
   with constants described as sound but loose;
3. the degree-three determinant argument and quantitative singular floor;
4. the symmetric right inverse and Neumann correction;
5. the robust commuting lower bound, including the overlap and Lagrange
   factors;
6. the scalar binomial-tail/Hoeffding obstruction to constant tolerance;
7. the exact Bernstein continuum certificate and `delta<7/1920` threshold;
8. the polynomial-to-palindromic FIR conversion and restrained empirical
   interpretation.

The model reported that its counterexample attempt—forcing a commuting
solution by asymmetric channel weights—still encountered the Lemma 2.2
full-spark root obstruction.

## Reservations retained

The report kept two substantive limitations rather than converting the score
into an unconditional novelty certificate:

1. The theorem is a frequency-domain/operator-norm separation and does not
   prove a runtime or deployed Krylov-recurrence advantage.
2. The exposition moves abruptly between approximation theory, MIMO FIR
   realization, and pump-vibration data.  This density was the main reason the
   manuscript score remained below nine.

The repository priority matrix remains authoritative about literature
uncertainty.  In particular, restricted full texts are not converted into
positive evidence of priority, and the model score does not establish novelty.

## Adversarial-control trace

Before the blind Pro run, a separate Flash Advanced conversation returned a
`5.8/10` report based on three objections that did not match the PDF: it replaced
Lemma 2.2's scalar root proof by parameter counting, replaced Lemma 2.4's
three-by-three determinant by an `O(N)` Cramer matrix, and omitted the explicit
factorial cancellation in Lemma 2.3.  When challenged to quote the displayed
steps, that model retracted all three objections, plus its float-arithmetic and
scale-class objections, and rescored the manuscript `8.6/10` scientific and
`8.9/10` editorial.  Neither of those conditional scores is used as the release
gate; the independent attached-PDF Gemini 3.1 Pro run above is the recorded
blind verdict.

This trace is retained because a favorable correction should not erase an
initial model hallucination.  The exact-arithmetic verifiers and human peer
review remain more probative than any model score.
