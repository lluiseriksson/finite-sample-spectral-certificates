# The Rarest Target Sets Passive Memory

Reproducible source for the routing-word memory paper.

## Replay

```powershell
python research/routing_word_memory_certificate.py --output results/routing_word_memory/certificate.json --figure paper_routing_word_memory/figures/routing_word_memory.pdf --trials 1024
python verification/verify_routing_word_memory.py
```

Compile `main.tex` with `pdflatex`, `bibtex`, `pdflatex`, `pdflatex`.
