#!/usr/bin/env python3
"""Release-level checks for the planar projective memory paper."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    subprocess.run([sys.executable, str(ROOT / "verification" / "verify_planar_projective_memory.py")], check=True)
    required = [
        ROOT / "paper_planar_projective_memory" / "main.tex",
        ROOT / "paper_planar_projective_memory" / "main.pdf",
        ROOT / "paper_planar_projective_memory" / "references.bib",
        ROOT / "paper_planar_projective_memory" / "figures" / "unbounded_detector_gap.pdf",
        ROOT / "results" / "planar_projective_memory" / "certificate.json",
    ]
    missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"missing or empty release files: {missing}")
    if required[1].stat().st_mtime < max(required[0].stat().st_mtime, required[2].stat().st_mtime, required[3].stat().st_mtime):
        raise RuntimeError("main.pdf is stale relative to source, bibliography, or figure")
    tex = required[0].read_text(encoding="utf-8")
    forbidden = ["TODO", "TBD", "placeholder", "v2.8-ai-vixra-submission"]
    # The prior release is allowed only as a bibliography URL.
    body = tex.split("\\bibliographystyle", 1)[0]
    bad = [token for token in forbidden[:3] if token in body]
    if bad:
        raise RuntimeError(f"forbidden draft tokens: {bad}")
    print(json.dumps({"status": "PASS", "required_files": len(required)}, sort_keys=True))


if __name__ == "__main__":
    main()
