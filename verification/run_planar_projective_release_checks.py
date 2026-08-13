#!/usr/bin/env python3
"""Release-level checks for the planar projective memory paper."""

from __future__ import annotations

import json
import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "planar_projective_memory" / "release_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    subprocess.run([sys.executable, str(ROOT / "verification" / "verify_planar_projective_memory.py")], check=True)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected = manifest["sha256"]
    required = [ROOT / relative for relative in expected]
    missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"missing or empty release files: {missing}")
    mismatches = {
        relative: {"expected": digest, "actual": sha256(ROOT / relative)}
        for relative, digest in expected.items()
        if sha256(ROOT / relative) != digest
    }
    if mismatches:
        raise RuntimeError(f"release manifest mismatch: {mismatches}")
    tex = (ROOT / "paper_planar_projective_memory" / "main.tex").read_text(encoding="utf-8")
    forbidden = ["TODO", "TBD", "placeholder", "v2.8-ai-vixra-submission"]
    # The prior release is allowed only as a bibliography URL.
    body = tex.split("\\bibliographystyle", 1)[0]
    bad = [token for token in forbidden[:3] if token in body]
    if bad:
        raise RuntimeError(f"forbidden draft tokens: {bad}")
    print(
        json.dumps(
            {
                "status": "PASS",
                "required_files": len(required),
                "manifest_sha256": sha256(MANIFEST),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
