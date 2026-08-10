#!/usr/bin/env python3
"""Run and hash-check the passive quantum reservoir paper release."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "programme" / "QUANTUM_RESERVOIR_ARTIFACT.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    if path.suffix in {".json", ".py"}:
        # Git checkouts and Python text writers use platform-native newlines.
        # Artifact identity for source/JSON is therefore defined after the
        # conventional LF normalization; PDFs remain byte-exact below.
        content = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        digest.update(content)
    else:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def run(*arguments: str) -> None:
    subprocess.run([sys.executable, *arguments], cwd=ROOT, check=True)


def check_hashes(expected: dict[str, str]) -> None:
    for relative, digest in expected.items():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != digest:
            raise RuntimeError(f"hash mismatch for {relative}: {observed} != {digest}")
        print(f"PASS hash {relative}: {observed}")


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected = manifest["artifacts"]
    check_hashes(expected)
    run("verification/verify_passive_quantum_filter.py")
    run("research/quantum_reservoir_filter.py")
    run("research/quantum_reservoir_robustness.py")
    check_hashes(expected)
    pdf = ROOT / "paper_quantum_reservoir" / "architecture_dependent_decoherence.pdf"
    if not pdf.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("release artifact is not a PDF")
    print("PASS quantum reservoir release gate")


if __name__ == "__main__":
    main()
