#!/usr/bin/env python3
"""Hash, schema, and independent-verification gate for the unified paper."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "programme" / "UNIFIED_ROUTING_TABLE_MEMORY_ARTIFACT.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    if path.suffix in {".bib", ".json", ".md", ".py", ".tex", ".yaml", ".yml"}:
        digest.update(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n"))
    else:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def check_hashes(artifacts: dict[str, str]) -> None:
    for relative, expected in artifacts.items():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"hash mismatch for {relative}: {observed} != {expected}")
        print(f"PASS hash {relative}: {observed}")


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest["schema"] != "unified-routing-table-memory-artifact-v1":
        raise RuntimeError("manifest schema mismatch")
    check_hashes(manifest["artifacts"])
    subprocess.run(
        [sys.executable, "verification/verify_unified_routing_table_memory.py"],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [sys.executable, "verification/verify_detector_rank_hierarchy.py"],
        cwd=ROOT,
        check=True,
    )
    canonical = ROOT / "paper_unified_routing_table_memory/exact_memory_finite_spectral_routing_tables.pdf"
    submission = ROOT / "output/pdf/exact_memory_finite_spectral_routing_tables.pdf"
    if not canonical.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("canonical paper is not a PDF")
    if canonical.read_bytes() != submission.read_bytes():
        raise RuntimeError("submission PDF is not byte-identical to the canonical paper")
    check_hashes(manifest["artifacts"])
    print("PASS unified routing-table release gate")


if __name__ == "__main__":
    main()
