#!/usr/bin/env python3
"""Replay and hash-check the generic fan-out memory release."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "programme" / "GENERIC_FANOUT_MEMORY_ARTIFACT.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    if path.suffix in {".bib", ".json", ".md", ".py", ".tex", ".txt", ".yaml", ".yml"}:
        digest.update(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n"))
    else:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def compare(expected: object, observed: object, location: str = "root") -> None:
    if isinstance(expected, bool) or isinstance(expected, str) or expected is None:
        if observed != expected:
            raise RuntimeError(f"replay mismatch at {location}: {observed!r} != {expected!r}")
    elif isinstance(expected, (int, float)):
        if not isinstance(observed, (int, float)):
            raise RuntimeError(f"numeric replay type mismatch at {location}")
        if not math.isclose(float(observed), float(expected), rel_tol=5.0e-8, abs_tol=3.0e-9):
            raise RuntimeError(f"numeric replay mismatch at {location}: {observed!r} != {expected!r}")
    elif isinstance(expected, list):
        if not isinstance(observed, list) or len(observed) != len(expected):
            raise RuntimeError(f"list mismatch at {location}")
        for index, (left, right) in enumerate(zip(expected, observed, strict=True)):
            compare(left, right, f"{location}[{index}]")
    elif isinstance(expected, dict):
        if not isinstance(observed, dict) or observed.keys() != expected.keys():
            raise RuntimeError(f"mapping mismatch at {location}")
        for key in expected:
            compare(expected[key], observed[key], f"{location}.{key}")
    else:
        raise TypeError(f"unsupported type at {location}: {type(expected)}")


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
    check_hashes(manifest["artifacts"])
    subprocess.run(
        [sys.executable, "verification/verify_generic_fanout_memory.py"],
        cwd=ROOT,
        check=True,
    )
    with tempfile.TemporaryDirectory(prefix="generic-fanout-replay-") as directory:
        temporary = Path(directory)
        output = temporary / "certificate.json"
        figure = temporary / "generic_fanout_memory.pdf"
        subprocess.run(
            [
                sys.executable,
                "research/generic_fanout_memory_certificate.py",
                "--output", str(output),
                "--figure", str(figure),
            ],
            cwd=ROOT,
            check=True,
            stdout=subprocess.DEVNULL,
        )
        frozen = json.loads(
            (ROOT / "results/generic_fanout_memory/certificate.json").read_text(encoding="utf-8")
        )
        replayed = json.loads(output.read_text(encoding="utf-8"))
        compare(frozen, replayed)
        if not figure.read_bytes().startswith(b"%PDF-"):
            raise RuntimeError("replayed figure is not a PDF")
        print("PASS deterministic numerical replay")
    pdf = ROOT / "paper_generic_fanout_memory/generic_fanout_memory.pdf"
    if not pdf.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("paper artifact is not a PDF")
    check_hashes(manifest["artifacts"])
    print("PASS generic fan-out memory release gate")


if __name__ == "__main__":
    main()
