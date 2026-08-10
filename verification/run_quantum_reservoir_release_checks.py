#!/usr/bin/env python3
"""Run and hash-check the passive quantum reservoir paper release."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "programme" / "QUANTUM_RESERVOIR_ARTIFACT.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    if path.suffix in {".bib", ".json", ".md", ".py", ".tex", ".txt", ".yaml", ".yml"}:
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
    subprocess.run(
        [sys.executable, *arguments],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
    )


def compare_replay(expected: object, observed: object, location: str = "root") -> None:
    if isinstance(expected, bool) or isinstance(expected, str) or expected is None:
        if observed != expected:
            raise RuntimeError(f"replay mismatch at {location}: {observed!r} != {expected!r}")
    elif isinstance(expected, (int, float)):
        if not isinstance(observed, (int, float)) or not math.isclose(
            float(observed), float(expected), rel_tol=1e-10, abs_tol=1e-12
        ):
            raise RuntimeError(f"numeric replay mismatch at {location}: {observed!r} != {expected!r}")
    elif isinstance(expected, list):
        if not isinstance(observed, list) or len(observed) != len(expected):
            raise RuntimeError(f"list replay mismatch at {location}")
        for index, (expected_item, observed_item) in enumerate(zip(expected, observed, strict=True)):
            compare_replay(expected_item, observed_item, f"{location}[{index}]")
    elif isinstance(expected, dict):
        if not isinstance(observed, dict) or observed.keys() != expected.keys():
            raise RuntimeError(f"mapping replay mismatch at {location}")
        for key in expected:
            compare_replay(expected[key], observed[key], f"{location}.{key}")
    else:
        raise TypeError(f"unsupported replay value at {location}: {type(expected)}")


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
    with tempfile.TemporaryDirectory(prefix="quantum-reservoir-replay-") as directory:
        temporary = Path(directory)
        scaling = temporary / "scaling.json"
        tolerance = temporary / "tolerance.json"
        robust = temporary / "robust.json"
        closed = temporary / "closed.json"
        run(
            "research/quantum_reservoir_filter.py",
            "--output", str(scaling),
            "--figure", str(temporary / "scaling.pdf"),
            "--architecture-figure", str(temporary / "architecture.pdf"),
        )
        run(
            "research/quantum_reservoir_robustness.py",
            "--output", str(tolerance),
            "--figure", str(temporary / "tolerance.pdf"),
        )
        run(
            "research/quantum_reservoir_robust_bound.py",
            "--output", str(robust),
        )
        run(
            "research/quantum_reservoir_closed_dephasing.py",
            "--scaling", str(scaling),
            "--output", str(closed),
            "--figure", str(temporary / "closed.pdf"),
        )
        compare_replay(
            json.loads((ROOT / "results/quantum_reservoir/passive_filter_scaling.json").read_text(encoding="utf-8")),
            json.loads(scaling.read_text(encoding="utf-8")),
            "scaling",
        )
        compare_replay(
            json.loads((ROOT / "results/quantum_reservoir/passive_tolerance.json").read_text(encoding="utf-8")),
            json.loads(tolerance.read_text(encoding="utf-8")),
            "tolerance",
        )
        compare_replay(
            json.loads((ROOT / "results/quantum_reservoir/robust_separation.json").read_text(encoding="utf-8")),
            json.loads(robust.read_text(encoding="utf-8")),
            "robust",
        )
        compare_replay(
            json.loads((ROOT / "results/quantum_reservoir/closed_dephasing.json").read_text(encoding="utf-8")),
            json.loads(closed.read_text(encoding="utf-8")),
            "closed",
        )
        print("PASS cross-platform numerical replay tolerance")
    check_hashes(expected)
    pdf = ROOT / "paper_quantum_reservoir" / "architecture_dependent_decoherence.pdf"
    if not pdf.read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("release artifact is not a PDF")
    print("PASS quantum reservoir release gate")


if __name__ == "__main__":
    main()
