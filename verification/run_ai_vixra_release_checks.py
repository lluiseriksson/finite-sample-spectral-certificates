"""Run the frozen numerical, exact-arithmetic, formal, and artifact checks.

Invoke from the repository root with

    python verification/run_ai_vixra_release_checks.py

The script deliberately captures verbose experiment output and prints one line per
gate.  On failure it exposes the complete captured transcript.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "paper_noncommuting" / "noncommuting_filters_draft.pdf"
EXPECTED_PDF_SHA256 = "ad711a211c3457ba0fef30c827b5fa543403023cc485ab16b34d0ead393c0a81"
EXPECTED_VERSIONS = {
    "numpy": "2.5.1",
    "scipy": "1.18.0",
    "cvxpy": "1.9.2",
}


def run_gate(label: str, command: list[str], marker: str, cwd: Path = ROOT) -> None:
    completed = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    transcript = completed.stdout + completed.stderr
    if completed.returncode or marker not in transcript:
        print(f"FAIL  {label}")
        print(transcript)
        raise SystemExit(completed.returncode or 1)
    print(f"PASS  {label}")


def check_environment() -> None:
    actual = {name: importlib.metadata.version(name) for name in EXPECTED_VERSIONS}
    if actual != EXPECTED_VERSIONS:
        raise SystemExit(f"environment version mismatch: expected {EXPECTED_VERSIONS}, got {actual}")
    print(
        "PASS  Python environment "
        + sys.version.split()[0]
        + " ("
        + ", ".join(f"{name}={version}" for name, version in actual.items())
        + ")"
    )


def check_lean_boundary() -> None:
    lean_files = sorted(
        path
        for path in (ROOT / "formal").rglob("*.lean")
        if ".lake" not in path.relative_to(ROOT / "formal").parts
    )
    forbidden = re.compile(r"(^|\s)(sorry|axiom)\b")
    offenders = []
    for path in lean_files:
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            code = line.split("--", 1)[0]
            if forbidden.search(code):
                offenders.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()}")
    if offenders:
        raise SystemExit("admitted or axiomatic Lean declaration found:\n" + "\n".join(offenders))
    print(f"PASS  Lean source boundary ({len(lean_files)} files; no sorry or user axiom)")


def check_pdf() -> None:
    digest = hashlib.sha256(PDF.read_bytes()).hexdigest()
    if digest != EXPECTED_PDF_SHA256:
        raise SystemExit(f"PDF SHA-256 mismatch: expected {EXPECTED_PDF_SHA256}, got {digest}")
    print(f"PASS  frozen 18-page PDF SHA-256 {digest}")


def main() -> None:
    check_environment()
    python = sys.executable
    replay_dir = ROOT / "tmp" / "ai-vixra-release-replay"
    replay_dir.mkdir(parents=True, exist_ok=True)
    gates = [
        (
            "asymptotic constant-channel arithmetic",
            [python, "verification/verify_constant_channel_scaling.py"],
            "robust constant-channel scaling arithmetic: PASS",
            ROOT,
        ),
        (
            "exact five-tap rational certificate",
            [python, "verification/verify_quadratic_filter_witness.py"],
            "strict separation survives every delta < 7/1920",
            ROOT,
        ),
        (
            "VBL-VA001 exact continuum certificate",
            [python, "verification/verify_vbl_va001_witness.py"],
            "VBL-VA001 exact witness: PASS",
            ROOT,
        ),
        (
            "deliberately negative graph-filter control",
            [python, "research/route_c_graph_filter_pilot.py", "--output", str(replay_dir / "graph_filter_pilot.json")],
            "natural graph pilot remains a negative result",
            ROOT,
        ),
        (
            "block-Krylov benchmark replay",
            [python, "research/route_c_block_krylov_benchmark.py", "--output", str(replay_dir / "block_krylov_benchmark.json")],
            '"exact_certificate_bound": 0.78125',
            ROOT,
        ),
        (
            "held-out VBL-VA001 replay",
            [python, "research/route_c_vibration_calibration.py", "--output", str(replay_dir / "vbl_va001_calibration.json")],
            '"maximum_held_out_residual": 0.00384905282359036',
            ROOT,
        ),
        (
            "Lean/mathlib kernel build",
            ["lake", "build"],
            "Build completed successfully",
            ROOT / "formal",
        ),
    ]
    for gate in gates:
        run_gate(*gate)
    check_lean_boundary()
    check_pdf()
    print("PASS  all ai.viXra release gates")


if __name__ == "__main__":
    main()
