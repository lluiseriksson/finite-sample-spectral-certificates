#!/usr/bin/env python3
"""Range-download only the twelve VBL-VA001 records used by Route C."""

import argparse
import hashlib
import json
from pathlib import Path

from remotezip import RemoteZip


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "route_c" / "vbl_va001_calibration.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "raw" / "vbl_va001")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    archive = manifest["dataset"]["archive"]
    selected = {entry["member"]: entry for entry in manifest["selection"]["files"]}
    args.output.mkdir(parents=True, exist_ok=True)

    with RemoteZip(archive["content_url"], headers={"User-Agent": "finite-sample-spectral-certificates/1.0"}) as remote:
        available = {entry.filename: entry for entry in remote.infolist()}
        for member, expected in selected.items():
            info = available[member]
            assert info.file_size == expected["size_bytes"]
            assert f"{info.CRC:08x}" == expected["crc32"]
            destination = args.output / Path(member).name
            if not destination.exists() or sha256(destination) != expected["sha256"]:
                with remote.open(member) as source, destination.open("wb") as target:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(block)
            actual = sha256(destination)
            if actual != expected["sha256"]:
                raise RuntimeError(f"SHA-256 mismatch for {destination.name}: {actual}")
            print(f"PASS {destination.name} {actual}")


if __name__ == "__main__":
    main()
