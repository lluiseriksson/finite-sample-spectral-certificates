#!/usr/bin/env python3
"""Re-extract the frozen VBL-VA001 cospectral calibration from raw CSVs."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, resample_poly


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "route_c" / "vbl_va001_calibration.json"


def cospectrum(paths, pipeline):
    nfft = pipeline["nfft"]
    hop = pipeline["hop"]
    window = np.hanning(nfft)
    result = np.zeros((nfft // 2 + 1, 3, 3), dtype=np.complex128)
    count = 0
    for path in paths:
        data = np.loadtxt(path, delimiter=",")[:, 1:4]
        data -= data.mean(axis=0)
        data = resample_poly(data, up=1, down=pipeline["decimation_factor"], axis=0)
        for start in range(0, len(data) - nfft + 1, hop):
            transform = np.fft.rfft(data[start : start + nfft] * window[:, None], axis=0)
            result += np.einsum("fi,fj->fij", transform.conj(), transform)
            count += 1
    return result / count


def leading_directions(spectral, indices):
    result = []
    for index in indices:
        _, eigenvectors = np.linalg.eigh(spectral[index].real)
        vector = eigenvectors[:, -1]
        if vector[np.argmax(np.abs(vector))] < 0:
            vector = -vector
        result.append(vector)
    return np.array(result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "raw" / "vbl_va001")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    pipeline = manifest["spectral_pipeline"]
    files = manifest["selection"]["files"]
    paths = {
        split: [args.input / Path(entry["member"]).name for entry in files if entry["split"] == split]
        for split in ("train", "held_out")
    }
    missing = [str(path) for split in paths.values() for path in split if not path.exists()]
    if missing:
        raise FileNotFoundError("Run acquire_vbl_va001_subset.py first; missing: " + ", ".join(missing))

    train = cospectrum(paths["train"], pipeline)
    held_out = cospectrum(paths["held_out"], pipeline)
    rate = pipeline["design_rate_hz"]
    frequencies = np.fft.rfftfreq(pipeline["nfft"], 1 / rate)
    power = np.trace(train.real, axis1=1, axis2=2)
    low, high = pipeline["peak_search_hz"]
    candidates = np.flatnonzero((frequencies >= low) & (frequencies <= high))
    peaks, _ = find_peaks(
        np.log(power[candidates] + 1e-30),
        distance=pipeline["peak_distance_bins"],
        prominence=pipeline["log_trace_prominence"],
    )
    indices = candidates[peaks]
    indices = indices[np.argsort(power[indices])[::-1]][:5]
    indices = indices[np.argsort(frequencies[indices])]
    train_directions = leading_directions(train, indices)
    held_out_directions = leading_directions(held_out, indices)
    frozen = manifest["calibration"]

    assert indices.tolist() == frozen["fft_bins"]
    assert np.max(np.abs(train_directions - np.array(frozen["training_directions"]))) < 5e-13
    assert np.max(np.abs(held_out_directions - np.array(frozen["held_out_directions"]))) < 5e-13
    assert np.max(np.abs(train[indices].real - np.array(frozen["training_cospectral_matrices"]))) < 5e-11
    assert np.max(np.abs(held_out[indices].real - np.array(frozen["held_out_cospectral_matrices"]))) < 5e-11
    angles = np.degrees(
        np.arccos(np.clip(np.abs(np.sum(train_directions * held_out_directions, axis=1)), 0, 1))
    )
    print("VBL-VA001 extraction replay: PASS")
    print("frequencies_hz", frequencies[indices].tolist())
    print("train_to_held_out_angles_degrees", angles.tolist())
    print("maximum_angle_degrees", float(angles.max()))
    print("selected cospectral matrices: PASS")


if __name__ == "__main__":
    main()
