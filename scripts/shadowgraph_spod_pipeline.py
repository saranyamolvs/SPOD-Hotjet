#!/usr/bin/env python3
"""Batch preprocessing and SPOD analysis for impinging-jet shadowgraph TIFF image sets.

Expected input layout by default:
    <input-root>/hD2/NPR2p5/*.tif
    <input-root>/hD2/NPR3p67/*.tif
    <input-root>/hD2/NPR5/*.tif
    ...

The script writes publication-oriented outputs under a separate results directory,
including enhanced image examples, mean/RMS fields, SPOD eigenvalue spectra,
mode images, coefficient time histories, and machine-readable CSV/NPZ files.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence


TIFF_EXTENSIONS = ("*.tif", "*.tiff", "*.TIF", "*.TIFF")


@dataclass(frozen=True)
class CaseInfo:
    h_over_d: str
    npr: str
    input_dir: str
    output_dir: str
    frame_count: int
    image_shape: tuple[int, int]
    dt: float
    preprocessing: dict[str, float | int | str]
    spod: dict[str, float | int | str]


def natural_key(path: Path) -> list[int | str]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", path.name)]


def discover_cases(input_root: Path, hds: Sequence[str] | None, nprs: Sequence[str] | None) -> list[Path]:
    hd_dirs = [p for p in input_root.iterdir() if p.is_dir() and p.name.lower().startswith("hd")]
    if hds:
        wanted = {normalize_hd(hd).lower() for hd in hds}
        hd_dirs = [p for p in hd_dirs if p.name.lower() in wanted]
    cases: list[Path] = []
    for hd_dir in sorted(hd_dirs, key=natural_key):
        npr_dirs = [p for p in hd_dir.iterdir() if p.is_dir() and p.name.lower().startswith(("npr", "npt"))]
        if nprs:
            wanted_nprs = {normalize_npr(npr).lower() for npr in nprs}
            npr_dirs = [p for p in npr_dirs if canonical_npr_name(p.name).lower() in wanted_nprs]
        cases.extend(sorted(npr_dirs, key=natural_key))
    return cases


def normalize_hd(value: str) -> str:
    value = value.strip()
    return value if value.lower().startswith("hd") else f"hD{value}"


def normalize_npr(value: str) -> str:
    value = value.strip().replace(".", "p")
    return canonical_npr_name(value if value.lower().startswith(("npr", "npt")) else f"NPR{value}")


def canonical_npr_name(name: str) -> str:
    return "NPR" + name[3:] if name.lower().startswith("npt") else name


def image_paths(case_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for pattern in TIFF_EXTENSIONS:
        paths.extend(case_dir.glob(pattern))
    return sorted(paths, key=natural_key)


def load_stack(paths: Sequence[Path], max_frames: int | None = None) -> np.ndarray:
    import imageio.v3 as iio
    import numpy as np

    selected = paths[:max_frames] if max_frames else paths
    frames = []
    for path in selected:
        img = iio.imread(path)
        if img.ndim == 3:
            img = img[..., :3].mean(axis=2)
        frames.append(img.astype(np.float32))
    if not frames:
        raise ValueError("No TIFF frames were loaded")
    return np.stack(frames, axis=0)


def percentile_normalize(img: np.ndarray, low: float, high: float) -> np.ndarray:
    import numpy as np

    lo, hi = np.percentile(img, [low, high])
    if not np.isfinite(hi - lo) or hi <= lo:
        return np.zeros_like(img, dtype=np.float32)
    return np.clip((img - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def preprocess_stack(stack: np.ndarray, clip_low: float, clip_high: float, flat_sigma: float, gamma: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    import numpy as np
    from scipy.ndimage import gaussian_filter

    raw_mean = stack.mean(axis=0)
    flat = gaussian_filter(raw_mean, sigma=flat_sigma) if flat_sigma > 0 else raw_mean
    corrected = (stack - flat) / (np.abs(flat).mean() + 1e-6)
    enhanced = np.empty_like(corrected, dtype=np.float32)
    for idx, frame in enumerate(corrected):
        norm = percentile_normalize(frame, clip_low, clip_high)
        enhanced[idx] = np.power(norm, gamma, dtype=np.float32)
    fluctuations = enhanced - enhanced.mean(axis=0, keepdims=True)
    return enhanced, fluctuations, raw_mean, flat


def compute_spod(fluctuations: np.ndarray, dt: float, nperseg: int, overlap: float, max_modes: int) -> dict[str, np.ndarray]:
    import numpy as np
    from scipy import signal

    nt, ny, nx = fluctuations.shape
    nperseg = min(nperseg, nt)
    if nperseg < 4:
        raise ValueError("SPOD requires at least 4 frames per segment")
    noverlap = int(round(nperseg * overlap))
    noverlap = min(max(noverlap, 0), nperseg - 1)
    step = nperseg - noverlap
    starts = list(range(0, nt - nperseg + 1, step))
    if not starts:
        starts = [0]
    window = signal.windows.hann(nperseg, sym=False).astype(np.float32)
    window_norm = np.sqrt(np.mean(window**2))
    freqs = np.fft.rfftfreq(nperseg, d=dt)
    n_freq = len(freqs)
    n_blocks = len(starts)
    n_pixels = ny * nx
    qhat = np.empty((n_freq, n_blocks, n_pixels), dtype=np.complex64)
    flat = fluctuations.reshape(nt, n_pixels)
    for block_idx, start in enumerate(starts):
        block = flat[start : start + nperseg]
        block = (block - block.mean(axis=0, keepdims=True)) * (window[:, None] / window_norm)
        qhat[:, block_idx, :] = np.fft.rfft(block, axis=0) / math.sqrt(nperseg)

    n_modes = min(max_modes, n_blocks)
    eigvals = np.zeros((n_freq, n_modes), dtype=np.float32)
    modes = np.zeros((n_freq, n_modes, ny, nx), dtype=np.float32)
    coeffs = np.zeros((n_freq, n_modes, n_blocks), dtype=np.complex64)
    for fi in range(n_freq):
        x = qhat[fi]
        csd = (x @ x.conj().T) / n_blocks
        vals, vecs = np.linalg.eigh(csd)
        order = vals.argsort()[::-1]
        vals = np.maximum(vals[order], 0.0)
        vecs = vecs[:, order]
        keep = min(n_modes, len(vals))
        eigvals[fi, :keep] = vals[:keep].real.astype(np.float32)
        for mi in range(keep):
            if vals[mi] > 0:
                mode = (x.conj().T @ vecs[:, mi]) / math.sqrt(vals[mi] * n_blocks)
                modes[fi, mi] = np.real(mode).reshape(ny, nx).astype(np.float32)
                coeffs[fi, mi] = np.sqrt(vals[mi]) * vecs[:, mi]
    return {"frequencies": freqs.astype(np.float32), "eigenvalues": eigvals, "modes": modes, "coefficients": coeffs, "starts": np.asarray(starts)}


def save_image(path: Path, data: np.ndarray, cmap: str = "magma", title: str | None = None) -> None:
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 5), constrained_layout=True)
    plt.imshow(data, cmap=cmap, origin="upper")
    plt.colorbar(fraction=0.046, pad=0.04)
    if title:
        plt.title(title)
    plt.axis("off")
    plt.savefig(path, dpi=300)
    plt.close()


def write_eigen_csv(path: Path, freqs: np.ndarray, eigvals: np.ndarray) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frequency_hz", *[f"mode_{i+1}_energy" for i in range(eigvals.shape[1])]])
        for freq, row in zip(freqs, eigvals):
            writer.writerow([freq, *row])


def plot_spectrum(path: Path, freqs: np.ndarray, eigvals: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(7, 5), constrained_layout=True)
    for idx in range(eigvals.shape[1]):
        plt.semilogy(freqs, eigvals[:, idx] + 1e-30, label=f"Mode {idx + 1}")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("SPOD eigenvalue / modal energy")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.savefig(path, dpi=300)
    plt.close()


def process_case(case_dir: Path, output_root: Path, args: argparse.Namespace) -> CaseInfo:
    import numpy as np

    paths = image_paths(case_dir)
    stack = load_stack(paths, args.max_frames)
    enhanced, fluctuations, raw_mean, flat = preprocess_stack(stack, args.clip_low, args.clip_high, args.flat_sigma, args.gamma)
    hd_name = case_dir.parent.name
    npr_name = canonical_npr_name(case_dir.name)
    out = output_root / hd_name / npr_name
    out.mkdir(parents=True, exist_ok=True)

    save_image(out / "mean_raw.png", raw_mean, "gray", "Raw mean shadowgraph")
    save_image(out / "flat_field_background.png", flat, "gray", "Estimated background")
    save_image(out / "mean_enhanced.png", enhanced.mean(axis=0), "magma", "Enhanced mean")
    save_image(out / "rms_fluctuation.png", fluctuations.std(axis=0), "viridis", "RMS fluctuation")
    for idx in sorted(set([0, len(enhanced) // 2, len(enhanced) - 1])):
        save_image(out / f"enhanced_frame_{idx:05d}.png", enhanced[idx], "magma", f"Enhanced frame {idx}")

    spod = compute_spod(fluctuations, args.dt, args.nperseg, args.overlap, args.max_modes)
    np.savez_compressed(out / "spod_results.npz", **spod)
    write_eigen_csv(out / "spod_eigenvalues.csv", spod["frequencies"], spod["eigenvalues"])
    plot_spectrum(out / "spod_spectrum.png", spod["frequencies"], spod["eigenvalues"])

    ranked = np.argsort(spod["eigenvalues"][:, 0])[::-1]
    for rank, fi in enumerate(ranked[: args.save_top_frequencies], start=1):
        freq_dir = out / f"top_frequency_{rank:02d}_{spod['frequencies'][fi]:.6g}Hz"
        freq_dir.mkdir(exist_ok=True)
        for mi in range(spod["modes"].shape[1]):
            save_image(freq_dir / f"spod_mode_{mi + 1}.png", spod["modes"][fi, mi], "RdBu_r", f"SPOD mode {mi + 1}, f={spod['frequencies'][fi]:.6g} Hz")

    info = CaseInfo(hd_name, npr_name, str(case_dir), str(out), int(stack.shape[0]), tuple(map(int, stack.shape[1:])), args.dt,
                    {"clip_low": args.clip_low, "clip_high": args.clip_high, "flat_sigma": args.flat_sigma, "gamma": args.gamma},
                    {"nperseg": min(args.nperseg, stack.shape[0]), "overlap": args.overlap, "max_modes": args.max_modes})
    (out / "case_summary.json").write_text(json.dumps(asdict(info), indent=2))
    return info


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preprocess shadowgraph TIFFs and compute SPOD for all h/D-NPR cases.")
    parser.add_argument("--input-root", type=Path, default=Path("RP"), help="Root containing hD*/NPR* TIFF folders.")
    parser.add_argument("--output-root", type=Path, default=Path("RP") / "SPODResultsCodex" / "results", help="Separate directory for all outputs.")
    parser.add_argument("--fps", type=float, default=22000.0, help="Camera frame rate in frames per second; default is 22000 fps.")
    parser.add_argument("--dt", type=float, help="Optional time between frames in seconds. Overrides --fps when provided.")
    parser.add_argument("--hds", nargs="*", help="Optional subset, e.g. 2 4 6 8 10 12 20 or hD2 hD4.")
    parser.add_argument("--nprs", nargs="*", help="Optional subset, e.g. 2.5 3.67 5.0 or NPR2p5.")
    parser.add_argument("--max-frames", type=int, help="Optional limit for trial runs.")
    parser.add_argument("--clip-low", type=float, default=1.0, help="Low percentile for contrast stretching.")
    parser.add_argument("--clip-high", type=float, default=99.0, help="High percentile for contrast stretching.")
    parser.add_argument("--flat-sigma", type=float, default=35.0, help="Gaussian sigma for mean-image flat-field/background estimate.")
    parser.add_argument("--gamma", type=float, default=0.8, help="Gamma for enhanced display images (<1 brightens weak structures).")
    parser.add_argument("--nperseg", type=int, default=256, help="Frames per Welch/SPOD block.")
    parser.add_argument("--overlap", type=float, default=0.5, help="Fractional block overlap, 0 to <1.")
    parser.add_argument("--max-modes", type=int, default=6, help="SPOD modes saved at every frequency.")
    parser.add_argument("--save-top-frequencies", type=int, default=5, help="Number of energetic frequencies for mode PNG exports.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.dt is None:
        args.dt = 1.0 / args.fps
    if not args.input_root.exists():
        raise SystemExit(f"Input root not found: {args.input_root}")
    cases = discover_cases(args.input_root, args.hds, args.nprs)
    if not cases:
        raise SystemExit(f"No hD*/NPR* cases found under {args.input_root}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    summaries = []
    for case_dir in cases:
        print(f"Processing {case_dir} ...", flush=True)
        summaries.append(asdict(process_case(case_dir, args.output_root, args)))
    (args.output_root / "run_summary.json").write_text(json.dumps(summaries, indent=2))
    print(f"Done. Results written to {args.output_root}")


if __name__ == "__main__":
    main()
