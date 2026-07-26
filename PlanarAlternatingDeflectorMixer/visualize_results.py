#!/usr/bin/env python3
"""Visualise BO results using the per-sample Python-rendered field PNGs."""

import argparse
import csv
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


CASE_ROOT = Path(__file__).resolve().parent


def load_n_init(default: int = 0) -> int:
    """Read the Sobol initialisation count from the sequential BO config."""
    config_path = CASE_ROOT / "bayes_optimize_sequential.yaml"
    if not config_path.exists():
        return default

    for raw_line in config_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or not line.startswith("n_init:"):
            continue
        try:
            return int(line.split(":", 1)[1].strip())
        except ValueError:
            return default
    return default


def pareto_mask(j1: np.ndarray, j2: np.ndarray) -> np.ndarray:
    """Boolean mask: True where a point is Pareto-optimal (minimise both)."""
    n = len(j1)
    dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if j != i and j1[j] <= j1[i] and j2[j] <= j2[i]:
                if j1[j] < j1[i] or j2[j] < j2[i]:
                    dominated[i] = True
                    break
    return ~dominated


def find_sample_png(
    sample_results_dir: Path,
    sample_id: str,
    extra_search_dirs: tuple[Path, ...] = (),
) -> Path | None:
    """Return the per-sample field PNG, or None if not found."""
    search_dirs = [*extra_search_dirs, sample_results_dir / "visualizations"]
    for directory in search_dirs:
        png_path = directory / f"{sample_id}_T.png"
        if png_path.exists() and png_path.stat().st_size > 0:
            return png_path
    return None


def sample_metrics(sample: dict) -> tuple[float | None, float | None]:
    """Return pressure-drop and mixing metrics when available."""
    try:
        pdrop_raw = sample.get("pdrop_pressure_drop_m2_s2")
        if pdrop_raw in (None, ""):
            # Backward compatibility for the original, dimensionally
            # mislabeled result files. simpleFoam stores kinematic pressure.
            pdrop_raw = sample["pdrop_pressure_drop_Pa"]
        return (
            float(pdrop_raw),
            float(sample["mixing_intensity_of_segregation"]),
        )
    except (KeyError, TypeError, ValueError):
        return None, None


def build_field_title(sample: dict) -> str:
    """Human-readable title for the field panel."""
    sid = sample["sample_id"]
    pdrop, j_mix = sample_metrics(sample)
    if pdrop is None or j_mix is None:
        return f"Concentration field T - sample {sid}"
    return (
        f"Concentration field T - sample {sid} | "
        f"J_dp = {pdrop:.4g} m^2/s^2 | I_s = {j_mix:.4f}"
    )


def draw_missing_panel(ax, title: str, message: str) -> None:
    """Draw a placeholder panel when no field PNG is available."""
    ax.axis("off")
    ax.set_title(title, fontsize=9)
    ax.text(
        0.5,
        0.5,
        message,
        transform=ax.transAxes,
        ha="center",
        va="center",
        color="gray",
    )


def collect_per_sample_pngs(
    samples: list[dict], results_root: Path, output_dir: Path
) -> None:
    """Copy already-rendered per-sample PNGs into the requested output dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for sample in samples:
        sid = sample["sample_id"]
        sample_dir = results_root / sid
        src = find_sample_png(sample_dir, sid)
        if src is None:
            print(
                f"  WARNING: {sid}: field PNG missing - animation panel will be blank",
                file=sys.stderr,
            )
            continue

        dst = output_dir / f"{sid}_T.png"
        if src.resolve() != dst.resolve():
            shutil.copy2(src, dst)
        print(f"  {sid}: saved {dst.name}")


def create_pareto_animation(
    samples: list[dict], results_root: Path, output_dir: Path
) -> None:
    """Animated GIF/MP4: Pareto front scatter plus rendered T-field panel."""
    output_dir.mkdir(parents=True, exist_ok=True)

    j_dp = np.array([sample_metrics(sample)[0] for sample in samples], dtype=float)
    j_mix = np.array([
        float(sample["mixing_intensity_of_segregation"])
        for sample in samples
    ])
    n_init = load_n_init()

    fig = plt.figure(figsize=(11, 8))
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 1.1], hspace=0.40)
    ax_pareto = fig.add_subplot(grid[0])
    ax_field = fig.add_subplot(grid[1])

    def update(frame: int) -> None:
        ax_pareto.cla()
        ax_field.cla()

        visible = frame + 1
        j_dp_k = j_dp[:visible]
        j_mix_k = j_mix[:visible]
        pareto = pareto_mask(j_dp_k, j_mix_k)
        point_colors = np.array(
            ["steelblue" if i < n_init else "darkorange" for i in range(visible)],
            dtype=object,
        )

        non_pareto = ~pareto
        if non_pareto.any():
            ax_pareto.scatter(
                j_dp_k[non_pareto],
                j_mix_k[non_pareto],
                c=point_colors[non_pareto].tolist(),
                alpha=0.6,
                s=40,
            )
        if pareto.any():
            ax_pareto.scatter(
                j_dp_k[pareto],
                j_mix_k[pareto],
                c=point_colors[pareto].tolist(),
                s=90,
                edgecolors="black",
                linewidths=1.4,
                zorder=5,
            )

        ax_pareto.scatter(
            j_dp_k[-1],
            j_mix_k[-1],
            c="red",
            s=100,
            marker="D",
            zorder=6,
        )
        if pareto.any():
            pidx = np.where(pareto)[0]
            order = pidx[np.argsort(j_dp_k[pidx])]
            ax_pareto.plot(j_dp_k[order], j_mix_k[order], color="black", lw=1.2, zorder=4)

        ax_pareto.set_xlabel(
            "Kinematic pressure drop  $J_{dp}$  [m$^2$/s$^2$, log scale]"
        )
        ax_pareto.set_ylabel("Intensity of segregation  $J_{mix}$  [ ]")
        ax_pareto.set_xscale("log")
        ax_pareto.set_title(
            f"Pareto front - sample {samples[frame]['sample_id']} ({visible} / {len(samples)})"
        )
        ax_pareto.grid(True, which="both", alpha=0.3)
        ax_pareto.legend(
            handles=[
                Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=8,
                       label=f"Sobol init ({min(n_init, visible)})"),
                Line2D([0], [0], marker="o", color="w", markerfacecolor="darkorange", markersize=8,
                       label=f"BO suggested ({max(0, visible - n_init)})"),
                Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
                       markeredgecolor="black", markeredgewidth=1.4, markersize=8,
                       label=f"Pareto ({int(pareto.sum())})"),
                Line2D([0], [0], marker="D", color="w", markerfacecolor="red", markersize=8,
                       label="Current"),
            ],
            loc="upper right",
            fontsize=8,
        )

        sample = samples[frame]
        sid = sample["sample_id"]
        sample_dir = results_root / sid
        title = build_field_title(sample)

        png_path = find_sample_png(sample_dir, sid, extra_search_dirs=(output_dir,))
        if png_path is not None:
            image = plt.imread(str(png_path))
            ax_field.imshow(image, aspect="auto")
            ax_field.axis("off")
            ax_field.set_title(title, fontsize=9)
        else:
            draw_missing_panel(ax_field, title, "Field PNG not available")

    ani = animation.FuncAnimation(fig, update, frames=len(samples), interval=700)

    gif_path = output_dir / "pareto_animation.gif"
    ani.save(str(gif_path), writer=animation.PillowWriter(fps=1.5))
    print(f"  Saved {gif_path}")

    mp4_path = output_dir / "pareto_animation.mp4"
    try:
        ani.save(str(mp4_path), writer=animation.FFMpegWriter(fps=1.5, bitrate=1800))
        print(f"  Saved {mp4_path}")
    except Exception as exc:
        print(
            f"  WARNING: MP4 export failed ({exc}) - install ffmpeg for MP4 output",
            file=sys.stderr,
        )

    plt.close(fig)


def copy_single_sample_png(results_root: Path, output_dir: Path, sample_id: str) -> None:
    """Copy one already-rendered field PNG into *output_dir* when present."""
    output_dir.mkdir(parents=True, exist_ok=True)
    src = find_sample_png(results_root, sample_id, extra_search_dirs=(output_dir,))
    if src is None:
        print(
            f"WARNING: {sample_id}: field PNG not found under {results_root / 'visualizations'}",
            file=sys.stderr,
        )
        return

    dst = output_dir / f"{sample_id}_T.png"
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    print(f"  {sample_id}: saved {dst.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Root results directory (contains all_samples.csv and sample subdirs)",
    )
    parser.add_argument(
        "--output-dir",
        help="Where to write visualizations (default: <results-dir>/visualizations)",
    )
    parser.add_argument(
        "--sample-id",
        help="Copy only the rendered field PNG for one sample ID and exit.",
    )
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else results_root / "visualizations"

    if args.sample_id:
        copy_single_sample_png(results_root, output_dir, args.sample_id)
        print()
        print(f"Done. Output in {output_dir}/")
        return

    all_csv = results_root / "all_samples.csv"
    if not all_csv.exists():
        sys.exit(f"ERROR: {all_csv} not found - run the BO loop first")

    with open(all_csv, newline="") as handle:
        samples = list(csv.DictReader(handle))

    samples.sort(key=lambda sample: sample["sample_id"])
    print(f"Loaded {len(samples)} samples from {all_csv}")

    print()
    print("[1/2] Collecting per-sample field PNGs ...")
    collect_per_sample_pngs(samples, results_root, output_dir)

    print()
    print("[2/2] Animated Pareto front ...")
    create_pareto_animation(samples, results_root, output_dir)

    print()
    print(f"Done. Output in {output_dir}/")


if __name__ == "__main__":
    main()
