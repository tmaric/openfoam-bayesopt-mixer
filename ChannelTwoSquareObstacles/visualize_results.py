#!/usr/bin/env python3
"""Visualize hydro-only BO results with pressure-drop history and pressure-field PNGs."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


CASE_ROOT = Path(__file__).resolve().parent


def find_sample_png(
    sample_results_dir: Path,
    sample_id: str,
    extra_search_dirs: tuple[Path, ...] = (),
) -> Path | None:
    search_dirs = [*extra_search_dirs, sample_results_dir / "visualizations"]
    for directory in search_dirs:
        png_path = directory / f"{sample_id}_p.png"
        if png_path.exists() and png_path.stat().st_size > 0:
            return png_path
    return None


def successful_samples(all_samples: list[dict]) -> list[dict]:
    successful = []
    for sample in all_samples:
        try:
            sample["pdrop_value"] = float(sample["pdrop_pressure_drop_Pa"])
        except (KeyError, TypeError, ValueError):
            continue
        successful.append(sample)
    successful.sort(key=lambda sample: int(sample["sample_id"]))
    return successful


def build_field_title(sample: dict) -> str:
    return f"Pressure field p - sample {sample['sample_id']} | J_dp = {sample['pdrop_value']:.4g} Pa"


def draw_missing_panel(ax, title: str, message: str) -> None:
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


def collect_per_sample_pngs(samples: list[dict], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for sample in samples:
        sid = sample["sample_id"]
        src = find_sample_png(Path(sample["results_dir"]), sid)
        if src is None:
            print(
                f"  WARNING: {sid}: pressure PNG missing - animation panel will be blank",
                file=sys.stderr,
            )
            continue

        dst = output_dir / f"{sid}_p.png"
        if src.resolve() != dst.resolve():
            shutil.copy2(src, dst)
        print(f"  {sid}: saved {dst.name}")


def create_pressure_animation(samples: list[dict], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = np.array([int(sample["sample_id"]) for sample in samples], dtype=int)
    pdrop = np.array([sample["pdrop_value"] for sample in samples], dtype=float)
    phases = np.array(
        [
            str(sample.get("phase", "")).strip().lower() or "bo"
            for sample in samples
        ],
        dtype=object,
    )
    colors = np.where(phases == "sobol", "steelblue", "darkorange")
    best_so_far = np.maximum.accumulate(pdrop)

    fig = plt.figure(figsize=(11, 8))
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 1.1], hspace=0.40)
    ax_history = fig.add_subplot(grid[0])
    ax_field = fig.add_subplot(grid[1])

    def update(frame: int) -> None:
        ax_history.cla()
        ax_field.cla()

        visible = frame + 1
        visible_ids = sample_ids[:visible]
        visible_pdrop = pdrop[:visible]
        visible_best = best_so_far[:visible]
        visible_colors = colors[:visible]

        sobol_mask = visible_colors == "steelblue"
        bo_mask = visible_colors == "darkorange"
        if sobol_mask.any():
            ax_history.scatter(
                visible_ids[sobol_mask],
                visible_pdrop[sobol_mask],
                c="steelblue",
                s=45,
                alpha=0.8,
            )
        if bo_mask.any():
            ax_history.scatter(
                visible_ids[bo_mask],
                visible_pdrop[bo_mask],
                c="darkorange",
                s=45,
                alpha=0.8,
            )

        ax_history.plot(
            visible_ids,
            visible_best,
            color="black",
            lw=1.5,
            label="Best so far",
            zorder=4,
        )
        ax_history.scatter(
            visible_ids[-1],
            visible_pdrop[-1],
            c="red",
            s=100,
            marker="D",
            zorder=5,
        )
        ax_history.set_xlabel("Sample index")
        ax_history.set_ylabel("Pressure drop  $J_{dp}$  [Pa]")
        ax_history.set_title(
            f"Pressure-drop history - sample {samples[frame]['sample_id']} ({visible} / {len(samples)})"
        )
        ax_history.grid(True, alpha=0.3)
        ax_history.legend(
            handles=[
                Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=8, label="Sobol init"),
                Line2D([0], [0], marker="o", color="w", markerfacecolor="darkorange", markersize=8, label="BO suggested"),
                Line2D([0], [0], color="black", lw=1.5, label="Best so far"),
                Line2D([0], [0], marker="D", color="w", markerfacecolor="red", markersize=8, label="Current"),
            ],
            loc="best",
            fontsize=8,
        )

        sample = samples[frame]
        title = build_field_title(sample)
        png_path = find_sample_png(Path(sample["results_dir"]), sample["sample_id"], extra_search_dirs=(output_dir,))
        if png_path is not None:
            image = plt.imread(str(png_path))
            ax_field.imshow(image, aspect="auto")
            ax_field.axis("off")
            ax_field.set_title(title, fontsize=9)
        else:
            draw_missing_panel(ax_field, title, "ParaView pressure PNG not available")

    ani = animation.FuncAnimation(fig, update, frames=len(samples), interval=700)

    gif_path = output_dir / "pressure_drop_animation.gif"
    ani.save(str(gif_path), writer=animation.PillowWriter(fps=1.5))
    print(f"  Saved {gif_path}")

    mp4_path = output_dir / "pressure_drop_animation.mp4"
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
    output_dir.mkdir(parents=True, exist_ok=True)
    src = find_sample_png(results_root, sample_id, extra_search_dirs=(output_dir,))
    if src is None:
        print(
            f"WARNING: {sample_id}: pressure PNG not found under {results_root / 'visualizations'}",
            file=sys.stderr,
        )
        return

    dst = output_dir / f"{sample_id}_p.png"
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    print(f"  {sample_id}: saved {dst.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, help="Root results directory")
    parser.add_argument("--output-dir", help="Where to write visualizations (default: <results-dir>/visualizations)")
    parser.add_argument("--sample-id", help="Copy only the rendered pressure PNG for one sample ID and exit.")
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
        all_samples = list(csv.DictReader(handle))

    samples = successful_samples(all_samples)
    if not samples:
        sys.exit("ERROR: no successful samples with pressure-drop values were found")

    print(f"Loaded {len(samples)} successful samples from {all_csv}")

    print()
    print("[1/2] Collecting per-sample ParaView PNGs ...")
    collect_per_sample_pngs(samples, output_dir)

    print()
    print("[2/2] Pressure-drop animation ...")
    create_pressure_animation(samples, output_dir)

    print()
    print(f"Done. Output in {output_dir}/")


if __name__ == "__main__":
    main()
