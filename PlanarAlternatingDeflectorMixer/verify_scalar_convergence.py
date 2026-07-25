#!/usr/bin/env python3
"""Verify convergence of outlet mixing statistics for a bounded scalar solve."""

import argparse
import csv
import math
from pathlib import Path


def _values(rows: list[dict], column: str) -> list[float]:
    try:
        values = [float(row[column]) for row in rows]
    except KeyError as exc:
        raise ValueError(f"missing required CSV column: {column}") from exc
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"non-finite value in CSV column: {column}")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument("--min-rows", type=int, default=200)
    parser.add_argument("--intensity-tol", type=float, default=1e-4)
    parser.add_argument("--mean-tol", type=float, default=1e-4)
    args = parser.parse_args()
    if args.window < 2 or args.min_rows < args.window:
        raise SystemExit("require min-rows >= window >= 2")
    if args.intensity_tol <= 0 or args.mean_tol <= 0:
        raise SystemExit("convergence tolerances must be positive")

    with args.csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) < max(args.min_rows, args.window):
        raise SystemExit(
            f"scalar objective not converged: only {len(rows)} rows available"
        )

    tail = rows[-args.window :]
    intensity = _values(tail, "flux_weighted_intensity_of_segregation")
    mean = _values(tail, "flux_weighted_mean_concentration")
    if min(intensity) < -1e-10 or max(intensity) > 1.0 + 1e-10:
        raise SystemExit("scalar objective invalid: flux intensity lies outside [0, 1]")
    if min(mean) < -1e-10 or max(mean) > 1.0 + 1e-10:
        raise SystemExit("scalar objective invalid: flux mean lies outside [0, 1]")
    intensity_span = max(intensity) - min(intensity)
    mean_span = max(mean) - min(mean)

    print(
        "[scalar-convergence] "
        f"rows={len(rows)} window={args.window} "
        f"flux_Is_span={intensity_span:.6g} flux_mean_span={mean_span:.6g}"
    )
    if intensity_span > args.intensity_tol or mean_span > args.mean_tol:
        raise SystemExit(
            "scalar objective not converged: final-window variation exceeds tolerance"
        )


if __name__ == "__main__":
    main()
