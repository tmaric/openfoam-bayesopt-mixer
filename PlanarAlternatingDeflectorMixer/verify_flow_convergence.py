#!/usr/bin/env python3
"""Verify flow convergence without a false failure from near-zero transverse U."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import yaml


def _last_pressure_values(path: Path, log: str, window: int) -> tuple[list[float], int]:
    if not path.exists():
        values = [
            float(value)
            for value in re.findall(r"\bdeltaP=([-+0-9.eE]+)\s+m2/s2", log)
        ]
        if not values:
            raise FileNotFoundError(path)
        return values[-window:], len(values)
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = [float(row["pressure_drop_m2_s2"]) for row in rows]
    return values[-window:], len(values)


def _last_initial_residual(log: str, component: str) -> float:
    pattern = re.compile(
        rf"Solving for {re.escape(component)}, Initial residual = ([-+0-9.eE]+)"
    )
    values = [float(value) for value in pattern.findall(log)]
    if not values:
        raise ValueError(f"no {component} residuals found in flow log")
    return values[-1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("pressure_csv", type=Path)
    parser.add_argument("geometry_yaml", type=Path)
    parser.add_argument("--window", type=int, default=100)
    parser.add_argument("--min-rows", type=int, default=400)
    parser.add_argument("--pressure-span-tol", type=float, default=1e-7)
    parser.add_argument("--residual-tol", type=float, default=1e-5)
    args = parser.parse_args()

    log = args.log.read_text(errors="replace")
    if "SIMPLE solution converged" in log:
        print("[flow-convergence] PASS: SIMPLE residualControl satisfied")
        return

    geometry = yaml.safe_load(args.geometry_yaml.read_text())
    if geometry.get("topology", "alternating_deflector") != "straight":
        raise RuntimeError(
            "SIMPLE residualControl was not satisfied for an obstructed mixer"
        )

    values, total_rows = _last_pressure_values(args.pressure_csv, log, args.window)
    if total_rows < args.min_rows or len(values) < args.window:
        raise RuntimeError(
            f"straight baseline has only {total_rows} pressure rows; need at least "
            f"{max(args.min_rows, args.window)}"
        )
    pressure_span = max(values) - min(values)
    ux_residual = _last_initial_residual(log, "Ux")
    p_residual = _last_initial_residual(log, "p")
    failures = []
    if pressure_span > args.pressure_span_tol:
        failures.append(
            f"final pressure-drop span {pressure_span:.3e} exceeds "
            f"{args.pressure_span_tol:.3e} m2/s2"
        )
    if ux_residual > args.residual_tol:
        failures.append(f"final Ux residual {ux_residual:.3e} exceeds {args.residual_tol:.3e}")
    if p_residual > args.residual_tol:
        failures.append(f"final p residual {p_residual:.3e} exceeds {args.residual_tol:.3e}")
    if failures:
        raise RuntimeError("straight-channel flow convergence failed:\n  " + "\n  ".join(failures))
    print(
        "[flow-convergence] PASS: straight baseline pressure is stable "
        f"(span={pressure_span:.3e}), Ux={ux_residual:.3e}, p={p_residual:.3e}; "
        "near-zero transverse relative residual excluded"
    )


if __name__ == "__main__":
    main()
