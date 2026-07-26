#!/usr/bin/env python3
"""Run corrected comparison cases sequentially and validate straight-channel dp."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


CASE_ROOT = Path(__file__).resolve().parent
BASELINE_CONFIG_DIR = CASE_ROOT / "research" / "baselines"
RESULTS_ROOT = CASE_ROOT / "results" / "corrected_boundary_v3_baselines"
CAD_CONFIG_NAME = "alternating_deflector_cad.yaml"
ORDER = ("straight", "symmetric_deflector", "strong_alternating")


def _read_single_row(path: Path) -> dict:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def _run_case(name: str, cores: int) -> None:
    source = BASELINE_CONFIG_DIR / f"{name}.yaml"
    sample_dir = RESULTS_ROOT / name
    sample_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, sample_dir / CAD_CONFIG_NAME)
    command = [
        "snakemake",
        "--snakefile",
        str(CASE_ROOT / "Snakefile"),
        "--directory",
        str(sample_dir),
        "--cores",
        str(cores),
        "--config",
        f"results_dir={sample_dir}",
        f"python_bin={Path(sys.executable).resolve()}",
    ]
    print(f"[baseline] running {name} (strictly sequential, {cores} CFD ranks)")
    subprocess.run(command, check=True)


def _write_summary() -> list[dict]:
    research = yaml.safe_load((CASE_ROOT / "research_config.yaml").read_text())
    operating = research["operating_point"]
    tolerance = float(
        research["validation"]["straight_pressure_drop_relative_tolerance"]
    )
    density = float(operating["fluid_density_kg_m3"])
    viscosity = float(operating["kinematic_viscosity_m2_s"])
    velocity = float(operating["mean_velocity_m_s"])

    rows = []
    for name in ORDER:
        path = RESULTS_ROOT / name / "objectives.csv"
        if not path.exists():
            continue
        row = _read_single_row(path)
        rows.append(
            {
                "baseline": name,
                "pressure_drop_m2_s2": row["pdrop_pressure_drop_m2_s2"],
                "pressure_drop_Pa": row["metric_pressure_drop_Pa"],
                "flux_weighted_intensity_of_segregation": row[
                    "mixing_flux_weighted_intensity_of_segregation"
                ],
                "flux_weighted_mixing_index": row[
                    "metric_flux_weighted_mixing_index"
                ],
                "flow_rate_m3_s": row["metric_flow_rate_m3_s"],
                "pumping_power_W": row["metric_pumping_power_W"],
            }
        )

    straight = next((row for row in rows if row["baseline"] == "straight"), None)
    if straight is not None:
        geometry = yaml.safe_load(
            (BASELINE_CONFIG_DIR / "straight.yaml").read_text()
        )
        height = float(geometry["H"]) * float(geometry["scale"])
        length = (
            2.0 * float(geometry["L0"])
            + int(geometry["N"]) * float(geometry["L_cell"])
        ) * float(geometry["scale"])
        analytical = 12.0 * viscosity * velocity * length / (height * height)
        computed = float(straight["pressure_drop_m2_s2"])
        relative_error = abs(computed - analytical) / analytical
        straight["analytical_pressure_drop_m2_s2"] = analytical
        straight["analytical_pressure_drop_Pa"] = density * analytical
        straight["analytical_relative_error"] = relative_error
        straight["analytical_check_passed"] = relative_error <= tolerance
        if relative_error > tolerance:
            raise RuntimeError(
                "straight-channel pressure validation failed: "
                f"CFD={computed:.8g}, analytical={analytical:.8g} m2/s2, "
                f"relative error={relative_error:.2%} > {tolerance:.2%}"
            )

        straight_pa = float(straight["pressure_drop_Pa"])
        for row in rows:
            row["pressure_ratio_to_straight"] = (
                float(row["pressure_drop_Pa"]) / straight_pa
            )

    if rows:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
        RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
        with open(RESULTS_ROOT / "baseline_summary.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
            writer.writeheader()
            writer.writerows(rows)
        print(f"[baseline] summary -> {RESULTS_ROOT / 'baseline_summary.csv'}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-new-evaluations", type=int, default=1)
    parser.add_argument("--cores", type=int, choices=(1, 2), default=2)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    if args.max_new_evaluations < 1:
        raise ValueError("--max-new-evaluations must be positive")

    if not args.summary_only:
        launched = 0
        for name in ORDER:
            if (RESULTS_ROOT / name / "objectives.csv").exists():
                continue
            if launched >= args.max_new_evaluations:
                break
            _run_case(name, args.cores)
            launched += 1
    rows = _write_summary()
    completed = {row["baseline"] for row in rows}
    print(f"[baseline] completed {len(completed)}/{len(ORDER)}: {sorted(completed)}")


if __name__ == "__main__":
    main()
