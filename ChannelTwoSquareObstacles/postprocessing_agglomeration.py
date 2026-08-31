#!/usr/bin/env python3
"""Agglomerate one hydro-only BO sample into a single-row objectives CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml


def last_row(csv_path: Path) -> dict:
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        row = None
        for row in reader:
            pass
    if row is None:
        raise ValueError(f"No data rows found in {csv_path}")
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--geometry-manifest",
        required=True,
        type=Path,
        help="Path to the fully expanded geometry manifest YAML.",
    )
    parser.add_argument("--pdrop", required=True, type=Path, help="Path to pressureDrop.csv.")
    parser.add_argument("--output", required=True, type=Path, help="Destination objectives.csv.")
    parser.add_argument("--results-dir", required=True, type=Path, help="Root results directory for this sample.")
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    sample_id = results_dir.name

    with open(args.geometry_manifest, encoding="utf-8") as handle:
        geometry_manifest = yaml.safe_load(handle)

    geo_params = {
        f"geo_{key}": value
        for key, value in geometry_manifest.items()
        if isinstance(value, (int, float, str, bool))
    }

    pdrop_row = last_row(args.pdrop)
    pdrop_cols = {f"pdrop_{key}": value for key, value in pdrop_row.items()}

    row = {
        "sample_id": sample_id,
        "results_dir": str(results_dir),
        "cad_mode": geometry_manifest.get("cad_mode", ""),
        "feasible": int(bool(geometry_manifest.get("feasible", True))),
        **geo_params,
        **pdrop_cols,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    print(f"[agglomerate] written: {args.output}")


if __name__ == "__main__":
    main()
