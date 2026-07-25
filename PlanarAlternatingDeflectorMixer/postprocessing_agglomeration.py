#!/usr/bin/env python3
"""
Agglomerate one BO sample into a single-row objectives CSV.

Columns produced:
    sample_id        -- name of the results directory (set by the BO loop)
    results_dir      -- absolute path to the sample results directory
    geo_<param>      -- geometry parameters from alternating_deflector_cad.yaml
    pdrop_<col>      -- all columns from pressureDrop.csv (last converged row)
    mixing_<col>     -- all columns from mixing.csv (last converged row)

Usage:
    python postprocessing_agglomeration.py \
        --yaml        <results>/FlowCase/alternating_deflector_cad.yaml \
        --pdrop       <results>/FlowCase/pressureDrop.csv                \
        --mixing      <results>/ScalarTransportCase/mixing.csv           \
        --output      <results>/objectives.csv                            \
        --results-dir <results>
"""

import argparse
import csv
from pathlib import Path

import yaml


def last_row(csv_path: Path) -> dict:
    """Return the last data row of a CSV as an ordered dict."""
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        row = None
        for row in reader:
            pass
    if row is None:
        raise ValueError(f"No data rows found in {csv_path}")
    return row


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Agglomerate geometry + CFD objectives for one BO sample."
    )
    parser.add_argument("--yaml",        required=True, type=Path,
                        help="Path to alternating_deflector_cad.yaml in the staged flow case.")
    parser.add_argument("--pdrop",       required=True, type=Path,
                        help="Path to pressureDrop.csv.")
    parser.add_argument("--mixing",      required=True, type=Path,
                        help="Path to mixing.csv.")
    parser.add_argument("--output",      required=True, type=Path,
                        help="Destination objectives.csv.")
    parser.add_argument("--results-dir", required=True, type=Path,
                        help="Root results directory for this sample.")
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    sample_id   = results_dir.name

    # --- geometry parameters ------------------------------------------------
    with open(args.yaml) as fh:
        raw = yaml.safe_load(fh)
    raw.setdefault("k", 0.0)
    # Historical sample YAMLs predate the linear-delta slope parameter.
    # Normalise them to k = 0.0 so re-agglomerated result tables stay schema-stable.
    geo_params = {
        f"geo_{k}": v for k, v in raw.items()
        if isinstance(v, (int, float, str, bool))
    }

    # --- last converged row from each solver CSV ----------------------------
    pdrop_row  = last_row(args.pdrop)
    mixing_row = last_row(args.mixing)

    pdrop_cols  = {f"pdrop_{k}":  v for k, v in pdrop_row.items()}
    mixing_cols = {f"mixing_{k}": v for k, v in mixing_row.items()}

    # --- assemble and write -------------------------------------------------
    row = {
        "sample_id":   sample_id,
        # The sample directory name is portable; callers already know the
        # results root and should not persist a machine-specific absolute path.
        "results_dir": sample_id,
        **geo_params,
        **pdrop_cols,
        **mixing_cols,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    print(f"[agglomerate] written: {args.output}")


if __name__ == "__main__":
    main()
