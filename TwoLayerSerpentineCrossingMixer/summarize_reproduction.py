#!/usr/bin/env python3
"""Summarize M10-source benchmark cases for the inspired reconstruction."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "reproduction"
REVIEW_TARGETS = {
    1.0: {"mixing_index": 0.915, "pressure_drop_Pa": 16.3},
    20.0: {"mixing_index": 0.901, "pressure_drop_Pa": 390.0},
}


def last_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"no rows in {path}")
    return rows[-1]


def collect() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for objective_path in RESULTS.glob("*/Re*/*/objectives.json"):
        fidelity = objective_path.parent.name
        reynolds_label = objective_path.parent.parent.name
        protocol = objective_path.parent.parent.parent.name
        objectives = json.loads(objective_path.read_text(encoding="utf-8"))
        reynolds = float(reynolds_label.removeprefix("Re"))
        mixing_path = objective_path.parent / "ScalarTransportCase" / "mixing.csv"
        mixing = last_csv_row(mixing_path)
        area_intensity = float(mixing["intensity_of_segregation"])
        record: dict[str, object] = {
            "protocol": protocol,
            "reynolds_number": reynolds,
            "fidelity": fidelity,
            "pressure_drop_Pa": float(objectives["pressure_drop_Pa"]),
            "flux_weighted_intensity_of_segregation": float(
                objectives["flux_weighted_intensity_of_segregation"]
            ),
            "flux_weighted_mixing_index": float(objectives["flux_weighted_mixing_index"]),
            "area_weighted_intensity_of_segregation": area_intensity,
            "area_weighted_mixing_index": 1.0 - math.sqrt(area_intensity),
        }
        target = (
            REVIEW_TARGETS.get(reynolds)
            if protocol in {"review", "review_second_order"}
            else None
        )
        if target:
            record["target_area_weighted_mixing_index"] = target["mixing_index"]
            record["area_weighted_mixing_index_error"] = (
                record["area_weighted_mixing_index"] - target["mixing_index"]
            )
            record["target_pressure_drop_Pa"] = target["pressure_drop_Pa"]
            record["pressure_drop_relative_error"] = (
                record["pressure_drop_Pa"] / target["pressure_drop_Pa"] - 1.0
            )
        records.append(record)
    return sorted(
        records,
        key=lambda row: (
            str(row["protocol"]),
            float(row["reynolds_number"]),
            0 if row["fidelity"] == "coarse" else 1,
        ),
    )


def markdown(records: list[dict[str, object]]) -> str:
    lines = [
        "| Protocol | Re | Fidelity | Pressure (Pa) | Area MI | Flux MI | Flux segregation |",
        "|---|---:|---|---:|---:|---:|---:|",
    ]
    for row in records:
        lines.append(
            "| {protocol} | {reynolds_number:g} | {fidelity} | {pressure_drop_Pa:.3f} | "
            "{area_weighted_mixing_index:.6f} | {flux_weighted_mixing_index:.6f} | "
            "{flux_weighted_intensity_of_segregation:.6g} |".format(**row)
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Print JSON instead of Markdown")
    args = parser.parse_args()
    records = collect()
    print(json.dumps(records, indent=2) if args.json else markdown(records))


if __name__ == "__main__":
    main()
