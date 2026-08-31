#!/usr/bin/env python3
"""Re-evaluate outlet metrics on existing fields after function-object updates."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess


def last_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"no rows in {path}")
    return rows[-1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", type=Path)
    args = parser.parse_args()

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    for supplied in args.cases:
        case = supplied.resolve()
        if not (case / "system" / "controlDict").is_file():
            raise FileNotFoundError(f"not an OpenFOAM case: {case}")
        log_path = case / "log.postProcess.areaWeighted"
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [
                    "postProcess",
                    "-case",
                    str(case),
                    "-latestTime",
                    "-fields",
                    "(T U phi)",
                ],
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        if completed.returncode != 0:
            tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-30:])
            raise RuntimeError(f"postProcess failed for {case}:\n{tail}")
        row = last_row(case / "mixing.csv")
        report = {
            "case": str(case),
            "time": float(row["time"]),
            "area_weighted_intensity_of_segregation": float(
                row["intensity_of_segregation"]
            ),
            "area_weighted_mixing_index": float(row["mixing_index_rsd"]),
            "flux_weighted_intensity_of_segregation": float(
                row["flux_weighted_intensity_of_segregation"]
            ),
            "flux_weighted_mixing_index": 1.0
            - float(row["flux_weighted_intensity_of_segregation"]) ** 0.5,
        }
        print(json.dumps(report))


if __name__ == "__main__":
    main()
