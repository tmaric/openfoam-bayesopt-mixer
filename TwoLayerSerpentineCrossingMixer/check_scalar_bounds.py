#!/usr/bin/env python3
"""Run OpenFOAM fieldMinMax on converged scalar cases without ParaView."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", type=Path)
    args = parser.parse_args()

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    pattern = re.compile(
        r"min\(T\)\s*=\s*([0-9.eE+-]+).*?max\(T\)\s*=\s*([0-9.eE+-]+)",
        flags=re.DOTALL,
    )
    for supplied in args.cases:
        case = supplied.resolve()
        completed = subprocess.run(
            [
                "postProcess",
                "-case",
                str(case),
                "-latestTime",
                "-field",
                "T",
                "-func",
                "fieldMinMax(T)",
            ],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        log_path = case / "log.postProcess.scalarBounds"
        log_path.write_text(completed.stdout, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(
                f"fieldMinMax failed for {case}:\n"
                + "\n".join(completed.stdout.splitlines()[-30:])
            )
        matches = pattern.findall(completed.stdout)
        if not matches:
            raise RuntimeError(f"T bounds not found in {log_path}")
        minimum, maximum = map(float, matches[-1])
        report = {
            "case": str(case),
            "minimum_T": minimum,
            "maximum_T": maximum,
            "undershoot": max(0.0, -minimum),
            "overshoot": max(0.0, maximum - 1.0),
        }
        (case.parent / "scalar_bounds.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(report))


if __name__ == "__main__":
    main()
