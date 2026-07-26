#!/usr/bin/env python3
"""Advance or inspect the corrected research sequence one bounded stage at a time."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


CASE_ROOT = Path(__file__).resolve().parent
BASELINE_SUMMARY = (
    CASE_ROOT / "results" / "corrected_boundary_v3_baselines" / "baseline_summary.csv"
)
CAMPAIGN_ROOT = CASE_ROOT / "results" / "corrected_boundary_v3"
REQUIRED_BASELINES = {"straight", "symmetric_deflector", "strong_alternating"}


def _baseline_names() -> set[str]:
    if not BASELINE_SUMMARY.exists():
        return set()
    with open(BASELINE_SUMMARY, newline="") as handle:
        return {row["baseline"] for row in csv.DictReader(handle)}


def _campaign_counts() -> tuple[int, int]:
    successes = 0
    failures = 0
    if not CAMPAIGN_ROOT.exists():
        return successes, failures
    for sample_dir in CAMPAIGN_ROOT.iterdir():
        path = sample_dir / "objectives.csv"
        if not sample_dir.name.isdigit() or not path.exists():
            continue
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            continue
        row = rows[-1]
        if (
            str(row.get("failed", "false")).lower() == "true"
            or row.get("mixing_flux_weighted_intensity_of_segregation", "") == ""
        ):
            failures += 1
        else:
            successes += 1
    return successes, failures


def _status() -> dict:
    baselines = _baseline_names()
    successes, failures = _campaign_counts()
    gate = None
    gate_path = CAMPAIGN_ROOT / "screening_gate.json"
    if gate_path.exists():
        gate = json.loads(gate_path.read_text())
    return {
        "baselines_complete": sorted(baselines),
        "baselines_missing": sorted(REQUIRED_BASELINES.difference(baselines)),
        "screening_successes": successes,
        "screening_failures": failures,
        "screening_target": 12,
        "screening_gate": gate,
    }


def _print_status(status: dict) -> None:
    print("Corrected PADM research sequence")
    print(f"  baselines complete: {status['baselines_complete']}")
    print(f"  baselines missing:  {status['baselines_missing']}")
    print(
        f"  screening: {status['screening_successes']}/"
        f"{status['screening_target']} successful, "
        f"{status['screening_failures']} failed"
    )
    gate = status["screening_gate"]
    if gate is not None:
        print(f"  gate: {'PASS' if gate['passed'] else 'NO-GO'} ({gate['decision']})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("status", "next", "optimization"))
    parser.add_argument("--max-new-evaluations", type=int, default=1)
    parser.add_argument("--cores", type=int, choices=(1, 2), default=2)
    parser.add_argument("--override-screening-gate", action="store_true")
    args = parser.parse_args()
    if args.max_new_evaluations < 1:
        raise ValueError("--max-new-evaluations must be positive")

    status = _status()
    _print_status(status)
    if args.command == "status":
        return

    if status["baselines_missing"]:
        if args.command == "optimization":
            raise RuntimeError("all corrected baselines must complete before optimization")
        command = [
            sys.executable,
            str(CASE_ROOT / "run_baselines.py"),
            "--max-new-evaluations",
            str(args.max_new_evaluations),
            "--cores",
            str(args.cores),
        ]
    else:
        stage = "optimization" if args.command == "optimization" else "screening"
        command = [
            sys.executable,
            str(CASE_ROOT / "bayes_optimize_sequential.py"),
            "--stage",
            stage,
            "--max-new-evaluations",
            str(args.max_new_evaluations),
        ]
        if args.override_screening_gate:
            command.append("--override-screening-gate")
    print("[sequence] launching one strictly sequential bounded action")
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
