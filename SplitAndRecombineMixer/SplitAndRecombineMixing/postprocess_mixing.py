#!/usr/bin/env python3
"""
Compute mixing quality for the SAR mixer scalar transport case.

Reads:
  postProcessing/TAvgOutlet/*/surfaceFieldValue.dat  – area-average T at outlet
    (written by the surfaceFieldValue function object in controlDict)
  ../SplitAndRecombineHydro/sar_mixer_cad.yaml       – geometry parameters

Writes:
  mixing.csv  with columns:
    <all YAML geometry params>  |  T_avg_outlet  |  mixing_quality
    mixing_quality = |0.5 - T_avg_outlet|  → 0 means perfectly mixed
"""

import csv
import pathlib

import yaml

# ---------------------------------------------------------------------------
CASE_DIR  = pathlib.Path(__file__).parent
YAML_PATH = CASE_DIR / ".." / "SplitAndRecombineHydro" / "sar_mixer_cad.yaml"


def read_last_patch_average(patch_name: str) -> float:
    """Return the final areaAverage(T) value from a surfaceFieldValue dat file."""
    dat_files = sorted(
        (CASE_DIR / "postProcessing" / patch_name).glob("*/surfaceFieldValue.dat")
    )
    if not dat_files:
        raise FileNotFoundError(
            f"No surfaceFieldValue.dat found under "
            f"postProcessing/{patch_name}/. "
            f"Make sure the surfaceFieldValue function objects are active in "
            f"controlDict and scalarTransportFoam has completed."
        )
    dat_path = dat_files[-1]
    last_value = None
    with open(dat_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                last_value = float(parts[1])   # column 1: areaAverage(T)
    if last_value is None:
        raise ValueError(
            f"{dat_path} contains no data rows. "
            f"Did scalarTransportFoam run with the function objects active?"
        )
    return last_value


# ---------------------------------------------------------------------------
# Load geometry parameters
# ---------------------------------------------------------------------------
with open(YAML_PATH) as f:
    params = yaml.safe_load(f)

# ---------------------------------------------------------------------------
# Outlet T average from function object output
# ---------------------------------------------------------------------------
T_avg_outlet = read_last_patch_average("TAvgOutlet")

# mixing_quality → 0 means perfectly mixed (T_avg = 0.5)
mixing_quality = abs(0.5 - T_avg_outlet)

print(f"T_avg_outlet  = {T_avg_outlet:.6f}")
print(f"mixing_quality = |0.5 - {T_avg_outlet:.6f}| = {mixing_quality:.6f}")
print(f"  (0 = perfectly mixed, 0.5 = completely unmixed)")

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------
row = dict(params)
row["T_avg_outlet"]  = T_avg_outlet
row["mixing_quality"] = mixing_quality

out_path = CASE_DIR / "mixing.csv"
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(row.keys()))
    writer.writeheader()
    writer.writerow(row)

print(f"Written: {out_path}")
