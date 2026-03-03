#!/usr/bin/env python3
"""
Compute pressure drop and mesh statistics for the SAR mixer case.

Reads:
  postProcessing/pAvgInlet/*/surfaceFieldValue.dat   – area-average p at inlet
  postProcessing/pAvgOutlet/*/surfaceFieldValue.dat  – area-average p at outlet
    (written by the surfaceFieldValue function objects in controlDict)
  0/V   – cell volumes written by 'postProcess -func writeCellVolumes -time 0'
  sar_mixer_cad.yaml  – geometry parameters (written as CSV columns)

Writes:
  pressureDrop.csv  with columns:
    <all YAML geometry params>  |  pressureDrop_Pa  |
    delta_x_min  delta_x_max  delta_x_mean   (if 0/V is available)
"""

import csv
import math
import pathlib

import numpy as np
import yaml

# ---------------------------------------------------------------------------
CASE_DIR  = pathlib.Path(__file__).parent
YAML_PATH = CASE_DIR / "sar_mixer_cad.yaml"


def read_last_patch_average(patch_name: str) -> float:
    """Return the final areaAverage(p) value from a surfaceFieldValue dat file."""
    # Function objects write to postProcessing/<name>/<start_time>/surfaceFieldValue.dat
    # There may be multiple start-time sub-directories if the run was restarted.
    dat_files = sorted(
        (CASE_DIR / "postProcessing" / patch_name).glob("*/surfaceFieldValue.dat")
    )
    if not dat_files:
        raise FileNotFoundError(
            f"No surfaceFieldValue.dat found under "
            f"postProcessing/{patch_name}/. "
            f"Make sure the surfaceFieldValue function objects are active in "
            f"controlDict and simpleFoam has completed."
        )
    # Use the last file (latest start time) and read its last data row
    dat_path = dat_files[-1]
    last_value = None
    with open(dat_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                last_value = float(parts[1])   # column 1: areaAverage(p)
    if last_value is None:
        raise ValueError(
            f"{dat_path} contains no data rows. "
            f"Did simpleFoam run with the function objects active?"
        )
    return last_value


# ---------------------------------------------------------------------------
# Load geometry parameters
# ---------------------------------------------------------------------------
with open(YAML_PATH) as f:
    params = yaml.safe_load(f)

# DEPTH = 0.01 * H * scale  (must match sar_mixer_cad.py)
DEPTH = params["scale"] * params["H"] * 0.01

# ---------------------------------------------------------------------------
# Pressure drop from function object output
# ---------------------------------------------------------------------------
p_inlet  = read_last_patch_average("pAvgInlet")
p_outlet = read_last_patch_average("pAvgOutlet")
delta_p  = p_inlet - p_outlet

print(f"p_inlet  = {p_inlet:.6e} Pa")
print(f"p_outlet = {p_outlet:.6e} Pa")
print(f"Pressure drop = {delta_p:.6e} Pa")

# ---------------------------------------------------------------------------
# Mesh statistics from cell volumes  (0/V written by writeCellVolumes)
#
# The mesh is one cell thick in z (DEPTH).  Characteristic xy cell size:
#   delta_x = sqrt(V_3d / DEPTH)
# ---------------------------------------------------------------------------
have_mesh_stats = False
delta_x_min = delta_x_max = delta_x_mean = math.nan

V_path = CASE_DIR / "0" / "V"
if V_path.exists():
    try:
        import foamlib
        V   = np.asarray(
            foamlib.FoamFieldFile(V_path).internal_field, dtype=float
        )
        dxs = np.sqrt(V / DEPTH)
        delta_x_min  = float(dxs.min())
        delta_x_max  = float(dxs.max())
        delta_x_mean = float(dxs.mean())
        have_mesh_stats = True
        print(f"Cell sizes: min={delta_x_min:.3e} m  "
              f"max={delta_x_max:.3e} m  mean={delta_x_mean:.3e} m")
    except Exception as exc:
        print(f"Warning: could not compute mesh statistics ({exc})")
else:
    print("Warning: 0/V not found; skipping mesh statistics. "
          "Run 'postProcess -func writeCellVolumes -time 0' to generate it.")

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------
row = dict(params)
row["pressureDrop_Pa"] = delta_p
if have_mesh_stats:
    row["delta_x_min"]  = delta_x_min
    row["delta_x_max"]  = delta_x_max
    row["delta_x_mean"] = delta_x_mean

out_path = CASE_DIR / "pressureDrop.csv"
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(row.keys()))
    writer.writeheader()
    writer.writerow(row)

print(f"Written: {out_path}")
