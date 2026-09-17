#!/usr/bin/env python3
"""Validate reconstructed inlet/outlet flux and write a machine-readable report."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import yaml


def _latest_time(case: Path) -> Path:
    times = []
    for child in case.iterdir():
        if not child.is_dir():
            continue
        try:
            value = float(child.name)
        except ValueError:
            continue
        if value > 0.0:
            times.append((value, child))
    if not times:
        raise ValueError(f"no reconstructed positive time directory found in {case}")
    return max(times)[1]


def _balanced_block(text: str, name: str) -> str:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*\{{", text)
    if not match:
        raise ValueError(f"patch {name!r} not found")
    start = text.find("{", match.start()) + 1
    depth = 1
    for index in range(start, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start:index]
    raise ValueError(f"unterminated patch block for {name!r}")


def _patch_flux(phi_text: str, patch_name: str) -> tuple[float, int]:
    block = _balanced_block(phi_text, patch_name)
    nonuniform = re.search(
        r"value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        block,
        flags=re.DOTALL,
    )
    if nonuniform:
        declared = int(nonuniform.group(1))
        values = [float(value) for value in nonuniform.group(2).split()]
        if len(values) != declared:
            raise ValueError(
                f"patch {patch_name} declares {declared} phi values but has {len(values)}"
            )
        return sum(values), declared
    uniform = re.search(r"value\s+uniform\s+([-+0-9.eE]+)\s*;", block)
    if uniform:
        raise ValueError(
            f"uniform phi on {patch_name} cannot be integrated without a face count"
        )
    raise ValueError(f"could not parse phi values for patch {patch_name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--research-config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    config = yaml.safe_load(args.research_config.read_text())
    operating = config["operating_point"]
    validation = config["validation"]
    latest = _latest_time(args.case)
    phi_text = (latest / "phi").read_text()
    inlet_flux, inlet_faces = _patch_flux(phi_text, "inlet")
    outlet_flux, outlet_faces = _patch_flux(phi_text, "outlet")

    area = float(manifest["patches"]["inlet"]["area_m2"])
    velocity = float(operating["mean_velocity_m_s"])
    expected = velocity * area
    inlet_magnitude_error = abs(abs(inlet_flux) - expected) / expected
    outlet_magnitude_error = abs(outlet_flux - expected) / expected
    imbalance = abs(inlet_flux + outlet_flux) / expected
    flow_tol = float(validation["flow_rate_relative_tolerance"])
    balance_tol = float(validation["mass_balance_relative_tolerance"])

    failures = []
    if inlet_magnitude_error > flow_tol:
        failures.append(
            f"inlet |Q| differs from U_mean*A by {inlet_magnitude_error:.3%}"
        )
    if outlet_magnitude_error > flow_tol:
        failures.append(
            f"outlet Q differs from U_mean*A by {outlet_magnitude_error:.3%}"
        )
    if imbalance > balance_tol:
        failures.append(f"relative mass imbalance is {imbalance:.3e}")

    report = {
        "schema_version": 1,
        "passed": not failures,
        "latest_time": latest.name,
        "expected_flow_rate_m3_s": expected,
        "inlet_flow_rate_m3_s": inlet_flux,
        "outlet_flow_rate_m3_s": outlet_flux,
        "inlet_face_count": inlet_faces,
        "outlet_face_count": outlet_faces,
        "inlet_relative_error": inlet_magnitude_error,
        "outlet_relative_error": outlet_magnitude_error,
        "mass_balance_relative_error": imbalance,
        "failures": failures,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if failures:
        raise RuntimeError("flow validation failed:\n  " + "\n  ".join(failures))
    print(
        f"[flow-validation] PASS: Q={outlet_flux:.8e} m3/s, "
        f"mass imbalance={imbalance:.3e}"
    )


if __name__ == "__main__":
    main()
