#!/usr/bin/env python3
"""Prepare and run one resource-limited M10-inspired OpenFOAM evaluation.

This runner is deliberately single-case. A BO controller invokes it once per
candidate and fidelity, preserving q=1 while OpenFOAM may use up to four MPI
ranks internally.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import yaml


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent
FLOW_TEMPLATE = ROOT / "FlowCase"
SCALAR_TEMPLATE = ROOT / "ScalarTransportCase"
RESEARCH_CONFIG = ROOT / "research_config.yaml"
FIDELITY_CONFIG = ROOT / "bayes_optimize_multifidelity.yaml"
CAD_SCRIPT_NAME = "two_layer_serpentine_crossing_cad.py"
CAD_CONFIG_NAME = "two_layer_serpentine_crossing_cad.yaml"


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def safe_results_path(path: Path) -> Path:
    resolved = path.resolve()
    results_root = (ROOT / "results").resolve()
    try:
        relative = resolved.relative_to(results_root)
    except ValueError as exc:
        raise ValueError(f"result directory must stay below {results_root}") from exc
    if not relative.parts:
        raise ValueError("refusing to use the complete results root as one case")
    return resolved


def replace_token(path: Path, token: str, value: str) -> None:
    content = path.read_text(encoding="utf-8")
    if token not in content:
        raise ValueError(f"token {token!r} not found in {path}")
    path.write_text(content.replace(token, value), encoding="utf-8")


def materialize_snappy_background(
    flow: Path, manifest: dict, cell_size: float, meshing_ranks: int = 1
) -> dict:
    """Materialize a uniform blockMesh background around the CAD surface."""

    patch_bounds = [patch["bounds_m"] for patch in manifest["patches"].values()]
    geometry_min = [
        min(float(bounds[f"{axis}min"]) for bounds in patch_bounds)
        for axis in "xyz"
    ]
    geometry_max = [
        max(float(bounds[f"{axis}max"]) for bounds in patch_bounds)
        for axis in "xyz"
    ]
    margin = 2.0 * cell_size
    lower = [
        math.floor((value - margin) / cell_size) * cell_size
        for value in geometry_min
    ]
    upper = [
        math.ceil((value + margin) / cell_size) * cell_size
        for value in geometry_max
    ]
    counts = [int(round((high - low) / cell_size)) for low, high in zip(lower, upper)]
    if any(count < 5 for count in counts):
        raise RuntimeError(f"invalid snappyHexMesh background counts: {counts}")

    location = [
        float(value)
        for value in manifest["derived"]["snappy_location_in_mesh_m"]
    ]
    for index, (coordinate, low) in enumerate(zip(location, lower)):
        cell_coordinate = (coordinate - low) / cell_size
        distance_to_plane = abs(cell_coordinate - round(cell_coordinate))
        if distance_to_plane < 1.0e-6:
            location[index] += 0.071 * cell_size

    block_mesh = flow / "system" / "blockMeshDict"
    for token, value in {
        "__XMIN__": lower[0],
        "__XMAX__": upper[0],
        "__YMIN__": lower[1],
        "__YMAX__": upper[1],
        "__ZMIN__": lower[2],
        "__ZMAX__": upper[2],
        "__NX__": counts[0],
        "__NY__": counts[1],
        "__NZ__": counts[2],
    }.items():
        replace_token(block_mesh, token, f"{value:.12g}" if isinstance(value, float) else str(value))

    snappy = flow / "system" / "snappyHexMeshDict"
    for token, value in zip(
        ("__LOCATION_X__", "__LOCATION_Y__", "__LOCATION_Z__"), location
    ):
        replace_token(snappy, token, f"{value:.12g}")

    return {
        "generator": "snappyHexMesh",
        "meshing_ranks": meshing_ranks,
        "background_cell_size_m": cell_size,
        "background_bounds_m": {"minimum": lower, "maximum": upper},
        "background_cell_counts": counts,
        "location_in_mesh_m": location,
        "surface_refinement_level": 0,
        "boundary_layers": 0,
    }


def run_command(command: list[str], cwd: Path, log_name: str, env: dict[str, str]) -> None:
    log_path = cwd / log_name
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-30:])
        raise RuntimeError(
            f"{' '.join(command)} failed with exit code {completed.returncode}; "
            f"tail of {log_path}:\n{tail}"
        )


def remove_generated_processor_directories(case: Path) -> None:
    """Remove only OpenFOAM processor directories generated in this case."""

    for path in case.iterdir():
        if path.is_dir() and re.fullmatch(r"processor\d+", path.name):
            shutil.rmtree(path)


def latest_time(case: Path) -> Path:
    candidates: list[tuple[float, Path]] = []
    for child in case.iterdir():
        if child.is_dir():
            try:
                value = float(child.name)
            except ValueError:
                continue
            if value > 0.0:
                candidates.append((value, child))
    if not candidates:
        raise RuntimeError(f"no positive time directory in {case}")
    return max(candidates)[1]


def last_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"no data rows in {path}")
    return rows[-1]


def fully_developed_rectangular_duct_pressure_drop_pa(
    dynamic_viscosity_pa_s: float,
    length_m: float,
    width_m: float,
    height_m: float,
    mean_velocity_m_s: float,
    odd_terms: int = 50,
) -> float:
    """Return the exact-series laminar pressure drop for a rectangular duct.

    The design-matched straight reference has the M10 outlet cross-section
    ``w x D``, the same total axial length, and the same mean velocity.  The
    series is symmetric with respect to swapping width and height; ``a`` is
    chosen as the longer side and ``b`` as the shorter side.
    """

    values = {
        "dynamic_viscosity_pa_s": dynamic_viscosity_pa_s,
        "length_m": length_m,
        "width_m": width_m,
        "height_m": height_m,
        "mean_velocity_m_s": mean_velocity_m_s,
    }
    if any(not math.isfinite(value) or value <= 0.0 for value in values.values()):
        raise ValueError(f"straight-duct inputs must be finite and positive: {values}")
    if odd_terms < 1:
        raise ValueError("odd_terms must be positive")

    a = max(width_m, height_m)
    b = min(width_m, height_m)
    series = sum(
        math.tanh(n * math.pi * a / (2.0 * b)) / n**5
        for n in range(1, 2 * odd_terms, 2)
    )
    correction = 1.0 - 192.0 * b * series / (math.pi**5 * a)
    if correction <= 0.0:
        raise RuntimeError(f"invalid rectangular-duct correction {correction}")
    return (
        12.0
        * dynamic_viscosity_pa_s
        * length_m
        * mean_velocity_m_s
        / (b**2 * correction)
    )


def validate_scalar_history(
    path: Path,
    window: int = 50,
    intensity_absolute_tolerance: float = 1.0e-8,
    intensity_relative_tolerance: float = 1.0e-2,
    mean_tolerance: float = 1.0e-6,
) -> dict:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < window:
        raise RuntimeError(f"scalar history has {len(rows)} rows; at least {window} required")
    tail = rows[-window:]
    intensity = [float(row["flux_weighted_intensity_of_segregation"]) for row in tail]
    mean = [float(row["flux_weighted_mean_concentration"]) for row in tail]
    if not all(math.isfinite(value) for value in intensity + mean):
        raise RuntimeError("scalar history contains non-finite outlet statistics")
    if min(intensity) < -1.0e-10 or max(intensity) > 1.0 + 1.0e-10:
        raise RuntimeError("flux-weighted segregation intensity lies outside [0, 1]")
    intensity_span = max(intensity) - min(intensity)
    mean_span = max(mean) - min(mean)
    intensity_tolerance = max(
        intensity_absolute_tolerance,
        intensity_relative_tolerance * max(abs(intensity[-1]), 1.0e-12),
    )
    if intensity_span > intensity_tolerance or mean_span > mean_tolerance:
        raise RuntimeError(
            "scalar outlet statistics are not stable: "
            f"intensity span={intensity_span:.3e}, mean span={mean_span:.3e}"
        )
    return {
        "history_rows": len(rows),
        "window": window,
        "intensity_span": intensity_span,
        "intensity_tolerance": intensity_tolerance,
        "intensity_relative_tolerance": intensity_relative_tolerance,
        "mean_span": mean_span,
        "mean_tolerance": mean_tolerance,
    }


def validate_scalar_bounds(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(
        r"min\(T\)\s*=\s*([0-9.eE+-]+).*?max\(T\)\s*=\s*([0-9.eE+-]+)",
        text,
        flags=re.DOTALL,
    )
    if not matches:
        raise RuntimeError(f"scalar bounds not found in {path}")
    minimum, maximum = map(float, matches[-1])
    undershoot = max(0.0, -minimum)
    overshoot = max(0.0, maximum - 1.0)
    tolerance = float(
        load_yaml(RESEARCH_CONFIG)["validation_gates"][
            "scalar_bound_excursion_tolerance"
        ]
    )
    if max(undershoot, overshoot) > tolerance:
        raise RuntimeError(
            "unbounded scalar excursion exceeds tolerance: "
            f"min={minimum:.8g}, max={maximum:.8g}, tolerance={tolerance:.3g}"
        )
    return {
        "minimum_T": minimum,
        "maximum_T": maximum,
        "undershoot": undershoot,
        "overshoot": overshoot,
        "excursion_tolerance": tolerance,
    }


def validate_mesh_quality(path: Path, fidelity: str) -> dict:
    """Validate hex-dominant composition and OpenFOAM quality criteria.

    A body-fitted snappyHexMesh grid necessarily contains a thin layer of
    prism/polyhedral cut cells at oblique walls.  It must nevertheless contain
    no tetrahedra, remain predominantly hexahedral, and pass every explicit
    operator-oriented criterion in ``meshQualityDict``.  The default
    ``checkMesh`` concavity/warpage diagnostics are retained as separately
    bounded populations instead of being silently ignored.
    """

    text = path.read_text(encoding="utf-8", errors="replace")

    def required_integer(pattern: str, label: str) -> int:
        match = re.search(pattern, text, flags=re.MULTILINE)
        if not match:
            raise RuntimeError(f"{label} not found in {path}")
        return int(match.group(1))

    faces = required_integer(r"^\s*faces:\s+(\d+)", "face count")
    cells = required_integer(r"^\s*cells:\s+(\d+)", "cell count")
    hexahedra = required_integer(r"^\s*hexahedra:\s+(\d+)", "hexahedron count")
    prisms = required_integer(r"^\s*prisms:\s+(\d+)", "prism count")
    tet_wedges = required_integer(r"^\s*tet wedges:\s+(\d+)", "tet-wedge count")
    tetrahedra = required_integer(r"^\s*tetrahedra:\s+(\d+)", "tetrahedron count")
    polyhedra = required_integer(r"^\s*polyhedra:\s+(\d+)", "polyhedron count")
    concave_match = re.search(r"Concave cells .* number of cells:\s*(\d+)", text)
    concave_cells = int(concave_match.group(1)) if concave_match else 0
    angle_match = re.search(r"Max concave angle\s*=\s*([0-9.eE+-]+)", text)
    max_concave_angle = float(angle_match.group(1)) if angle_match else 0.0
    warped_match = re.search(
        r"There are\s+(\d+)\s+faces with ratio between projected and actual area",
        text,
    )
    warped_faces = int(warped_match.group(1)) if warped_match else 0
    flatness_match = re.search(
        r"Face flatness .*?min\s*=\s*([0-9.eE+-]+)", text
    )
    minimum_face_flatness = float(flatness_match.group(1)) if flatness_match else 1.0
    failed_match = re.search(r"Failed\s+(\d+)\s+mesh checks", text)
    default_failed_checks = int(failed_match.group(1)) if failed_match else 0

    marker = "Checking faces in error :"
    if marker not in text:
        raise RuntimeError("checkMesh did not evaluate meshQualityDict")
    error_block = text.split(marker, 1)[1].split("\n\n", 1)[0]
    explicit_error_counts = [
        int(value)
        for value in re.findall(r":\s*(\d+)\s*$", error_block, flags=re.MULTILINE)
    ]
    if len(explicit_error_counts) < 9:
        raise RuntimeError("incomplete meshQualityDict error report")
    if any(explicit_error_counts):
        raise RuntimeError(
            "mesh violates explicit transport-oriented criteria: "
            f"error counts={explicit_error_counts}"
        )

    policy = load_yaml(RESEARCH_CONFIG)["validation_gates"]["mesh_quality"][fidelity]
    hex_fraction = hexahedra / cells
    tetra_fraction = tetrahedra / cells
    tet_wedge_fraction = tet_wedges / cells
    concave_fraction = concave_cells / cells
    warped_face_fraction = warped_faces / faces
    limits = {
        "minimum_hexahedral_cell_fraction": float(policy["minimum_hexahedral_cell_fraction"]),
        "maximum_tetrahedral_cell_fraction": float(policy["maximum_tetrahedral_cell_fraction"]),
        "maximum_tet_wedge_cell_fraction": float(policy["maximum_tet_wedge_cell_fraction"]),
        "maximum_concave_cell_fraction": float(policy["maximum_concave_cell_fraction"]),
        "maximum_concave_angle_degrees": float(policy["maximum_concave_angle_degrees"]),
        "maximum_warped_face_fraction": float(policy["maximum_warped_face_fraction"]),
        "minimum_face_flatness": float(policy["minimum_face_flatness"]),
        "maximum_default_failed_checks": int(policy["maximum_default_failed_checks"]),
    }
    failures: list[str] = []
    if hex_fraction < limits["minimum_hexahedral_cell_fraction"]:
        failures.append(f"hexahedral fraction {hex_fraction:.6g}")
    if tetra_fraction > limits["maximum_tetrahedral_cell_fraction"]:
        failures.append(f"tetrahedral fraction {tetra_fraction:.6g}")
    if tet_wedge_fraction > limits["maximum_tet_wedge_cell_fraction"]:
        failures.append(f"tet-wedge fraction {tet_wedge_fraction:.6g}")
    if concave_fraction > limits["maximum_concave_cell_fraction"]:
        failures.append(f"concave fraction {concave_fraction:.6g}")
    if max_concave_angle > limits["maximum_concave_angle_degrees"]:
        failures.append(f"maximum concave angle {max_concave_angle:.6g}")
    if warped_face_fraction > limits["maximum_warped_face_fraction"]:
        failures.append(f"warped-face fraction {warped_face_fraction:.6g}")
    if minimum_face_flatness < limits["minimum_face_flatness"]:
        failures.append(f"minimum face flatness {minimum_face_flatness:.6g}")
    if default_failed_checks > limits["maximum_default_failed_checks"]:
        failures.append(f"default failed checks {default_failed_checks}")
    if failures:
        raise RuntimeError("mesh composition gate failed: " + ", ".join(failures))

    return {
        "faces": faces,
        "cells": cells,
        "hexahedra": hexahedra,
        "prisms": prisms,
        "tet_wedges": tet_wedges,
        "tetrahedra": tetrahedra,
        "polyhedra": polyhedra,
        "hexahedral_cell_fraction": hex_fraction,
        "tetrahedral_cell_fraction": tetra_fraction,
        "tet_wedge_cell_fraction": tet_wedge_fraction,
        "concave_cells": concave_cells,
        "concave_cell_fraction": concave_fraction,
        "maximum_concave_angle_degrees": max_concave_angle,
        "warped_faces": warped_faces,
        "warped_face_fraction": warped_face_fraction,
        "minimum_face_flatness": minimum_face_flatness,
        "default_failed_checks": default_failed_checks,
        "explicit_quality_error_counts": explicit_error_counts,
        "limits": limits,
    }


def validate_boundary_topology(path: Path) -> dict:
    """Require the three solver-facing patches and no leaked background."""

    text = path.read_text(encoding="utf-8", errors="replace")
    patch_faces: dict[str, int] = {}
    for match in re.finditer(
        r"(?m)^\s*([A-Za-z0-9_.:-]+)\s*\n\s*\{([^{}]*)\}", text
    ):
        faces = re.search(r"\bnFaces\s+(\d+)\s*;", match.group(2))
        if faces:
            patch_faces[match.group(1)] = int(faces.group(1))
    required = ("inlet", "outlet", "walls")
    missing = [name for name in required if patch_faces.get(name, 0) <= 0]
    if missing:
        raise RuntimeError(f"missing or empty mesh boundary patches: {missing}")
    leaked = {
        name: count
        for name, count in patch_faces.items()
        if name not in required and count > 0
    }
    if leaked:
        raise RuntimeError(f"unexpected non-empty mesh boundary patches: {leaked}")
    return {"patch_face_counts": patch_faces}


def reuse_validated_mesh(
    source_result: Path, target_flow: Path, target_metadata: dict
) -> str:
    source_result = safe_results_path(source_result)
    source_flow = source_result / "FlowCase"
    source_metadata_path = source_result / "case_metadata.json"
    source_mesh = source_flow / "constant" / "polyMesh"
    source_validation = source_result / "mesh_validation.json"
    for required in (source_metadata_path, source_mesh, source_validation):
        if not required.exists():
            raise FileNotFoundError(f"validated mesh source is incomplete: {required}")

    source_metadata = json.loads(source_metadata_path.read_text(encoding="utf-8"))
    if source_metadata["fidelity"] != target_metadata["fidelity"]:
        raise ValueError("mesh reuse requires the same fidelity")
    if not math.isclose(
        float(source_metadata["cell_size_m"]),
        float(target_metadata["cell_size_m"]),
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise ValueError("mesh reuse requires the same nominal cell size")
    if source_metadata["geometry_manifest"] != target_metadata["geometry_manifest"]:
        raise ValueError("mesh reuse requires an identical geometry manifest")

    target_mesh = target_flow / "constant" / "polyMesh"
    shutil.copytree(source_mesh, target_mesh)
    relative_source = str(source_result.relative_to(ROOT))
    return relative_source


def _balanced_patch_block(text: str, patch_name: str) -> str:
    match = re.search(rf"(?m)^\s*{re.escape(patch_name)}\s*\{{", text)
    if not match:
        raise RuntimeError(f"patch {patch_name!r} not found in phi")
    start = text.find("{", match.start()) + 1
    depth = 1
    for index in range(start, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start:index]
    raise RuntimeError(f"unterminated phi block for patch {patch_name!r}")


def _integrated_patch_flux(phi_text: str, patch_name: str) -> float:
    block = _balanced_patch_block(phi_text, patch_name)
    match = re.search(
        r"value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        block,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"nonuniform flux values not found for patch {patch_name!r}")
    count = int(match.group(1))
    values = [float(value) for value in match.group(2).split()]
    if len(values) != count:
        raise RuntimeError(f"patch {patch_name!r} declares {count} fluxes but contains {len(values)}")
    return sum(values)


def validate_flow_balance(flow: Path, metadata: dict) -> dict:
    phi_text = (latest_time(flow) / "phi").read_text(encoding="utf-8")
    inlet_flux = _integrated_patch_flux(phi_text, "inlet")
    outlet_flux = _integrated_patch_flux(phi_text, "outlet")
    patches = metadata["geometry_manifest"]["patches"]
    inlet_area = float(patches["inlet1"]["area_m2"]) + float(patches["inlet2"]["area_m2"])
    expected = float(metadata["inlet_mean_velocity_m_s"]) * inlet_area
    inlet_error = abs(abs(inlet_flux) - expected) / expected
    outlet_error = abs(outlet_flux - expected) / expected
    imbalance = abs(inlet_flux + outlet_flux) / expected
    gates = load_yaml(RESEARCH_CONFIG)["validation_gates"]
    flow_tolerance = float(
        gates["cad_to_mesh_inlet_flow_relative_tolerance"][metadata["fidelity"]]
    )
    balance_tolerance = float(gates["mass_balance_relative_tolerance"])
    if inlet_error > flow_tolerance or outlet_error > flow_tolerance or imbalance > balance_tolerance:
        raise RuntimeError(
            "flow balance failed: "
            f"inlet error={inlet_error:.3e}, outlet error={outlet_error:.3e}, "
            f"imbalance={imbalance:.3e}"
        )
    return {
        "expected_flow_rate_m3_s": expected,
        "inlet_flow_rate_m3_s": inlet_flux,
        "outlet_flow_rate_m3_s": outlet_flux,
        "inlet_relative_error": inlet_error,
        "outlet_relative_error": outlet_error,
        "mass_balance_relative_error": imbalance,
        "cad_to_mesh_flow_relative_tolerance": flow_tolerance,
        "mass_balance_relative_tolerance": balance_tolerance,
    }


def prepare_geometry_config(
    protocol: str,
    source: Path,
    destination: Path,
    inlet_lead_mm: float | None = None,
    outlet_lead_mm: float | None = None,
) -> dict:
    config = load_yaml(source)
    if protocol == "review":
        config["geometry"]["number_of_units"] = 6
        # Six-unit core is 4.89 mm. Symmetric 0.08 mm leads give the review's
        # standardized 5.05 mm axial length.
        config["geometry"]["inlet_lead_length"] = 0.08
        config["geometry"]["outlet_lead_length"] = 0.08
    if (inlet_lead_mm is None) != (outlet_lead_mm is None):
        raise ValueError("inlet and outlet lead overrides must be supplied together")
    if inlet_lead_mm is not None and outlet_lead_mm is not None:
        if inlet_lead_mm <= 0.0 or outlet_lead_mm <= 0.0:
            raise ValueError("lead overrides must be positive")
        config["geometry"]["inlet_lead_length"] = inlet_lead_mm
        config["geometry"]["outlet_lead_length"] = outlet_lead_mm
    with destination.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return config


def prepare_case(args: argparse.Namespace) -> tuple[Path, Path, dict, dict[str, str]]:
    result_dir = safe_results_path(args.results_dir)
    if result_dir.exists():
        if not args.force:
            raise FileExistsError(f"{result_dir} exists; pass --force to replace this case")
        shutil.rmtree(result_dir)
    flow = result_dir / "FlowCase"
    scalar = result_dir / "ScalarTransportCase"
    shutil.copytree(FLOW_TEMPLATE, flow)
    shutil.copytree(SCALAR_TEMPLATE, scalar)

    geometry_source = args.geometry_config.resolve() if args.geometry_config else FLOW_TEMPLATE / CAD_CONFIG_NAME
    geometry_config = flow / CAD_CONFIG_NAME
    prepare_geometry_config(
        args.protocol,
        geometry_source,
        geometry_config,
        args.inlet_lead_mm,
        args.outlet_lead_mm,
    )
    tri_surface = flow / "constant" / "triSurface"
    tri_surface.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    wm_options = env.get("WM_OPTIONS")
    if wm_options:
        local_lib = REPO_ROOT / "platforms" / wm_options / "lib"
        env["FOAM_USER_LIBBIN"] = str(local_lib)
        env["LD_LIBRARY_PATH"] = f"{local_lib}:{env.get('LD_LIBRARY_PATH', '')}"

    run_command(
        [sys.executable, CAD_SCRIPT_NAME, "--config", str(geometry_config), "--output-dir", str(tri_surface)],
        flow,
        "log.cad",
        env,
    )
    manifest = json.loads((tri_surface / "geometry_manifest.json").read_text(encoding="utf-8"))

    research = load_yaml(RESEARCH_CONFIG)
    fidelity_cfg = load_yaml(FIDELITY_CONFIG)["fidelities"][args.fidelity]
    fluid = research["literature_reproduction"]["fluid"]
    rho = float(fluid["density_kg_m3"])
    mu = float(fluid["dynamic_viscosity_Pa_s"])
    nu = mu / rho
    diffusivity = (
        float(research["review_matched_benchmark"]["scalar_diffusivity_m2_s"])
        if args.protocol == "review"
        else float(fluid["scalar_diffusivity_m2_s"])
    )
    hydraulic_diameter = float(manifest["dimensions_m"]["inlet_hydraulic_diameter"])
    inlet_speed = float(args.reynolds) * nu / hydraulic_diameter
    cell_size = float(
        args.cell_size_m
        if args.cell_size_m is not None
        else fidelity_cfg["nominal_cell_size_m"]
    )

    replace_token(flow / "0" / "U", "__INLET_SPEED__", f"{inlet_speed:.12g}")
    replace_token(flow / "constant" / "transportProperties", "__KINEMATIC_VISCOSITY__", f"{nu:.12g}")
    mesh = materialize_snappy_background(flow, manifest, cell_size, args.cores)
    replace_token(flow / "system" / "decomposeParDict", "__CORES__", str(args.cores))
    replace_token(scalar / "system" / "decomposeParDict", "__CORES__", str(args.cores))
    replace_token(scalar / "constant" / "transportProperties", "__SCALAR_DIFFUSIVITY__", f"{diffusivity:.12g}")

    metadata = {
        "schema_version": 1,
        "scientific_framing": load_yaml(FIDELITY_CONFIG)["scientific_framing"],
        "protocol": args.protocol,
        "fidelity": args.fidelity,
        "fidelity_coordinate": float(fidelity_cfg["coordinate"]),
        "cell_size_m": cell_size,
        "cell_size_source": "command_line_override" if args.cell_size_m is not None else "fidelity_config",
        "mesh": mesh,
        "reynolds_number": float(args.reynolds),
        "density_kg_m3": rho,
        "dynamic_viscosity_Pa_s": mu,
        "kinematic_viscosity_m2_s": nu,
        "scalar_diffusivity_m2_s": diffusivity,
        "scalar_transport_numerics": research["scalar_transport_numerics"],
        "straight_pressure_reference": research["straight_pressure_reference"],
        "inlet_hydraulic_diameter_m": hydraulic_diameter,
        "inlet_mean_velocity_m_s": inlet_speed,
        "cores": int(args.cores),
        "openfoam": {
            "version": env.get("WM_PROJECT_VERSION"),
            "wm_options": env.get("WM_OPTIONS"),
        },
        "geometry_manifest": manifest,
    }
    (result_dir / "case_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return flow, scalar, metadata, env


def run_case(args: argparse.Namespace) -> None:
    if not 1 <= args.cores <= 4:
        raise ValueError("cores must stay in [1, 4]")
    flow, scalar, metadata, env = prepare_case(args)
    result_dir = flow.parent
    print(f"prepared {result_dir}")
    if args.prepare_only:
        return

    required = [
        "checkMesh",
        "decomposePar",
        "simpleFoam",
        "reconstructPar",
        "scalarTransportFoam",
        "postProcess",
    ]
    if args.reuse_mesh_from is None:
        required.extend(
            [
                "surfaceCheck",
                "blockMesh",
                "surfaceFeatureExtract",
                "snappyHexMesh",
                "createPatch",
            ]
        )
        if args.cores > 1:
            required.extend(["mpirun", "reconstructParMesh"])
    missing = [command for command in required if shutil.which(command, path=env.get("PATH")) is None]
    if missing:
        raise RuntimeError(f"OpenFOAM environment is not sourced; missing commands: {', '.join(missing)}")
    required_version = str(load_yaml(RESEARCH_CONFIG)["software"]["openfoam_version"])
    if env.get("WM_PROJECT_VERSION") != required_version:
        raise RuntimeError(
            f"this study is qualified for OpenFOAM-{required_version}; "
            "source its etc/bashrc first "
            f"(found {env.get('WM_PROJECT_VERSION')!r})"
        )

    if args.reuse_mesh_from is not None:
        mesh_source = reuse_validated_mesh(
            args.reuse_mesh_from, flow, metadata
        )
        metadata["mesh_reused_from"] = mesh_source
        (result_dir / "case_metadata.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        print(f"reused validated mesh from {mesh_source}")
    else:
        surface_file = "constant/triSurface/two_layer_serpentine_crossing_mixer.stl"
        run_command(
            ["surfaceCheck", "-checkSelfIntersection", surface_file],
            flow,
            "log.surfaceCheck",
            env,
        )
        run_command(["blockMesh"], flow, "log.blockMesh", env)
        run_command(
            ["surfaceFeatureExtract"], flow, "log.surfaceFeatureExtract", env
        )
        if args.cores == 1:
            run_command(
                ["snappyHexMesh", "-overwrite"], flow, "log.snappyHexMesh", env
            )
        else:
            run_command(
                ["decomposePar", "-force", "-no-fields"],
                flow,
                "log.decomposeParMesh",
                env,
            )
            run_command(
                [
                    "mpirun",
                    "-np",
                    str(args.cores),
                    "snappyHexMesh",
                    "-overwrite",
                    "-parallel",
                ],
                flow,
                "log.snappyHexMesh",
                env,
            )
            run_command(
                ["reconstructParMesh", "-constant", "-mergeTol", "1e-6"],
                flow,
                "log.reconstructParMesh",
                env,
            )
            remove_generated_processor_directories(flow)
        run_command(["createPatch", "-overwrite"], flow, "log.createPatch", env)
    # Recheck both fresh and reused meshes under this checkout's fixed policy.
    run_command(
        ["checkMesh", "-allGeometry", "-allTopology", "-meshQuality"],
        flow,
        "log.checkMesh",
        env,
    )
    mesh_validation = validate_mesh_quality(
        flow / "log.checkMesh", metadata["fidelity"]
    )
    mesh_validation["boundary_topology"] = validate_boundary_topology(
        flow / "constant" / "polyMesh" / "boundary"
    )
    (result_dir / "mesh_validation.json").write_text(
        json.dumps(mesh_validation, indent=2) + "\n", encoding="utf-8"
    )
    print(f"mesh complete at {metadata['cell_size_m'] * 1.0e6:.1f} um")
    if args.mesh_only:
        return

    run_command(["decomposePar", "-force"], flow, "log.decomposePar", env)
    flow_solver = ["simpleFoam"] if args.cores == 1 else ["mpirun", "-np", str(args.cores), "simpleFoam", "-parallel"]
    run_command(flow_solver, flow, "log.simpleFoam", env)
    run_command(["reconstructPar", "-latestTime"], flow, "log.reconstructPar", env)

    flow_time = latest_time(flow)
    shutil.copytree(flow / "constant" / "polyMesh", scalar / "constant" / "polyMesh")
    shutil.copy2(flow_time / "U", scalar / "0" / "U")
    shutil.copy2(flow_time / "phi", scalar / "0" / "phi")
    run_command(["decomposePar", "-force"], scalar, "log.decomposePar", env)
    scalar_solver = ["scalarTransportFoam"] if args.cores == 1 else ["mpirun", "-np", str(args.cores), "scalarTransportFoam", "-parallel"]
    run_command(scalar_solver, scalar, "log.scalarTransportFoam", env)
    run_command(["reconstructPar", "-latestTime"], scalar, "log.reconstructPar", env)
    run_command(
        ["postProcess", "-latestTime", "-field", "T", "-func", "fieldMinMax(T)"],
        scalar,
        "log.scalarBounds",
        env,
    )

    flow_validation = validate_flow_balance(flow, metadata)
    scalar_validation = validate_scalar_history(scalar / "mixing.csv")
    scalar_bounds = validate_scalar_bounds(scalar / "log.scalarBounds")
    pressure_row = last_csv_row(flow / "pressureDrop.csv")
    mixing_row = last_csv_row(scalar / "mixing.csv")
    pressure_kinematic = float(pressure_row["pressure_drop_m2_s2"])
    pressure_pa = pressure_kinematic * float(metadata["density_kg_m3"])
    dimensions = metadata["geometry_manifest"]["dimensions_m"]
    straight_pressure_pa = fully_developed_rectangular_duct_pressure_drop_pa(
        dynamic_viscosity_pa_s=float(metadata["dynamic_viscosity_Pa_s"]),
        length_m=float(dimensions["total_axial_length"]),
        width_m=float(dimensions["w"]),
        height_m=float(dimensions["D"]),
        mean_velocity_m_s=float(metadata["inlet_mean_velocity_m_s"]),
        odd_terms=int(metadata["straight_pressure_reference"]["odd_series_terms"]),
    )
    intensity = float(mixing_row["flux_weighted_intensity_of_segregation"])
    area_intensity = float(mixing_row["intensity_of_segregation"])
    if not 0.0 <= intensity <= 1.0:
        raise RuntimeError(f"invalid flux-weighted segregation intensity {intensity}")
    if not 0.0 <= area_intensity <= 1.0:
        raise RuntimeError(f"invalid area-weighted segregation intensity {area_intensity}")
    objectives = {
        "pressure_drop_m2_s2": pressure_kinematic,
        "pressure_drop_Pa": pressure_pa,
        "straight_channel_pressure_drop_Pa": straight_pressure_pa,
        "pressure_ratio_to_straight": pressure_pa / straight_pressure_pa,
        "flux_weighted_intensity_of_segregation": intensity,
        "flux_weighted_mixing_index": 1.0 - math.sqrt(intensity),
        "flux_weighted_mean_concentration": float(mixing_row["flux_weighted_mean_concentration"]),
        "area_weighted_intensity_of_segregation": area_intensity,
        "area_weighted_mixing_index": 1.0 - math.sqrt(area_intensity),
        "area_weighted_mean_concentration": float(mixing_row["mean_concentration"]),
        "fidelity": metadata["fidelity"],
        "fidelity_coordinate": metadata["fidelity_coordinate"],
        "reynolds_number": metadata["reynolds_number"],
        "protocol": metadata["protocol"],
        "mass_balance_relative_error": flow_validation["mass_balance_relative_error"],
        "scalar_final_window_intensity_span": scalar_validation["intensity_span"],
        "scalar_final_window_mean_span": scalar_validation["mean_span"],
        "minimum_scalar": scalar_bounds["minimum_T"],
        "maximum_scalar": scalar_bounds["maximum_T"],
    }
    (result_dir / "flow_validation.json").write_text(
        json.dumps(flow_validation, indent=2) + "\n", encoding="utf-8"
    )
    (result_dir / "scalar_validation.json").write_text(
        json.dumps(scalar_validation, indent=2) + "\n", encoding="utf-8"
    )
    (result_dir / "scalar_bounds.json").write_text(
        json.dumps(scalar_bounds, indent=2) + "\n", encoding="utf-8"
    )
    (result_dir / "objectives.json").write_text(json.dumps(objectives, indent=2) + "\n", encoding="utf-8")
    with (result_dir / "objectives.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(objectives))
        writer.writeheader()
        writer.writerow(objectives)
    print(json.dumps(objectives, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", choices=("original", "review"), default="original")
    parser.add_argument("--fidelity", choices=("coarse", "fine"), required=True)
    parser.add_argument("--reynolds", type=float, required=True)
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument(
        "--cell-size-m",
        type=float,
        help="Recorded mesh-qualification override; omitted in production BO runs",
    )
    parser.add_argument(
        "--reuse-mesh-from",
        type=Path,
        help="Reuse an identical, revalidated result mesh below this study's results/ directory",
    )
    parser.add_argument("--geometry-config", type=Path)
    parser.add_argument(
        "--inlet-lead-mm",
        type=float,
        help="Geometry-audit override; requires --outlet-lead-mm",
    )
    parser.add_argument(
        "--outlet-lead-mm",
        type=float,
        help="Geometry-audit override; requires --inlet-lead-mm",
    )
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.reynolds <= 0.0:
        parser.error("--reynolds must be positive")
    if args.cell_size_m is not None and args.cell_size_m <= 0.0:
        parser.error("--cell-size-m must be positive")
    if (args.inlet_lead_mm is None) != (args.outlet_lead_mm is None):
        parser.error("--inlet-lead-mm and --outlet-lead-mm must be supplied together")
    if args.inlet_lead_mm is not None and (
        args.inlet_lead_mm <= 0.0 or args.outlet_lead_mm <= 0.0
    ):
        parser.error("lead overrides must be positive")
    return args


if __name__ == "__main__":
    run_case(parse_args())
