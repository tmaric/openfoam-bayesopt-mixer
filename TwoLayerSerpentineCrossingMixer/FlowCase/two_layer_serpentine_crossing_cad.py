#!/usr/bin/env python3
"""Build an M10-inspired two-layer serpentine-crossing fluid domain.

The construction follows the topology in Hossain et al. (2017): an inverse-N
channel in the upper layer and an N channel in the lower layer. Repeated
vertical segments join successive diagonals. The face-adjacent layers connect
where the projected channels overlap: at every X node and vertical segment.

The paper specifies H, w, b, P, d, unit count, inlet/outlet sections, and the
locations where layers interconnect. It does not publish CAD/mask data or
independent aperture, lead-transition, and corner-treatment dimensions. The
manifest therefore records this geometry as an inspired reconstruction rather
than an exact reproduction.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import cadquery as cq
import yaml


MM = 1.0e-3


def _positive(name: str, value: float) -> float:
    value = float(value)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive, received {value}")
    return value


def _segment(
    start: tuple[float, float],
    end: tuple[float, float],
    width: float,
    z0: float,
    depth: float,
) -> cq.Workplane:
    """Return a rectangular channel prism centred on an XY line segment."""
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = math.hypot(dx, dy)
    if length <= 0.0:
        raise ValueError("channel segment has zero length")
    angle_deg = math.degrees(math.atan2(dy, dx))
    return (
        cq.Workplane("XY")
        .box(length, width, depth, centered=(False, True, False))
        .rotate((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), angle_deg)
        .translate((start[0], start[1], z0))
    )


def _box(
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    z0: float,
    z1: float,
) -> cq.Workplane:
    return (
        cq.Workplane("XY")
        .box(x1 - x0, y1 - y0, z1 - z0, centered=(False, False, False))
        .translate((x0, y0, z0))
    )


def _fuse_all(solids: list[cq.Workplane]) -> cq.Workplane:
    if not solids:
        raise ValueError("cannot fuse an empty solid list")
    fused = solids[0]
    for solid in solids[1:]:
        fused = fused.union(solid)
    return fused.clean()


def _bounds(face: cq.Face) -> dict[str, float]:
    bb = face.BoundingBox()
    return {
        "xmin": bb.xmin,
        "xmax": bb.xmax,
        "ymin": bb.ymin,
        "ymax": bb.ymax,
        "zmin": bb.zmin,
        "zmax": bb.zmax,
    }


def _combined_bounds(faces: list[cq.Face]) -> dict[str, float]:
    boxes = [_bounds(face) for face in faces]
    return {
        key: (min(box[key] for box in boxes) if key.endswith("min") else max(box[key] for box in boxes))
        for key in boxes[0]
    }


def _stl_block(
    faces: list[cq.Face],
    solid_name: str,
    tolerance: float,
    angular_tolerance: float = 0.1,
) -> str:
    if not faces:
        raise RuntimeError(f"no faces classified for patch {solid_name!r}")
    compound = cq.Compound.makeCompound(faces)
    vertices, triangles = compound.tessellate(tolerance, angular_tolerance)
    lines = [f"solid {solid_name}"]
    for triangle in triangles:
        v0, v1, v2 = (vertices[index] for index in triangle)
        e1 = (v1.x - v0.x, v1.y - v0.y, v1.z - v0.z)
        e2 = (v2.x - v0.x, v2.y - v0.y, v2.z - v0.z)
        nx = e1[1] * e2[2] - e1[2] * e2[1]
        ny = e1[2] * e2[0] - e1[0] * e2[2]
        nz = e1[0] * e2[1] - e1[1] * e2[0]
        magnitude = math.sqrt(nx * nx + ny * ny + nz * nz)
        if magnitude > 1.0e-30:
            nx, ny, nz = nx / magnitude, ny / magnitude, nz / magnitude
        lines.extend(
            [
                f"  facet normal {nx:.9e} {ny:.9e} {nz:.9e}",
                "    outer loop",
                f"      vertex {v0.x:.9e} {v0.y:.9e} {v0.z:.9e}",
                f"      vertex {v1.x:.9e} {v1.y:.9e} {v1.z:.9e}",
                f"      vertex {v2.x:.9e} {v2.y:.9e} {v2.z:.9e}",
                "    endloop",
                "  endfacet",
            ]
        )
    lines.append(f"endsolid {solid_name}")
    return "\n".join(lines) + "\n"


def build_geometry(config: dict) -> tuple[cq.Workplane, dict]:
    source = config["geometry"]
    H = _positive("main_span_H", source["main_span_H"]) * MM
    w = _positive("diagonal_width_w", source["diagonal_width_w"]) * MM
    b = _positive("vertical_segment_width_b", source["vertical_segment_width_b"]) * MM
    P = _positive("clear_pitch_P", source["clear_pitch_P"]) * MM
    d = _positive("single_layer_depth_d", source["single_layer_depth_d"]) * MM
    n_units = int(source["number_of_units"])
    inlet_lead = _positive("inlet_lead_length", source["inlet_lead_length"]) * MM
    outlet_lead = _positive("outlet_lead_length", source["outlet_lead_length"]) * MM
    inset_ratio = float(source.get("diagonal_end_inset_over_w", 0.5))
    phase_ratio = float(source.get("crossing_phase_over_b", 0.0))
    gap = float(source.get("interlayer_gap", 0.0)) * MM

    if n_units < 1:
        raise ValueError("number_of_units must be at least one")
    if w >= H:
        raise ValueError("diagonal_width_w must be smaller than main_span_H")
    if not 0.5 <= inset_ratio <= 1.5:
        raise ValueError("diagonal_end_inset_over_w must lie in [0.5, 1.5]")
    if not -0.4 <= phase_ratio <= 0.4:
        raise ValueError("crossing_phase_over_b must lie in [-0.4, 0.4]")
    if gap < 0.0:
        raise ValueError("interlayer_gap cannot be negative")
    if gap > 0.0:
        raise ValueError(
            "the validated baseline requires face-adjacent layers; a positive "
            "gap needs explicit via solids and is not implemented"
        )

    inset = inset_ratio * w
    y_low = -0.5 * H + inset
    y_high = 0.5 * H - inset
    if y_high <= y_low:
        raise ValueError("diagonal endpoint inset closes the transverse span")

    core_x0 = 0.0
    bar_centres = [0.5 * b + index * (P + b) for index in range(n_units + 1)]
    core_x1 = (n_units + 1) * b + n_units * P
    inlet_x = core_x0 - inlet_lead
    outlet_x = core_x1 + outlet_lead
    z_bottom = -d
    z_mid = 0.0
    z_top = d
    phase = phase_ratio * b

    upper: list[cq.Workplane] = []
    lower: list[cq.Workplane] = []
    upper_verticals: list[cq.Workplane] = []
    lower_verticals: list[cq.Workplane] = []
    for centre in bar_centres:
        upper_verticals.append(
            _box(centre - 0.5 * b, centre + 0.5 * b, -0.5 * H, 0.5 * H, z_mid, z_top)
        )
        lower_verticals.append(
            _box(centre - 0.5 * b, centre + 0.5 * b, -0.5 * H, 0.5 * H, z_bottom, z_mid)
        )
    upper.extend(upper_verticals)
    lower.extend(lower_verticals)

    for index in range(n_units):
        left = bar_centres[index]
        right = bar_centres[index + 1]
        upper.append(_segment((left + phase, y_high), (right + phase, y_low), w, z_mid, d))
        lower.append(_segment((left - phase, y_low), (right - phase, y_high), w, z_bottom, d))

    upper_inlet_y = 0.25 * H
    lower_inlet_y = -0.25 * H
    upper.append(_box(inlet_x, core_x0 + b, upper_inlet_y - 0.5 * w, upper_inlet_y + 0.5 * w, z_mid, z_top))
    lower.append(_box(inlet_x, core_x0 + b, lower_inlet_y - 0.5 * w, lower_inlet_y + 0.5 * w, z_bottom, z_mid))
    outlet = _box(core_x1 - b, outlet_x, -0.5 * w, 0.5 * w, z_bottom, z_top)

    upper_solid = _fuse_all(upper)
    lower_solid = _fuse_all(lower)
    # All channel side walls are normal to the layer interface. Translating the
    # upper solid downward by a thin probe depth therefore gives the exact open
    # interface area as intersection volume divided by depth. The outlet is
    # excluded so this measures only the repeated mixing-core connections.
    interface_probe_depth = 0.1 * d
    interface_open_area = (
        upper_solid.translate((0.0, 0.0, -interface_probe_depth))
        .intersect(lower_solid)
        .val()
        .Volume()
        / interface_probe_depth
    )
    vertical_open_area = (
        _fuse_all(upper_verticals)
        .translate((0.0, 0.0, -interface_probe_depth))
        .intersect(_fuse_all(lower_verticals))
        .val()
        .Volume()
        / interface_probe_depth
    )
    crossing_open_area = max(interface_open_area - vertical_open_area, 0.0)

    fluid = _fuse_all([upper_solid, lower_solid, outlet])
    all_faces = list(fluid.val().Faces())
    position_tolerance = max(1.0e-12, (outlet_x - inlet_x) * 1.0e-9)

    def on_x_plane(face: cq.Face, x_value: float) -> bool:
        bounds = face.BoundingBox()
        return abs(bounds.xmin - x_value) <= position_tolerance and abs(bounds.xmax - x_value) <= position_tolerance

    inlet_faces = [face for face in all_faces if on_x_plane(face, inlet_x)]
    outlet_faces = [face for face in all_faces if on_x_plane(face, outlet_x)]
    if len(inlet_faces) != 2:
        raise RuntimeError(f"expected two inlet faces, found {len(inlet_faces)}")
    if len(outlet_faces) != 1:
        raise RuntimeError(f"expected one outlet face, found {len(outlet_faces)}")

    inlet_faces.sort(key=lambda face: face.Center().z, reverse=True)
    inlet1_faces = [inlet_faces[0]]
    inlet2_faces = [inlet_faces[1]]
    claimed_ids = {id(face) for face in inlet_faces + outlet_faces}
    wall_faces = [face for face in all_faces if id(face) not in claimed_ids]
    patch_faces = {
        "inlet1": inlet1_faces,
        "inlet2": inlet2_faces,
        "outlet": outlet_faces,
        "walls": wall_faces,
    }
    if sum(len(faces) for faces in patch_faces.values()) != len(all_faces):
        raise RuntimeError("each CAD face must belong to exactly one patch")

    areas = {name: sum(face.Area() for face in faces) for name, faces in patch_faces.items()}
    expected_inlet_area = w * d
    expected_outlet_area = w * (2.0 * d)
    for name in ("inlet1", "inlet2"):
        if not math.isclose(areas[name], expected_inlet_area, rel_tol=1.0e-7, abs_tol=1.0e-16):
            raise RuntimeError(f"{name} area does not match w*d")
    if not math.isclose(areas["outlet"], expected_outlet_area, rel_tol=1.0e-7, abs_tol=1.0e-16):
        raise RuntimeError("outlet area does not match w*(2d)")

    hydraulic_diameter = 2.0 * w * d / (w + d)
    diagonal_span = y_high - y_low
    crossing_angle_deg = math.degrees(math.atan2(diagonal_span, P + b))
    manifest = {
        "schema_version": 1,
        "topology": "m10_inspired_two_layer_serpentine_crossing",
        "source": {
            "citation": "Hossain et al., Chemical Engineering Journal 327 (2017) 268-277",
            "doi": "10.1016/j.cej.2017.06.106",
            "classification": "m10_inspired_reconstruction",
            "exact_reproduction_claim": False,
            "reconstruction_assumptions": [
                "published sources do not provide CAD, mask, or aperture dimensions",
                "short axial lead lengths are not uniquely specified by the paper",
                "rectangular sharp-cornered channel unions implement the published schematic",
                "complete projected overlaps are open at X nodes and vertical segments",
            ],
        },
        "dimensions_m": {
            "H": H,
            "w": w,
            "b": b,
            "P": P,
            "d": d,
            "D": 2.0 * d,
            "inlet_lead": inlet_lead,
            "outlet_lead": outlet_lead,
            "core_length": core_x1 - core_x0,
            "total_axial_length": outlet_x - inlet_x,
            "inlet_hydraulic_diameter": hydraulic_diameter,
        },
        "derived": {
            "number_of_units": n_units,
            "unit_axial_length": P + b,
            "crossing_angle_deg": crossing_angle_deg,
            "crossing_phase_m": phase,
            "diagonal_end_inset_m": inset,
            # Strictly inside the upper inlet lead and deliberately offset
            # from symmetry planes.  run_case.py perturbs it further if a
            # background blockMesh plane happens to coincide with it.
            "snappy_location_in_mesh_m": [
                inlet_x + 0.43 * inlet_lead,
                upper_inlet_y + 0.07 * w,
                0.37 * d,
            ],
        },
        "interlayer_connections": {
            "definition": "projected_open_area_at_layer_interface_excluding_outlet",
            "probe_depth_m": interface_probe_depth,
            "total_open_area_m2": interface_open_area,
            "vertical_segment_open_area_m2": vertical_open_area,
            "crossing_open_area_outside_vertical_segments_m2": crossing_open_area,
            "vertical_segment_count": n_units + 1,
            "crossing_unit_count": n_units,
            "total_open_area_over_one_inlet_area": interface_open_area / expected_inlet_area,
        },
        "patches": {
            name: {
                "cad_face_count": len(faces),
                "area_m2": areas[name],
                "bounds_m": _combined_bounds(faces),
            }
            for name, faces in patch_faces.items()
        },
    }
    manifest["_patch_faces"] = patch_faces
    return fluid, manifest


def write_outputs(config_path: Path, output_override: Path | None = None) -> dict:
    with config_path.open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    fluid, manifest = build_geometry(config)
    patch_faces = manifest.pop("_patch_faces")
    output_cfg = config["output"]
    configured = Path(output_cfg["directory"])
    output_dir = output_override or (config_path.resolve().parents[1] / configured)
    output_dir.mkdir(parents=True, exist_ok=True)

    stl_path = output_dir / output_cfg["stl_name"]
    step_path = output_dir / output_cfg["step_name"]
    manifest_path = output_dir / output_cfg["manifest_name"]
    tolerance = min(manifest["dimensions_m"]["w"], manifest["dimensions_m"]["d"]) / 20.0
    with stl_path.open("w", encoding="ascii") as stream:
        for patch_name, faces in patch_faces.items():
            stream.write(_stl_block(faces, patch_name, tolerance=tolerance))
    cq.exporters.export(fluid, str(step_path))
    manifest["outputs"] = {
        "stl": stl_path.name,
        "step": step_path.name,
        "manifest": manifest_path.name,
    }
    with manifest_path.open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_suffix(".yaml"),
        help="geometry YAML (default: file beside this script)",
    )
    parser.add_argument("--output-dir", type=Path, help="override the YAML output directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = write_outputs(args.config.resolve(), args.output_dir.resolve() if args.output_dir else None)
    dimensions = manifest["dimensions_m"]
    print("Generated M10-inspired two-layer serpentine-crossing fluid domain")
    print(f"  units: {manifest['derived']['number_of_units']}")
    print(f"  axial length: {dimensions['total_axial_length'] * 1.0e3:.3f} mm")
    print(f"  inlet hydraulic diameter: {dimensions['inlet_hydraulic_diameter'] * 1.0e6:.1f} um")
    print(f"  crossing angle: {manifest['derived']['crossing_angle_deg']:.2f} deg")


if __name__ == "__main__":
    main()
