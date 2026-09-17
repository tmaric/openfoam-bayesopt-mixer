#!/usr/bin/env python3
"""Shared CAD utilities for the TN-05 two-square channel benchmark."""

from __future__ import annotations

from pathlib import Path

import cadquery as cq
import yaml


CFG_NAME = "channel_two_square_obstacles.yaml"
GEOMETRY_MANIFEST_NAME = "channel_two_square_obstacles_geometry_manifest.yaml"
STL_NAME = "channel_two_square_obstacles.stl"
TOL = 1e-6


def load_config(config_path: Path | None = None) -> dict:
    path = config_path or Path(__file__).resolve().with_name(CFG_NAME)
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def _stl_block(faces, solid_name: str, tolerance: float = 1.0e-4, angular_tol: float = 0.1) -> str:
    if not faces:
        raise ValueError(f"No faces found for patch '{solid_name}'")

    compound = cq.Compound.makeCompound(faces)
    verts, tris = compound.tessellate(tolerance, angular_tol)
    lines = [f"solid {solid_name}"]
    for tri in tris:
        v0, v1, v2 = verts[tri[0]], verts[tri[1]], verts[tri[2]]
        e1x, e1y, e1z = v1.x - v0.x, v1.y - v0.y, v1.z - v0.z
        e2x, e2y, e2z = v2.x - v0.x, v2.y - v0.y, v2.z - v0.z
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        nl = (nx * nx + ny * ny + nz * nz) ** 0.5
        if nl > 1.0e-30:
            nx, ny, nz = nx / nl, ny / nl, nz / nl
        lines += [
            f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}",
            "    outer loop",
            f"      vertex {v0.x:.6e} {v0.y:.6e} {v0.z:.6e}",
            f"      vertex {v1.x:.6e} {v1.y:.6e} {v1.z:.6e}",
            f"      vertex {v2.x:.6e} {v2.y:.6e} {v2.z:.6e}",
            "    endloop",
            "  endfacet",
        ]
    lines.append(f"endsolid {solid_name}")
    return "\n".join(lines) + "\n"


def compute_resolved_geometry(raw: dict, mode: str) -> dict:
    scale = float(raw["scale"])
    h_norm = float(raw["H"])
    inlet_buffer = float(raw["inlet_buffer"])
    outlet_buffer = float(raw["outlet_buffer"])
    a_max = float(raw["a_max_design"])
    d_max = float(raw["d_max_design"])
    depth_ratio = float(raw["depth_ratio"])
    spacing_margin = float(raw["spacing_margin_H"])

    a = float(raw["a"])
    if not (0.0 < a <= a_max):
        raise ValueError(f"Obstacle size a={a:.6f} must satisfy 0 < a <= {a_max:.6f}")

    d_ratio = None
    if mode == "constrained":
        d_ratio = _clamp01(raw["d_ratio"])
        d_slack_max = d_max - a - spacing_margin
        if d_slack_max < -1.0e-12:
            raise ValueError(
                f"Constrained spacing interval is empty for a={a:.6f}: "
                f"d_max={d_max:.6f}, spacing_margin_H={spacing_margin:.6f}"
            )
        d = a + spacing_margin + d_ratio * max(d_slack_max, 0.0)
    elif mode == "unconstrained":
        d = float(raw["d"])
    else:
        raise ValueError(f"Unsupported CAD mode '{mode}'")

    if not (0.0 < d <= d_max):
        raise ValueError(f"Obstacle spacing d={d:.6f} must satisfy 0 < d <= {d_max:.6f}")

    h_phys = h_norm * scale
    side_phys = a * h_phys
    depth_phys = depth_ratio * h_phys

    x1_center = inlet_buffer + 0.5 * a_max
    x2_center = x1_center + d
    y_center = 0.5 * h_norm
    total_length = inlet_buffer + a_max + d_max + outlet_buffer

    resolved = {
        "cad_mode": mode,
        "scale": scale,
        "H": h_norm,
        "a": a,
        "d": d,
        "feasible": bool(d > a),
        "inlet_buffer": inlet_buffer,
        "outlet_buffer": outlet_buffer,
        "a_max_design": a_max,
        "d_max_design": d_max,
        "depth_ratio": depth_ratio,
        "spacing_margin_H": spacing_margin,
        "x1_center": x1_center,
        "x2_center": x2_center,
        "y_center": y_center,
        "total_length": total_length,
        "a_phys": side_phys,
        "d_phys": d * h_phys,
        "H_phys": h_phys,
        "depth_phys": depth_phys,
        "total_length_phys": total_length * h_phys,
        "x1_center_phys": x1_center * h_phys,
        "x2_center_phys": x2_center * h_phys,
        "y_center_phys": y_center * h_phys,
    }
    if d_ratio is not None:
        resolved["d_ratio"] = d_ratio
    return resolved


def _rect_points(x0: float, x1: float, y0: float, y1: float) -> list[tuple[float, float]]:
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def _make_extruded_polygon(points_2d: list[tuple[float, float]], depth: float, eps: float):
    wire = (
        cq.Workplane("XY")
        .workplane(offset=-2.0 * eps)
        .polyline(points_2d)
        .close()
    )
    return wire.extrude(depth + 4.0 * eps)


def _is_inlet(face) -> bool:
    return face.normalAt().x < -(1.0 - TOL)


def _is_outlet(face) -> bool:
    return face.normalAt().x > (1.0 - TOL)


def _is_wall(face) -> bool:
    normal = face.normalAt()
    return not (abs(normal.x) > (1.0 - TOL)) and not (abs(normal.z) > (1.0 - TOL))


def write_geometry_manifest(resolved: dict) -> Path:
    manifest_path = Path(__file__).resolve().parent / GEOMETRY_MANIFEST_NAME
    with open(manifest_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(resolved, handle, default_flow_style=False, sort_keys=False)
    return manifest_path


def export_geometry(mode: str, config_path: Path | None = None) -> dict:
    raw = load_config(config_path)
    resolved = compute_resolved_geometry(raw, mode)
    manifest_path = write_geometry_manifest(resolved)

    if mode == "unconstrained" and not resolved["feasible"]:
        print(f"Written geometry manifest: {manifest_path}")
        raise ValueError(
            f"Invalid TN-05 benchmark geometry: d={resolved['d']:.6f} <= a={resolved['a']:.6f}. "
            "The unconstrained generator must fail in the infeasible overlap region."
        )

    h_phys = resolved["H_phys"]
    depth_phys = resolved["depth_phys"]
    total_length_phys = resolved["total_length_phys"]
    side_phys = resolved["a_phys"]
    y_center_phys = resolved["y_center_phys"]
    x1_center_phys = resolved["x1_center_phys"]
    x2_center_phys = resolved["x2_center_phys"]

    eps = depth_phys * 0.01
    out_dir = Path(__file__).resolve().parent / "constant" / "triSurface"
    out_dir.mkdir(parents=True, exist_ok=True)

    obstacle_solids = []
    for x_center in [x1_center_phys, x2_center_phys]:
        x0 = x_center - 0.5 * side_phys
        x1 = x_center + 0.5 * side_phys
        y0 = y_center_phys - 0.5 * side_phys
        y1 = y_center_phys + 0.5 * side_phys
        obstacle_solids.append(
            _make_extruded_polygon(_rect_points(x0, x1, y0, y1), depth_phys, eps)
        )

    channel_box = (
        cq.Workplane("XY")
        .box(total_length_phys, h_phys, depth_phys + 2.0 * eps)
        .translate((total_length_phys / 2.0, h_phys / 2.0, depth_phys / 2.0))
    )

    fluid = channel_box
    for obstacle in obstacle_solids:
        fluid = fluid.cut(obstacle)

    all_faces = fluid.val().Faces()
    patch_faces = {
        "inlet": [face for face in all_faces if _is_inlet(face)],
        "outlet": [face for face in all_faces if _is_outlet(face)],
        "walls": [face for face in all_faces if _is_wall(face)],
    }

    stl_path = out_dir / STL_NAME
    with open(stl_path, "w", encoding="utf-8") as handle:
        for patch_name, faces in patch_faces.items():
            handle.write(_stl_block(faces, patch_name))

    print(f"Written STL: {stl_path}")
    print(f"Written geometry manifest: {manifest_path}")
    print(
        "Geometry summary: "
        f"mode={mode}, a={resolved['a']:.4f}, d={resolved['d']:.4f}, "
        f"x1={resolved['x1_center']:.4f}, x2={resolved['x2_center']:.4f}, "
        f"L={resolved['total_length']:.4f}"
    )
    return resolved
