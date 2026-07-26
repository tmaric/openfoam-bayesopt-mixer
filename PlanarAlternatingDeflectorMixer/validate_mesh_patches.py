#!/usr/bin/env python3
"""Validate inlet/outlet location and area directly from an ASCII polyMesh.

This check deliberately does not rely on patch normals.  It reconstructs the
boundary polygons from ``points``, ``faces`` and ``boundary`` and verifies
that inlet and outlet are confined to the two exterior x-planes recorded by
the CAD geometry manifest.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


def _without_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*?$", "", text, flags=re.MULTILINE)


def _counted_list_body(path: Path) -> str:
    text = _without_comments(path.read_text())
    match = re.search(r"\b\d+\s*\n?\s*\(", text)
    if not match:
        raise ValueError(f"could not find counted list in {path}")
    start = text.find("(", match.start()) + 1
    depth = 1
    for index in range(start, len(text)):
        char = text[index]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                return text[start:index]
    raise ValueError(f"unterminated counted list in {path}")


def _read_points(path: Path) -> list[tuple[float, float, float]]:
    body = _counted_list_body(path)
    points = []
    for values in re.findall(r"\(([^()]*)\)", body):
        xyz = values.split()
        if len(xyz) == 3:
            points.append(tuple(float(value) for value in xyz))
    if not points:
        raise ValueError(f"no points parsed from {path}")
    return points


def _read_boundary(path: Path) -> dict[str, dict[str, int | str]]:
    body = _counted_list_body(path)
    patches: dict[str, dict[str, int | str]] = {}
    for match in re.finditer(r"([A-Za-z_][\w.]*)\s*\{([^{}]*)\}", body, re.DOTALL):
        name, content = match.groups()
        n_faces = re.search(r"\bnFaces\s+(\d+)\s*;", content)
        start_face = re.search(r"\bstartFace\s+(\d+)\s*;", content)
        patch_type = re.search(r"\btype\s+([^;\s]+)\s*;", content)
        if n_faces and start_face and patch_type:
            patches[name] = {
                "nFaces": int(n_faces.group(1)),
                "startFace": int(start_face.group(1)),
                "type": patch_type.group(1),
            }
    if not patches:
        raise ValueError(f"no boundary patches parsed from {path}")
    return patches


def _read_selected_faces(path: Path, wanted: set[int]) -> dict[int, list[int]]:
    body = _counted_list_body(path)
    selected: dict[int, list[int]] = {}
    face_pattern = re.compile(r"(\d+)\s*\(([^()]*)\)")
    for face_index, match in enumerate(face_pattern.finditer(body)):
        if face_index not in wanted:
            continue
        expected = int(match.group(1))
        labels = [int(value) for value in match.group(2).split()]
        if len(labels) != expected:
            raise ValueError(
                f"face {face_index} declares {expected} points but contains {len(labels)}"
            )
        selected[face_index] = labels
    missing = wanted.difference(selected)
    if missing:
        preview = ", ".join(str(value) for value in sorted(missing)[:8])
        raise ValueError(f"could not parse requested face indices: {preview}")
    return selected


def _subtract(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _polygon_area(vertices) -> float:
    anchor = vertices[0]
    area = 0.0
    for index in range(1, len(vertices) - 1):
        cross = _cross(
            _subtract(vertices[index], anchor),
            _subtract(vertices[index + 1], anchor),
        )
        area += 0.5 * math.sqrt(sum(component * component for component in cross))
    return area


def _patch_metrics(patch, faces, points) -> dict:
    start = int(patch["startFace"])
    end = start + int(patch["nFaces"])
    vertex_ids = []
    area = 0.0
    for face_index in range(start, end):
        labels = faces[face_index]
        vertices = [points[label] for label in labels]
        vertex_ids.extend(labels)
        area += _polygon_area(vertices)
    vertices = [points[label] for label in set(vertex_ids)]
    return {
        "n_faces": int(patch["nFaces"]),
        "area_m2": area,
        "bounds_m": {
            "xmin": min(point[0] for point in vertices),
            "xmax": max(point[0] for point in vertices),
            "ymin": min(point[1] for point in vertices),
            "ymax": max(point[1] for point in vertices),
            "zmin": min(point[2] for point in vertices),
            "zmax": max(point[2] for point in vertices),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--research-config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import yaml

    manifest = json.loads(args.manifest.read_text())
    config = yaml.safe_load(args.research_config.read_text())
    validation = config["validation"]
    height = float(manifest["channel_height_m"])
    length = float(manifest["total_length_m"])
    position_tol = height * float(validation["mesh_boundary_position_tolerance_H"])
    area_rel_tol = float(validation["mesh_patch_area_relative_tolerance"])

    mesh = args.case / "constant" / "polyMesh"
    points = _read_points(mesh / "points")
    boundary = _read_boundary(mesh / "boundary")
    required = {"inlet", "outlet", "walls", "frontAndBack"}
    missing = required.difference(boundary)
    if missing:
        raise RuntimeError(f"mesh is missing required patches: {sorted(missing)}")
    if boundary["walls"]["type"] != "wall":
        raise RuntimeError("walls patch is not type wall")
    if boundary["frontAndBack"]["type"] != "empty":
        raise RuntimeError("frontAndBack patch is not type empty")

    wanted = set()
    for name in ("inlet", "outlet"):
        start = int(boundary[name]["startFace"])
        wanted.update(range(start, start + int(boundary[name]["nFaces"])))
    faces = _read_selected_faces(mesh / "faces", wanted)
    metrics = {
        name: _patch_metrics(boundary[name], faces, points)
        for name in ("inlet", "outlet")
    }

    checks = []
    for name, target_x in (("inlet", 0.0), ("outlet", length)):
        bounds = metrics[name]["bounds_m"]
        expected_area = float(manifest["patches"][name]["area_m2"])
        checks.extend(
            [
                (
                    abs(bounds["xmin"] - target_x) <= position_tol
                    and abs(bounds["xmax"] - target_x) <= position_tol,
                    f"{name} must be confined to x={target_x:.12g} m; bounds={bounds}",
                ),
                (
                    math.isclose(
                        metrics[name]["area_m2"],
                        expected_area,
                        rel_tol=area_rel_tol,
                        abs_tol=1e-16,
                    ),
                    f"{name} area {metrics[name]['area_m2']:.12e} m2 differs from "
                    f"CAD area {expected_area:.12e} m2",
                ),
            ]
        )

    failures = [message for passed, message in checks if not passed]
    report = {
        "schema_version": 1,
        "passed": not failures,
        "position_tolerance_m": position_tol,
        "area_relative_tolerance": area_rel_tol,
        "patches": metrics,
        "failures": failures,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if failures:
        raise RuntimeError("mesh boundary validation failed:\n  " + "\n  ".join(failures))
    print(
        "[mesh-validation] PASS: inlet/outlet occupy only x=0 and x=L; "
        "areas agree with CAD"
    )


if __name__ == "__main__":
    main()
