#!/usr/bin/env python3
"""Render the latest reconstructed OpenFOAM mixing field without ParaView.

foamlib discovers and validates the latest field.  VTK's native OpenFOAM
reader supplies the polyhedral mesh and cell data, while Pillow rasterizes the
top faces directly.  No OpenGL context or display server is required.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import vtk
from foamlib import FoamCase
from PIL import Image, ImageDraw, ImageFont
from vtk.util.numpy_support import vtk_to_numpy


FIELD_MIN = 0.0
FIELD_MAX = 1.0
DEFAULT_WIDTH = 16_000
GEOMETRY_PADDING = 6


def load_font(size: int, bold: bool = False):
    """Load a common Linux font with a Pillow fallback."""
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    candidates = (
        Path("/usr/share/fonts/truetype/dejavu") / name,
        Path("/usr/share/fonts/dejavu") / name,
    )
    for candidate in candidates:
        try:
            return ImageFont.truetype(str(candidate), size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def field_display_name(field: str) -> str:
    """Return the label used below the scalar color bar."""
    return "Temperature" if field == "T" else field


def rainbow_lookup(size: int = 4096) -> np.ndarray:
    """Return the legacy blue-cyan-green-yellow-red transfer function."""
    values = np.linspace(0.0, 1.0, size)
    stops = np.array([0.0, 0.25, 0.50, 0.75, 1.0])
    red = np.interp(values, stops, [0.0, 0.0, 0.0, 1.0, 1.0])
    green = np.interp(values, stops, [0.0, 1.0, 1.0, 1.0, 0.0])
    blue = np.interp(values, stops, [1.0, 1.0, 0.0, 0.0, 0.0])
    return np.rint(255.0 * np.stack((red, green, blue), axis=1)).astype(np.uint8)


def append_colorbar(
    png_path: Path,
    field: str,
    field_min: float = FIELD_MIN,
    field_max: float = FIELD_MAX,
) -> None:
    """Append a colorbar below the geometry using Pillow only."""
    geometry = Image.open(png_path).convert("RGB")
    width = geometry.width

    bar_height_px = 340
    colorbar = Image.new("RGB", (width, bar_height_px), "white")
    draw = ImageDraw.Draw(colorbar)

    margin_x = max(120, int(0.08 * width))
    gradient_width = max(2, width - 2 * margin_x)
    gradient_height = 86
    gradient_top = 38

    lookup = rainbow_lookup(gradient_width)
    gradient = np.tile(lookup[np.newaxis, :, :], (gradient_height, 1, 1))
    colorbar.paste(Image.fromarray(gradient, mode="RGB"), (margin_x, gradient_top))
    draw.rectangle(
        (
            margin_x,
            gradient_top,
            margin_x + gradient_width - 1,
            gradient_top + gradient_height,
        ),
        outline="black",
        width=3,
    )

    tick_font = load_font(38)
    label_font = load_font(52, bold=True)
    tick_values = np.linspace(field_min, field_max, 5)
    tick_y = gradient_top + gradient_height
    for fraction, tick in zip(np.linspace(0.0, 1.0, 5), tick_values):
        x = margin_x + int(fraction * (gradient_width - 1))
        draw.line((x, tick_y, x, tick_y + 14), fill="black", width=3)
        label = f"{tick:.2f}"
        bbox = draw.textbbox((0, 0), label, font=tick_font)
        label_width = bbox[2] - bbox[0]
        draw.text((x - label_width / 2, tick_y + 18), label, fill="black", font=tick_font)

    title = field_display_name(field)
    bbox = draw.textbbox((0, 0), title, font=label_font)
    title_width = bbox[2] - bbox[0]
    draw.text(((width - title_width) / 2, 245), title, fill="black", font=label_font)

    gap = 24
    combined = Image.new("RGB", (width, geometry.height + gap + colorbar.height), "white")
    combined.paste(geometry, (0, 0))
    combined.paste(colorbar, (0, geometry.height + gap))
    combined.save(png_path)


def latest_field(case_dir: Path, field: str) -> tuple[float, np.ndarray]:
    """Use foamlib to locate and validate the latest cell field."""
    case = FoamCase(case_dir)
    latest = case[-1]
    values = np.asarray(latest[field].internal_field)
    if values.ndim != 1:
        raise ValueError(
            f"Expected scalar field {field!r} at time {latest.name}, "
            f"got shape {values.shape}."
        )
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Field {field!r} at time {latest.name} contains non-finite values.")
    return float(latest.name), values


def internal_mesh(reader: vtk.vtkOpenFOAMReader):
    """Return the internal-mesh dataset from vtkOpenFOAMReader output."""
    output = reader.GetOutput()
    fallback = None

    def visit(block):
        nonlocal fallback
        if block is None:
            return None
        if block.IsA("vtkDataSet"):
            fallback = fallback or block
            return None
        if not block.IsA("vtkMultiBlockDataSet"):
            return None
        for index in range(block.GetNumberOfBlocks()):
            child = block.GetBlock(index)
            name = None
            if block.HasMetaData(index):
                name = block.GetMetaData(index).Get(vtk.vtkCompositeDataSet.NAME())
            if name == "internalMesh" and child is not None and child.IsA("vtkDataSet"):
                return child
            match = visit(child)
            if match is not None:
                return match
        return None

    match = visit(output)
    if match is not None:
        return match
    if fallback is not None:
        return fallback
    raise RuntimeError("VTK did not produce an internal OpenFOAM mesh.")


def read_vtk_mesh(case_foam: Path, field: str, time_value: float):
    """Read the requested reconstructed OpenFOAM time with Python VTK."""
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(str(case_foam))
    reader.UpdateInformation()
    reader.DisableAllCellArrays()
    reader.SetCellArrayStatus(field, 1)
    reader.SetTimeValue(time_value)
    reader.Update()

    mesh = internal_mesh(reader)
    vtk_field = mesh.GetCellData().GetArray(field)
    if vtk_field is None:
        available = [
            mesh.GetCellData().GetArrayName(i)
            for i in range(mesh.GetCellData().GetNumberOfArrays())
        ]
        raise KeyError(f"Cell field {field!r} not found by VTK; available fields: {available}")
    return mesh


def surface_with_cell_field(mesh):
    """Extract boundary faces while retaining their originating cell values."""
    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputData(mesh)
    surface_filter.Update()
    surface = surface_filter.GetOutput()
    if surface.GetNumberOfCells() == 0:
        raise RuntimeError("VTK surface extraction produced no cells.")
    return surface


def render_top_faces(
    mesh,
    surface,
    output_png: Path,
    field: str,
    width: int,
    field_min: float,
    field_max: float,
) -> int:
    """Rasterize the top faces of the thin x-y mesh with Pillow."""
    if width < 800:
        raise ValueError("Image width must be at least 800 pixels.")
    if field_max <= field_min:
        raise ValueError("field_max must be greater than field_min.")

    points = vtk_to_numpy(surface.GetPoints().GetData())
    cell_values = vtk_to_numpy(surface.GetCellData().GetArray(field))
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.GetBounds()
    x_span = xmax - xmin
    y_span = ymax - ymin
    z_span = zmax - zmin
    if x_span <= 0.0 or y_span <= 0.0 or z_span <= 0.0:
        raise ValueError(f"Invalid mesh bounds: {mesh.GetBounds()}")

    usable_width = width - 2 * GEOMETRY_PADDING
    geometry_height = max(
        160,
        int(round(usable_width * y_span / x_span)) + 2 * GEOMETRY_PADDING,
    )
    usable_height = geometry_height - 2 * GEOMETRY_PADDING
    image = Image.new("RGB", (width, geometry_height), "white")
    draw = ImageDraw.Draw(image)

    lookup = rainbow_lookup()
    denominator = field_max - field_min
    top_tolerance = max(1.0e-3 * z_span, 1.0e-12)
    face_count = 0

    for cell_id in range(surface.GetNumberOfCells()):
        cell = surface.GetCell(cell_id)
        point_ids = [cell.GetPointId(index) for index in range(cell.GetNumberOfPoints())]
        xyz = points[point_ids]
        # Select only the front/top empty-patch face.  An average-z test also
        # admits some side polygons after face triangulation, whereas every
        # vertex of the desired face lies on z_max.
        if float(np.min(xyz[:, 2])) < zmax - top_tolerance:
            continue

        vertices = [
            (
                GEOMETRY_PADDING + int(round((x - xmin) / x_span * usable_width)),
                GEOMETRY_PADDING + int(round((ymax - y) / y_span * usable_height)),
            )
            for x, y in xyz[:, :2]
        ]
        normalized = np.clip((float(cell_values[cell_id]) - field_min) / denominator, 0.0, 1.0)
        color = tuple(int(channel) for channel in lookup[int(normalized * (len(lookup) - 1))])
        draw.polygon(vertices, fill=color, outline=color)
        face_count += 1

    if face_count == 0:
        raise RuntimeError("No top faces were found in the thin x-y mesh.")

    output_png.parent.mkdir(parents=True, exist_ok=True)
    image.save(output_png)
    append_colorbar(output_png, field, field_min, field_max)
    return face_count


def render_latest_field(
    case_foam: Path,
    output_png: Path,
    field: str,
    width: int = DEFAULT_WIDTH,
    field_min: float = FIELD_MIN,
    field_max: float = FIELD_MAX,
) -> None:
    """Render the latest reconstructed field to a PNG without OpenGL."""
    if not case_foam.is_file():
        raise FileNotFoundError(case_foam)

    time_value, foamlib_values = latest_field(case_foam.parent, field)
    mesh = read_vtk_mesh(case_foam, field, time_value)
    if mesh.GetNumberOfCells() != foamlib_values.size:
        raise ValueError(
            "foamlib/VTK cell-count mismatch: "
            f"{foamlib_values.size} field values versus {mesh.GetNumberOfCells()} mesh cells."
        )

    surface = surface_with_cell_field(mesh)
    face_count = render_top_faces(
        mesh,
        surface,
        output_png,
        field,
        width,
        field_min,
        field_max,
    )
    print(
        f"Rendered {field} at time {time_value:g}: {mesh.GetNumberOfCells()} cells, "
        f"{face_count} top faces -> {output_png}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-foam", required=True, help="Path to the OpenFOAM .foam marker")
    parser.add_argument("--output", required=True, help="PNG output path")
    parser.add_argument("--field", default="T", help="Scalar cell field (default: T)")
    parser.add_argument("--width", type=int, default=DEFAULT_WIDTH, help="Image width in pixels")
    parser.add_argument("--field-min", type=float, default=FIELD_MIN)
    parser.add_argument("--field-max", type=float, default=FIELD_MAX)
    args = parser.parse_args()

    render_latest_field(
        Path(args.case_foam).resolve(),
        Path(args.output).resolve(),
        args.field,
        args.width,
        args.field_min,
        args.field_max,
    )


if __name__ == "__main__":
    main()
