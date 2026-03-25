#!/usr/bin/env pvpython
"""Render the latest reconstructed mixing field with ParaView in batch mode."""

import argparse
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy
from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    OpenFOAMReader,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    UpdatePipeline,
    _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

FIELD_MIN = 0.0
FIELD_MAX = 1.0


def latest_time(reader) -> float | None:
    """Return the latest available timestep, or None when the case is steady-only."""
    reader.UpdatePipelineInformation()
    timesteps = list(getattr(reader, 'TimestepValues', []) or [])
    return timesteps[-1] if timesteps else None


def set_rainbow_transfer_function(lut, field_min: float = 0.0, field_max: float = 1.0) -> None:
    """Use a simple rainbow transfer function over the expected T range [0, 1]."""
    lut.ColorSpace = 'RGB'
    lut.RGBPoints = [
        field_min, 0.0, 0.0, 1.0,
        field_min + 0.25 * (field_max - field_min), 0.0, 1.0, 1.0,
        field_min + 0.50 * (field_max - field_min), 0.0, 1.0, 0.0,
        field_min + 0.75 * (field_max - field_min), 1.0, 1.0, 0.0,
        field_max, 1.0, 0.0, 0.0,
    ]
    lut.NanColor = [0.8, 0.8, 0.8]
    lut.RescaleTransferFunction(field_min, field_max)


def rainbow_colormap() -> LinearSegmentedColormap:
    """Matplotlib colormap matching the ParaView rainbow transfer function."""
    return LinearSegmentedColormap.from_list(
        "sar_rainbow",
        [
            (0.0, 0.0, 1.0),
            (0.0, 1.0, 1.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
            (1.0, 0.0, 0.0),
        ],
        N=2048,
    )


def set_xy_camera(render_view, bounds: tuple[float, ...]) -> None:
    """Look down the z-axis so the thin mixer is shown in the x-y plane."""
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    x_span = max(xmax - xmin, 1e-9)
    y_span = max(ymax - ymin, 1e-9)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)

    render_view.InteractionMode = '2D'
    render_view.CameraParallelProjection = 1
    render_view.CameraFocalPoint = [cx, cy, cz]
    render_view.CameraPosition = [cx, cy, cz + max(x_span, y_span, 1.0)]
    render_view.CameraViewUp = [0.0, 1.0, 0.0]


def fit_xy_camera(render_view, bounds: tuple[float, ...]) -> None:
    """Fit the x-y slice tightly in parallel projection without ParaView's loose reset."""
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    x_span = max(xmax - xmin, 1e-9)
    y_span = max(ymax - ymin, 1e-9)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)

    set_xy_camera(render_view, bounds)

    geometry_width_fraction = 0.98
    geometry_height_fraction = 0.78
    padding_factor = 1.01
    aspect = render_view.ViewSize[0] / render_view.ViewSize[1]
    scale = padding_factor * max(
        0.5 * y_span / geometry_height_fraction,
        0.5 * x_span / (aspect * geometry_width_fraction),
    )
    free_vertical_half_span = max(scale - 0.5 * y_span, 0.0)
    vertical_shift = 0.55 * free_vertical_half_span

    render_view.CameraParallelScale = scale
    render_view.CameraFocalPoint = [cx, cy - vertical_shift, cz]
    render_view.CameraPosition = [cx, cy - vertical_shift, cz + max(x_span, y_span, 1.0)]


def compute_view_size(bounds: tuple[float, ...]) -> list[int]:
    """Choose an ultra-wide render size so the whole channel remains legible."""
    xmin, xmax, ymin, ymax, _zmin, _zmax = bounds
    x_span = max(xmax - xmin, 1e-9)
    y_span = max(ymax - ymin, 1e-9)

    channel_aspect = x_span / y_span
    desired_channel_height_fraction = 0.78
    target_view_aspect = max(1.0, desired_channel_height_fraction * channel_aspect)

    view_width = 16000
    view_height = max(760, round(view_width / target_view_aspect))
    return [view_width, view_height]


def enable_offscreen_rendering(render_view) -> None:
    """Force ParaView to render off-screen so batch screenshots never pop up a window."""
    for attr in ("UseOffscreenRendering", "UseOffscreenRenderingForScreenshots"):
        # Newer ParaView versions may remove these properties entirely and
        # raise from proxy attribute lookup. Treat them as optional.
        try:
            setattr(render_view, attr, 1)
        except Exception:
            pass


def autocrop_png(png_path: Path, tolerance: int = 2, padding: int = 6) -> None:
    """Crop away the uniform white border from a saved screenshot."""
    reader = vtk.vtkPNGReader()
    reader.SetFileName(str(png_path))
    reader.Update()
    image = reader.GetOutput()

    scalars = image.GetPointData().GetScalars()
    if scalars is None:
        return

    width, height, _ = image.GetDimensions()
    components = scalars.GetNumberOfComponents()
    pixels = vtk_to_numpy(scalars).reshape(height, width, components)
    rgb = pixels[..., :3].astype(np.int16)
    background = np.array([255, 255, 255], dtype=np.int16)
    mask = np.any(np.abs(rgb - background) > tolerance, axis=2)
    if components >= 4:
        mask |= pixels[..., 3] > tolerance

    if not mask.any():
        return

    ys, xs = np.where(mask)
    xmin = max(int(xs.min()) - padding, 0)
    xmax = min(int(xs.max()) + padding, width - 1)
    ymin = max(int(ys.min()) - padding, 0)
    ymax = min(int(ys.max()) + padding, height - 1)

    clip = vtk.vtkImageClip()
    clip.SetInputData(image)
    clip.SetOutputWholeExtent(xmin, xmax, ymin, ymax, 0, 0)
    clip.ClipDataOn()
    clip.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(png_path))
    writer.SetInputData(clip.GetOutput())
    writer.Write()


def field_display_name(field: str) -> str:
    """Human-readable field label for annotations."""
    return "Temperature" if field == "T" else field


def append_colorbar(
    png_path: Path,
    field: str,
    field_min: float = FIELD_MIN,
    field_max: float = FIELD_MAX,
) -> None:
    """Append a crisp matplotlib colorbar below the rendered geometry image."""
    geometry = Image.open(png_path).convert("RGB")
    width, _height = geometry.size

    dpi = 200
    bar_height_px = 420
    fig = plt.figure(figsize=(width / dpi, bar_height_px / dpi), dpi=dpi, facecolor="white")
    ax = fig.add_axes([0.10, 0.32, 0.80, 0.28])

    sm = plt.cm.ScalarMappable(norm=Normalize(vmin=field_min, vmax=field_max), cmap=rainbow_colormap())
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=ax, orientation="horizontal")
    cbar.set_label(field_display_name(field), fontsize=52, fontweight="bold", labelpad=18)
    cbar.ax.tick_params(labelsize=38, width=1.8, length=10, pad=8)
    cbar.outline.set_linewidth(1.8)
    ticks = np.linspace(field_min, field_max, 5)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{tick:.2f}" for tick in ticks])

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, facecolor="white")
    plt.close(fig)
    buffer.seek(0)

    colorbar = Image.open(buffer).convert("RGB")
    if colorbar.width != width:
        colorbar = colorbar.resize((width, colorbar.height), Image.Resampling.LANCZOS)

    gap = 24
    combined = Image.new("RGB", (width, geometry.height + gap + colorbar.height), "white")
    combined.paste(geometry, (0, 0))
    combined.paste(colorbar, (0, geometry.height + gap))
    combined.save(png_path)


def render_latest_field(case_foam: Path, output_png: Path, field: str) -> None:
    """Render the latest reconstructed field from an OpenFOAM case to PNG."""
    output_png.parent.mkdir(parents=True, exist_ok=True)

    reader = OpenFOAMReader(FileName=str(case_foam))
    reader.MeshRegions = ['internalMesh']
    reader.CellArrays = [field]

    latest = latest_time(reader)
    if latest is None:
        UpdatePipeline(proxy=reader)
    else:
        UpdatePipeline(time=latest, proxy=reader)

    reader_bounds = reader.GetDataInformation().GetBounds()
    xmin, xmax, ymin, ymax, zmin, zmax = reader_bounds
    slice_filter = Slice(Input=reader)
    slice_filter.SliceType = 'Plane'
    slice_filter.SliceOffsetValues = [0.0]
    slice_filter.SliceType.Origin = [
        0.5 * (xmin + xmax),
        0.5 * (ymin + ymax),
        0.5 * (zmin + zmax),
    ]
    slice_filter.SliceType.Normal = [0.0, 0.0, 1.0]

    if latest is None:
        UpdatePipeline(proxy=slice_filter)
    else:
        UpdatePipeline(time=latest, proxy=slice_filter)

    bounds = slice_filter.GetDataInformation().GetBounds()
    render_view = GetActiveViewOrCreate('RenderView')
    render_view.ViewSize = compute_view_size(bounds)
    render_view.Background = [1.0, 1.0, 1.0]
    render_view.OrientationAxesVisibility = 0
    enable_offscreen_rendering(render_view)

    display = Show(slice_filter, render_view)
    display.Representation = 'Surface'
    try:
        ColorBy(display, ('CELLS', field))
    except Exception:
        ColorBy(display, ('POINTS', field))

    lut = GetColorTransferFunction(field)
    pwf = GetOpacityTransferFunction(field)
    set_rainbow_transfer_function(lut, FIELD_MIN, FIELD_MAX)
    pwf.RescaleTransferFunction(FIELD_MIN, FIELD_MAX)
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(render_view, False)

    fit_xy_camera(render_view, bounds)
    Render()

    SaveScreenshot(str(output_png), render_view, ImageResolution=render_view.ViewSize)
    autocrop_png(output_png)
    append_colorbar(output_png, field)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-foam', required=True, help='Path to the OpenFOAM .foam file')
    parser.add_argument('--output', required=True, help='PNG output path')
    parser.add_argument('--field', default='T', help='Scalar field to render (default: T)')
    args = parser.parse_args()

    render_latest_field(Path(args.case_foam).resolve(), Path(args.output).resolve(), args.field)


if __name__ == '__main__':
    main()
