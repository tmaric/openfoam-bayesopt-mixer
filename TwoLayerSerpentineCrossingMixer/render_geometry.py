#!/usr/bin/env python3
"""Render the generated M10 STL without ParaView or a display server."""

from __future__ import annotations

import argparse
from pathlib import Path

import vtk


ROOT = Path(__file__).resolve().parent


def render(stl_path: Path, output_path: Path) -> None:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(stl_path))
    reader.MergingOn()
    reader.Update()
    if reader.GetOutput().GetNumberOfCells() == 0:
        raise RuntimeError(f"no triangles read from {stl_path}")

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.19, 0.58, 0.82)
    actor.GetProperty().SetSpecular(0.25)
    actor.GetProperty().SetSpecularPower(24.0)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.97, 0.98, 1.0)
    renderer.AddActor(actor)
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Azimuth(-32.0)
    camera.Elevation(28.0)
    camera.Zoom(1.25)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(1800, 700)
    window.AddRenderer(renderer)
    window.Render()

    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetScale(1)
    capture.SetInputBufferTypeToRGBA()
    capture.ReadFrontBufferOff()
    capture.Update()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()
    window.Finalize()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stl",
        type=Path,
        default=ROOT / "generated/reference/two_layer_serpentine_crossing_mixer.stl",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "generated/reference/geometry.png",
    )
    args = parser.parse_args()
    render(args.stl.resolve(), args.output.resolve())
    print(args.output.resolve())


if __name__ == "__main__":
    main()
