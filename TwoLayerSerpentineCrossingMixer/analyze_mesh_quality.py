#!/usr/bin/env python3
"""Summarize OpenFOAM ``checkMesh -writeSets vtk`` cell-set output.

This keeps mesh-quality diagnostics headless and reproducible.  It reads the
VTP written by OpenFOAM and reports the spatially disconnected bad-cell
clusters in millimetres, without requiring ParaView.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import vtk


def _bounds_mm(bounds: tuple[float, ...]) -> list[float]:
    return [round(value * 1.0e3, 6) for value in bounds]


def summarize(vtp_path: Path) -> dict[str, object]:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(vtp_path))
    reader.Update()
    surface = reader.GetOutput()

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.ColorRegionsOn()
    connectivity.Update()

    regions: list[dict[str, object]] = []
    for region_id in range(connectivity.GetNumberOfExtractedRegions()):
        extractor = vtk.vtkPolyDataConnectivityFilter()
        extractor.SetInputData(surface)
        extractor.SetExtractionModeToSpecifiedRegions()
        extractor.AddSpecifiedRegion(region_id)
        extractor.Update()
        # The connectivity filter retains points that are not referenced by
        # the selected region.  Cleaning removes those unused points so the
        # reported bounds describe the selected bad-cell cluster only.
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputConnection(extractor.GetOutputPort())
        cleaner.Update()
        region = cleaner.GetOutput()
        bounds = region.GetBounds()
        regions.append(
            {
                "region": region_id,
                "center_mm": [
                    round(0.5 * (bounds[0] + bounds[1]) * 1.0e3, 6),
                    round(0.5 * (bounds[2] + bounds[3]) * 1.0e3, 6),
                    round(0.5 * (bounds[4] + bounds[5]) * 1.0e3, 6),
                ],
                "bounds_mm": _bounds_mm(bounds),
                "surface_cells": region.GetNumberOfCells(),
            }
        )

    regions.sort(key=lambda item: item["center_mm"])
    return {
        "source": str(vtp_path.resolve()),
        "points": surface.GetNumberOfPoints(),
        "surface_cells": surface.GetNumberOfCells(),
        "bounds_mm": _bounds_mm(surface.GetBounds()),
        "connected_regions": len(regions),
        "regions": regions,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtp", type=Path, help="VTP file written by checkMesh")
    parser.add_argument("--output", type=Path, help="Optional JSON output path")
    args = parser.parse_args()

    report = summarize(args.vtp)
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
