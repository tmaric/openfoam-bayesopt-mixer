#!/usr/bin/env python3
"""
Patch already-computed BO samples so that the samplingPlane post-processing
function object is available for every existing SplitAndRecombineMixing case.

What this script does for each results/<sample_id>/SplitAndRecombineMixing/:
  1. Copies system/samplingPlane from the template mixing case.
  2. Appends the samplingPlane function-object block to system/controlDict
     (idempotent – skipped if already present).

Run from the SplitAndRecombineMixer directory:
    python3 update_existing_results.py [--results-dir results]
"""

import argparse
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE_MIXING = SCRIPT_DIR / "SplitAndRecombineMixing"

SAMPLING_PLANE_BLOCK = """
    // Horizontal cut plane at z = DEPTH/2 = 5e-6 m.
    // Writes VTK surface with U and T fields.
    // Executed via: postProcess -func samplingPlane -latestTime
    samplingPlane
    {
        type            surfaces;
        libs            ("libsampling.so");

        writeControl    onEnd;

        surfaceFormat   vtk;
        interpolationScheme cellPoint;

        fields          (U T);

        surfaces
        (
            midPlane
            {
                type         cuttingPlane;
                planeType    pointAndNormal;
                pointAndNormalDict
                {
                    basePoint    (0 0 5e-6);
                    normalVector (0 0 1);
                }
                interpolate  true;
            }
        );
    }
"""


def patch_controldict(controldict_path: Path) -> bool:
    """
    Add samplingPlane block inside the functions{} section.
    Returns True if the file was modified, False if already patched.
    """
    text = controldict_path.read_text()
    if "samplingPlane" in text:
        return False  # already patched

    # Find the closing brace of the functions{} section.
    # Walk character-by-character from the 'functions' keyword.
    idx = text.find("functions")
    if idx == -1:
        print(f"  WARNING: no 'functions' entry in {controldict_path}", file=sys.stderr)
        return False

    start = text.find("{", idx)
    if start == -1:
        return False

    depth = 0
    pos = start
    close_pos = -1
    while pos < len(text):
        if text[pos] == "{":
            depth += 1
        elif text[pos] == "}":
            depth -= 1
            if depth == 0:
                close_pos = pos
                break
        pos += 1

    if close_pos == -1:
        print(f"  WARNING: could not find closing brace of functions in {controldict_path}",
              file=sys.stderr)
        return False

    # Insert the block just before the closing brace.
    new_text = text[:close_pos] + SAMPLING_PLANE_BLOCK + "\n" + text[close_pos:]
    controldict_path.write_text(new_text)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        default=str(SCRIPT_DIR / "results"),
        help="Root results directory (default: %(default)s)",
    )
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    if not results_root.is_dir():
        sys.exit(f"ERROR: results directory not found: {results_root}")

    template_samplingplane = TEMPLATE_MIXING / "system" / "samplingPlane"
    if not template_samplingplane.exists():
        sys.exit(f"ERROR: template samplingPlane not found: {template_samplingplane}")

    # Iterate over numeric sample directories.
    sample_dirs = sorted(
        d for d in results_root.iterdir()
        if d.is_dir() and d.name.isdigit()
    )

    if not sample_dirs:
        print(f"No sample directories found under {results_root}")
        return

    patched = 0
    for sample_dir in sample_dirs:
        mixing_system = sample_dir / "SplitAndRecombineMixing" / "system"
        if not mixing_system.is_dir():
            print(f"  SKIP {sample_dir.name}: no SplitAndRecombineMixing/system/")
            continue

        # 1. Copy system/samplingPlane
        dest = mixing_system / "samplingPlane"
        if not dest.exists():
            shutil.copy2(str(template_samplingplane), str(dest))
            print(f"  {sample_dir.name}: copied system/samplingPlane")
        else:
            print(f"  {sample_dir.name}: system/samplingPlane already present")

        # 2. Patch controlDict
        controldict = mixing_system / "controlDict"
        if not controldict.exists():
            print(f"  SKIP {sample_dir.name}: no system/controlDict")
            continue
        if patch_controldict(controldict):
            print(f"  {sample_dir.name}: patched system/controlDict")
            patched += 1
        else:
            print(f"  {sample_dir.name}: controlDict already has samplingPlane")

    print(f"\nDone. {patched}/{len(sample_dirs)} controlDicts patched.")


if __name__ == "__main__":
    main()
