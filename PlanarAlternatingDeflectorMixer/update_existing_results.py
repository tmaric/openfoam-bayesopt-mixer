#!/usr/bin/env python3
"""
Migrate already-computed BO samples to the current portable study layout and
ensure the samplingPlane post-processing function object is available.

What this script does for each results/<sample_id>/:
  1. Renames legacy flow/scalar case directories and generated artifacts.
  2. Rewrites objectives.csv with portable sample-relative paths and correct
     kinematic-pressure units.
  3. Copies system/samplingPlane from the template mixing case.
  4. Appends the samplingPlane function-object block to system/controlDict
     (idempotent – skipped if already present).

Run from the PlanarAlternatingDeflectorMixer directory:
    python3 update_existing_results.py [--results-dir results]
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE_MIXING = SCRIPT_DIR / "ScalarTransportCase"

LEGACY_COLUMN_NAMES = {
    "pdrop_patch1_average_Pa": "pdrop_patch1_average_m2_s2",
    "pdrop_patch2_average_Pa": "pdrop_patch2_average_m2_s2",
    "pdrop_pressure_drop_Pa": "pdrop_pressure_drop_m2_s2",
}

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


def rename_if_present(source: Path, destination: Path) -> bool:
    """Rename source when it exists and destination does not."""
    if not source.exists():
        return False
    if destination.exists():
        print(f"  WARNING: keeping both {source.name} and {destination.name}")
        return False
    source.rename(destination)
    return True


def replace_in_file(path: Path, replacements: dict[str, str]) -> None:
    """Apply fixed-string replacements to a text file when it exists."""
    if not path.exists():
        return
    text = path.read_text()
    updated = text
    for old, new in replacements.items():
        updated = updated.replace(old, new)
    if updated != text:
        path.write_text(updated)


def migrate_case_layout(sample_dir: Path) -> None:
    """Move legacy result-case names to the current study vocabulary."""
    flow_case = sample_dir / "FlowCase"
    scalar_case = sample_dir / "ScalarTransportCase"

    if rename_if_present(sample_dir / "SplitAndRecombineHydro", flow_case):
        print(f"  {sample_dir.name}: renamed flow case directory")
    if rename_if_present(sample_dir / "SplitAndRecombineMixing", scalar_case):
        print(f"  {sample_dir.name}: renamed scalar case directory")

    if flow_case.is_dir():
        renames = (
            (flow_case / "SplitAndRecombineHydro.foam", flow_case / "FlowCase.foam"),
            (flow_case / "sar_mixer_cad.py", flow_case / "alternating_deflector_cad.py"),
            (flow_case / "sar_mixer_cad.yaml", flow_case / "alternating_deflector_cad.yaml"),
            (flow_case / "log.sar_mixer_cad", flow_case / "log.alternating_deflector_cad"),
            (
                flow_case / "constant" / "triSurface" / "sar_mixer.stl",
                flow_case / "constant" / "triSurface" / "alternating_deflector_mixer.stl",
            ),
        )
        for source, destination in renames:
            rename_if_present(source, destination)

        flow_replacements = {
            "sar_mixer_cad.py": "alternating_deflector_cad.py",
            "sar_mixer_cad.yaml": "alternating_deflector_cad.yaml",
            "sar_mixer.stl": "alternating_deflector_mixer.stl",
            "log.sar_mixer_cad": "log.alternating_deflector_cad",
            "SplitAndRecombineHydro.foam": "FlowCase.foam",
        }
        for path in (
            flow_case / "Allrun",
            flow_case / "Allclean",
            flow_case / "alternating_deflector_cad.py",
            flow_case / "system" / "meshDict",
        ):
            replace_in_file(path, flow_replacements)

    if scalar_case.is_dir():
        rename_if_present(
            scalar_case / "SplitAndRecombineMixing.foam",
            scalar_case / "ScalarTransportCase.foam",
        )
        replace_in_file(
            scalar_case / "Allrun",
            {
                "../SplitAndRecombineHydro": "../FlowCase",
                "SplitAndRecombineMixing.foam": "ScalarTransportCase.foam",
            },
        )


def migrate_objectives_csv(objectives_path: Path, sample_id: str) -> bool:
    """Normalize archived objective metadata without changing values."""
    if not objectives_path.exists():
        return False

    with objectives_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if not rows:
        return False

    changed = False
    migrated_fields = []
    for field in fieldnames:
        replacement = LEGACY_COLUMN_NAMES.get(field, field)
        migrated_fields.append(replacement)
        changed |= replacement != field

    migrated_rows = []
    for row in rows:
        migrated = {
            LEGACY_COLUMN_NAMES.get(key, key): value
            for key, value in row.items()
        }
        if migrated.get("results_dir") != sample_id:
            migrated["results_dir"] = sample_id
            changed = True

        # Reproject the archived physical geometry through the current latent
        # transform. This matters when mesh-safety constraints tighten while
        # retaining the already-computed objective values.
        try:
            from bayes_optimize_sequential import GEO_PARAM_NAMES, geo_to_bo

            geo = {
                name: float(migrated[f"geo_{name}"])
                for name in GEO_PARAM_NAMES
            }
            bo = geo_to_bo(geo)
        except (ImportError, KeyError, TypeError, ValueError):
            bo = None
        if bo is not None:
            for name, value in bo.items():
                field = f"bo_{name}"
                formatted = f"{value:.17g}"
                if field not in migrated_fields:
                    migrated_fields.append(field)
                    changed = True
                if migrated.get(field) != formatted:
                    migrated[field] = formatted
                    changed = True
        migrated_rows.append(migrated)

    if changed:
        with objectives_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=migrated_fields)
            writer.writeheader()
            writer.writerows(migrated_rows)

    return changed


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
        migrate_case_layout(sample_dir)
        if migrate_objectives_csv(sample_dir / "objectives.csv", sample_dir.name):
            print(f"  {sample_dir.name}: migrated objectives.csv")

        mixing_system = sample_dir / "ScalarTransportCase" / "system"
        if not mixing_system.is_dir():
            print(f"  SKIP {sample_dir.name}: no scalar-transport system directory")
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
