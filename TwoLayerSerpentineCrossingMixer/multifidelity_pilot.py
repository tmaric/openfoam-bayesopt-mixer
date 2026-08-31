#!/usr/bin/env python3
"""Run and assess the bounded M10-inspired coarse/fine correlation pilot.

The pilot evaluates six paired designs drawn from the eventual 24-point Sobol
initialization. Evaluations are launched one at a time (q=1), while each
OpenFOAM case may use up to four MPI ranks. Completed pilot observations are
therefore reusable as campaign initialization data rather than throwaway
screening runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import warnings

from scipy.stats import qmc, spearmanr
import yaml

from multifidelity_design import (
    CAMPAIGN_PATH,
    ROOT,
    load_yaml,
    materialize_design,
    reference_parameters,
    validate_campaign,
)
from run_case import validate_boundary_topology, validate_mesh_quality


RUN_CASE = ROOT / "run_case.py"


def campaign_config() -> dict:
    config = load_yaml(CAMPAIGN_PATH)
    validate_campaign(config)
    pilot = config["pilot"]
    if int(pilot["paired_designs"]) != int(
        config["initialization"]["paired_fine_anchors"]
    ):
        raise ValueError("pilot pairs must equal the planned paired fine anchors")
    if int(pilot["paired_designs"]) > int(
        config["initialization"]["coarse_sobol_designs"]
    ):
        raise ValueError("pilot cannot select more pairs than the coarse Sobol pool")
    return config


def results_root(config: dict) -> Path:
    path = (ROOT / config["pilot"]["results_dir"]).resolve()
    path.relative_to((ROOT / "results").resolve())
    return path


def sobol_pool(config: dict) -> list[dict]:
    bounds = config["design"]["parameters"]
    names = list(bounds)
    count = int(config["initialization"]["coarse_sobol_designs"])
    seed = int(config["initialization"]["sobol_seed"])
    sampler = qmc.Sobol(d=len(names), scramble=True, seed=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        normalized = sampler.random(count)
    designs = []
    for index, point in enumerate(normalized):
        parameters = {}
        for coordinate, name in zip(point, names):
            lower = float(bounds[name]["lower"])
            upper = float(bounds[name]["upper"])
            parameters[name] = lower + float(coordinate) * (upper - lower)
        designs.append(
            {
                "sobol_index": index,
                "design_id": f"sobol_{index:03d}",
                "normalized": [float(value) for value in point],
                "parameters": parameters,
            }
        )
    return designs


def paired_anchor_indices(config: dict, pool: list[dict]) -> list[int]:
    bounds = config["design"]["parameters"]
    names = list(bounds)
    reference = reference_parameters()
    normalized_reference = []
    for name in names:
        lower = float(bounds[name]["lower"])
        upper = float(bounds[name]["upper"])
        normalized_reference.append((reference[name] - lower) / (upper - lower))

    def squared_distance(left: list[float], right: list[float]) -> float:
        return sum((a - b) ** 2 for a, b in zip(left, right))

    first = min(
        range(len(pool)),
        key=lambda index: (
            squared_distance(pool[index]["normalized"], normalized_reference),
            index,
        ),
    )
    selected = [first]
    target = int(config["pilot"]["paired_designs"])
    while len(selected) < target:
        candidates = [index for index in range(len(pool)) if index not in selected]
        next_index = max(
            candidates,
            key=lambda index: (
                min(
                    squared_distance(
                        pool[index]["normalized"], pool[chosen]["normalized"]
                    )
                    for chosen in selected
                ),
                -index,
            ),
        )
        selected.append(next_index)
    return selected


def manifest(config: dict) -> dict:
    pool = sobol_pool(config)
    anchors = paired_anchor_indices(config, pool)
    return {
        "schema_version": 1,
        "campaign": config["campaign"],
        "scientific_framing": config["scientific_framing"],
        "operating_reynolds_number": float(
            config["design"]["operating_reynolds_number"]
        ),
        "sobol_seed": int(config["initialization"]["sobol_seed"]),
        "coarse_sobol_pool_size": len(pool),
        "selection": config["pilot"]["selection"],
        "paired_anchor_indices": anchors,
        "paired_designs": [pool[index] for index in anchors],
    }


def ensure_manifest(config: dict) -> tuple[Path, dict]:
    root = results_root(config)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "pilot_manifest.json"
    expected = manifest(config)
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != expected:
            raise RuntimeError("existing pilot manifest differs from campaign configuration")
    else:
        path.write_text(json.dumps(expected, indent=2) + "\n", encoding="utf-8")
    return root, expected


def evaluation_sequence(config: dict, data: dict) -> list[tuple[dict, str]]:
    if int(config["resources"]["q"]) != 1:
        raise RuntimeError("pilot execution requires q=1")
    return [
        (design, fidelity)
        for design in data["paired_designs"]
        for fidelity in ("coarse", "fine")
    ]


def fidelity_cell_sizes(
    config: dict,
    coarse_cell_size_m: float | None = None,
    fine_cell_size_m: float | None = None,
) -> dict[str, float]:
    sizes = {
        "coarse": float(config["fidelities"]["coarse"]["nominal_cell_size_m"]),
        "fine": float(config["fidelities"]["fine"]["nominal_cell_size_m"]),
    }
    if coarse_cell_size_m is not None:
        sizes["coarse"] = float(coarse_cell_size_m)
    if fine_cell_size_m is not None:
        sizes["fine"] = float(fine_cell_size_m)
    if sizes["coarse"] <= sizes["fine"]:
        raise ValueError("coarse spacing must be larger than fine spacing")
    if any(value <= 0.0 for value in sizes.values()):
        raise ValueError("preflight spacings must be positive")
    return sizes


def preflight_root(config: dict, cell_sizes: dict[str, float] | None = None) -> Path:
    sizes = cell_sizes or fidelity_cell_sizes(config)
    coarse_um = round(sizes["coarse"] * 1.0e6)
    fine_um = round(sizes["fine"] * 1.0e6)
    path = (
        ROOT
        / "results"
        / "mesh_qualification"
        / f"snappy_v2606_pilot_preflight_{coarse_um}um_{fine_um}um"
    ).resolve()
    path.relative_to((ROOT / "results").resolve())
    return path


def revalidate_existing_preflight(target: Path, fidelity: str) -> dict:
    """Apply the current fixed policy to an already generated preflight mesh."""

    check_log = target / "FlowCase" / "log.checkMesh"
    boundary = target / "FlowCase" / "constant" / "polyMesh" / "boundary"
    validation = validate_mesh_quality(check_log, fidelity)
    validation["boundary_topology"] = validate_boundary_topology(boundary)
    (target / "mesh_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    return validation


def run_preflight(
    coarse_cell_size_m: float | None = None,
    fine_cell_size_m: float | None = None,
) -> dict:
    config = campaign_config()
    pilot_root, data = ensure_manifest(config)
    cell_sizes = fidelity_cell_sizes(
        config, coarse_cell_size_m, fine_cell_size_m
    )
    root = preflight_root(config, cell_sizes)
    root.mkdir(parents=True, exist_ok=True)
    evaluations = []
    all_passed = True
    for design, fidelity in evaluation_sequence(config, data):
        target = root / design["design_id"] / fidelity
        validation_path = target / "mesh_validation.json"
        if target.exists():
            try:
                validation = revalidate_existing_preflight(target, fidelity)
            except (FileNotFoundError, RuntimeError):
                evaluations.append(
                    {
                        "design_id": design["design_id"],
                        "fidelity": fidelity,
                        "status": "failed_or_incomplete",
                    }
                )
                all_passed = False
                continue
            evaluations.append(
                {
                    "design_id": design["design_id"],
                    "fidelity": fidelity,
                    "status": "passed",
                    "mesh_validation": validation,
                }
            )
            continue
        geometry = write_geometry_config(config, pilot_root, design)
        command = [
            sys.executable,
            str(RUN_CASE),
            "--protocol",
            "original",
            "--fidelity",
            fidelity,
            "--reynolds",
            str(config["design"]["operating_reynolds_number"]),
            "--cores",
            str(config["resources"]["cores"]),
            "--cell-size-m",
            str(cell_sizes[fidelity]),
            "--geometry-config",
            str(geometry),
            "--results-dir",
            str(target),
            "--mesh-only",
        ]
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        print(f"preflighting {design['design_id']} {fidelity}", flush=True)
        require_no_other_mpi_launcher()
        completed = subprocess.run(command, cwd=ROOT, env=env, check=False)
        if completed.returncode == 0 and validation_path.exists():
            status = "passed"
            validation = json.loads(validation_path.read_text(encoding="utf-8"))
        else:
            status = "failed"
            validation = None
            all_passed = False
        row = {
            "design_id": design["design_id"],
            "fidelity": fidelity,
            "status": status,
        }
        if validation is not None:
            row["mesh_validation"] = validation
        evaluations.append(row)
    summary = {
        "schema_version": 1,
        "coarse_cell_size_m": cell_sizes["coarse"],
        "fine_cell_size_m": cell_sizes["fine"],
        "all_meshes_passed": all_passed,
        "evaluations": evaluations,
    }
    (root / "preflight_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return summary


def require_passed_preflight(config: dict) -> None:
    path = preflight_root(config) / "preflight_summary.json"
    if not path.exists():
        raise RuntimeError("run the mesh-only preflight before pilot transport")
    summary = json.loads(path.read_text(encoding="utf-8"))
    if not bool(summary["all_meshes_passed"]):
        raise RuntimeError("mesh-only preflight has not passed for all pilot designs")


def objective_path(root: Path, design: dict, fidelity: str) -> Path:
    return root / design["design_id"] / fidelity / "objectives.json"


def case_status(root: Path, design: dict, fidelity: str) -> str:
    case = root / design["design_id"] / fidelity
    if (case / "objectives.json").exists():
        return "complete"
    if case.exists():
        return "incomplete"
    return "pending"


def write_geometry_config(config: dict, root: Path, design: dict) -> Path:
    path = root / design["design_id"] / "geometry.yaml"
    materialized = materialize_design(design["parameters"], "coarse")
    geometry_config = {
        "schema_version": 1,
        "scientific_framing": config["scientific_framing"],
        "bo_parameters": design["parameters"],
        "geometry": materialized["geometry"],
        "output": materialized["output"],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = yaml.safe_dump(geometry_config, sort_keys=False)
    if path.exists() and path.read_text(encoding="utf-8") != serialized:
        raise RuntimeError(f"existing geometry configuration differs: {path}")
    path.write_text(serialized, encoding="utf-8")
    return path


def acquire_lock(root: Path) -> tuple[int, Path]:
    lock = root / ".pilot.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        owner = lock.read_text(encoding="utf-8", errors="replace").strip()
        raise RuntimeError(f"pilot is already locked by {owner or 'an unknown process'}") from exc
    os.write(descriptor, f"pid={os.getpid()}\n".encode())
    return descriptor, lock


def require_no_other_mpi_launcher() -> None:
    """Avoid adding a four-rank case while another MPI job is active."""

    completed = subprocess.run(
        ["pgrep", "-a", "-x", "mpirun"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    if completed.returncode == 1:
        return
    if completed.returncode != 0:
        raise RuntimeError("could not inspect active MPI launchers with pgrep")
    launchers = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if launchers:
        raise RuntimeError(
            "refusing to launch another OpenFOAM case while MPI is active: "
            + "; ".join(launchers)
        )


def run_next(max_new_evaluations: int) -> None:
    if max_new_evaluations < 1:
        raise ValueError("max_new_evaluations must be positive")
    config = campaign_config()
    require_passed_preflight(config)
    root, data = ensure_manifest(config)
    descriptor, lock = acquire_lock(root)
    completed_now = 0
    try:
        for design, fidelity in evaluation_sequence(config, data):
            status = case_status(root, design, fidelity)
            if status == "complete":
                continue
            if status == "incomplete":
                raise RuntimeError(
                    f"incomplete case requires inspection before retry: "
                    f"{root / design['design_id'] / fidelity}"
                )
            geometry = write_geometry_config(config, root, design)
            result = root / design["design_id"] / fidelity
            command = [
                sys.executable,
                str(RUN_CASE),
                "--protocol",
                "original",
                "--fidelity",
                fidelity,
                "--reynolds",
                str(config["design"]["operating_reynolds_number"]),
                "--cores",
                str(config["resources"]["cores"]),
                "--geometry-config",
                str(geometry),
                "--results-dir",
                str(result),
            ]
            env = dict(os.environ)
            env["OMP_NUM_THREADS"] = "1"
            env["MKL_NUM_THREADS"] = "1"
            env["OPENBLAS_NUM_THREADS"] = "1"
            print(
                f"running {design['design_id']} {fidelity} "
                f"on {config['resources']['cores']} MPI ranks",
                flush=True,
            )
            require_no_other_mpi_launcher()
            subprocess.run(command, cwd=ROOT, env=env, check=True)
            completed_now += 1
            if completed_now >= max_new_evaluations:
                break
    finally:
        os.close(descriptor)
        lock.unlink(missing_ok=True)
    print(f"completed {completed_now} new sequential pilot evaluation(s)")


def status_summary(config: dict, root: Path, data: dict) -> dict:
    counts = {"complete": 0, "pending": 0, "incomplete": 0}
    evaluations = []
    for design, fidelity in evaluation_sequence(config, data):
        status = case_status(root, design, fidelity)
        counts[status] += 1
        evaluations.append(
            {
                "design_id": design["design_id"],
                "fidelity": fidelity,
                "status": status,
            }
        )
    return {"counts": counts, "evaluations": evaluations}


def paired_objectives(config: dict, root: Path, data: dict) -> list[dict]:
    pairs = []
    for design in data["paired_designs"]:
        row = {"design_id": design["design_id"], **design["parameters"]}
        for fidelity in ("coarse", "fine"):
            path = objective_path(root, design, fidelity)
            if not path.exists():
                raise RuntimeError(f"pilot is incomplete: {path}")
            objectives = json.loads(path.read_text(encoding="utf-8"))
            row[f"{fidelity}_pressure_ratio_to_straight"] = float(
                objectives["pressure_ratio_to_straight"]
            )
            row[f"{fidelity}_flux_weighted_intensity_of_segregation"] = float(
                objectives["flux_weighted_intensity_of_segregation"]
            )
            row[f"{fidelity}_flux_weighted_mixing_index"] = float(
                objectives["flux_weighted_mixing_index"]
            )
        pairs.append(row)
    return pairs


def pearson(left: list[float], right: list[float]) -> float:
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum(
        (a - left_mean) * (b - right_mean) for a, b in zip(left, right)
    )
    left_norm = math.sqrt(sum((value - left_mean) ** 2 for value in left))
    right_norm = math.sqrt(sum((value - right_mean) ** 2 for value in right))
    if left_norm == 0.0 or right_norm == 0.0:
        raise RuntimeError("correlation is undefined for a constant objective")
    return numerator / (left_norm * right_norm)


def objective_diagnostic(
    pairs: list[dict], objective: str, transform: str, minimum_spearman: float
) -> dict:
    coarse = [float(row[f"coarse_{objective}"]) for row in pairs]
    fine = [float(row[f"fine_{objective}"]) for row in pairs]
    if transform == "ln":
        coarse_transformed = [math.log(value) for value in coarse]
        fine_transformed = [math.log(value) for value in fine]
        inverse_mean_bias_ratio = math.exp
    elif transform == "log10":
        coarse_transformed = [math.log10(value) for value in coarse]
        fine_transformed = [math.log10(value) for value in fine]
        inverse_mean_bias_ratio = lambda value: 10.0**value
    else:
        raise ValueError(f"unknown transform {transform}")
    rank = float(spearmanr(coarse_transformed, fine_transformed).statistic)
    differences = [
        fine_value - coarse_value
        for coarse_value, fine_value in zip(coarse_transformed, fine_transformed)
    ]
    mean_bias = sum(differences) / len(differences)
    bias_std = math.sqrt(
        sum((value - mean_bias) ** 2 for value in differences)
        / (len(differences) - 1)
    )
    return {
        "transform": transform,
        "spearman_rank_correlation": rank,
        "pearson_correlation_transformed": pearson(
            coarse_transformed, fine_transformed
        ),
        "mean_fine_minus_coarse_transformed": mean_bias,
        "standard_deviation_fine_minus_coarse_transformed": bias_std,
        "geometric_mean_fine_to_coarse_ratio": inverse_mean_bias_ratio(mean_bias),
        "minimum_required_spearman": minimum_spearman,
        "gate_passed": rank >= minimum_spearman,
    }


def summarize() -> dict:
    config = campaign_config()
    root, data = ensure_manifest(config)
    pairs = paired_objectives(config, root, data)
    threshold = float(config["pilot"]["minimum_spearman_rank_correlation"])
    diagnostics = {
        "pressure_ratio_to_straight": objective_diagnostic(
            pairs, "pressure_ratio_to_straight", "ln", threshold
        ),
        "flux_weighted_intensity_of_segregation": objective_diagnostic(
            pairs,
            "flux_weighted_intensity_of_segregation",
            "log10",
            threshold,
        ),
    }
    summary = {
        "schema_version": 1,
        "campaign": config["campaign"],
        "scientific_framing": config["scientific_framing"],
        "paired_design_count": len(pairs),
        "operating_reynolds_number": float(
            config["design"]["operating_reynolds_number"]
        ),
        "diagnostics": diagnostics,
        "multifidelity_gate_passed": all(
            diagnostic["gate_passed"] for diagnostic in diagnostics.values()
        ),
    }
    (root / "pilot_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    with (root / "paired_objectives.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pairs[0]))
        writer.writeheader()
        writer.writerows(pairs)
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("designs")
    subparsers.add_parser("status")
    preflight_parser = subparsers.add_parser("preflight")
    preflight_parser.add_argument("--coarse-cell-size-m", type=float)
    preflight_parser.add_argument("--fine-cell-size-m", type=float)
    next_parser = subparsers.add_parser("next")
    next_parser.add_argument("--max-new-evaluations", type=int, default=1)
    subparsers.add_parser("summarize")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = campaign_config()
    root, data = ensure_manifest(config)
    if args.command == "designs":
        print(json.dumps(data, indent=2))
    elif args.command == "status":
        print(json.dumps(status_summary(config, root, data), indent=2))
    elif args.command == "preflight":
        run_preflight(args.coarse_cell_size_m, args.fine_cell_size_m)
    elif args.command == "next":
        run_next(args.max_new_evaluations)
    elif args.command == "summarize":
        summarize()


if __name__ == "__main__":
    main()
