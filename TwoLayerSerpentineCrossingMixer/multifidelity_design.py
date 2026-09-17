#!/usr/bin/env python3
"""Validate and materialize M10-inspired designs at an explicit fidelity.

The fidelity coordinate is part of the statistical input, not a label applied
after a conventional BO run. Both levels evaluate the same physical design and
objectives; they differ only in discretization. Fine (s=1) is the target
fidelity, while coarse (s=0) learns correlated low-cost information.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parent
CAMPAIGN_PATH = ROOT / "bayes_optimize_multifidelity.yaml"
REFERENCE_PATH = ROOT / "FlowCase" / "two_layer_serpentine_crossing_cad.yaml"


def load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def validate_campaign(config: dict) -> None:
    framing = config["scientific_framing"]
    if framing["classification"] != "m10_inspired_reconstruction":
        raise ValueError("campaign must record the M10-inspired reconstruction framing")
    if bool(framing["exact_reproduction_claim"]):
        raise ValueError("the published geometry does not support an exact reproduction claim")

    resources = config["resources"]
    if int(resources["q"]) != 1:
        raise ValueError("M10 BO must remain strictly sequential (q=1)")
    if not 1 <= int(resources["cores"]) <= 4:
        raise ValueError("cores must remain in [1, 4]")
    if int(resources["torch_threads"]) != 1:
        raise ValueError("torch_threads must remain 1")

    fidelities = config["fidelities"]
    if set(fidelities) != {"coarse", "fine"}:
        raise ValueError("exactly coarse and fine fidelities are required")
    if float(fidelities["coarse"]["coordinate"]) != 0.0:
        raise ValueError("coarse fidelity coordinate must be 0")
    if float(fidelities["fine"]["coordinate"]) != 1.0:
        raise ValueError("fine fidelity coordinate must be 1")
    if float(fidelities["fine"]["nominal_relative_cost"]) <= float(
        fidelities["coarse"]["nominal_relative_cost"]
    ):
        raise ValueError("fine fidelity must cost more than coarse fidelity")
    if config["model"]["class"] != "SingleTaskMultiFidelityGP":
        raise ValueError("the campaign must use an explicit multifidelity GP")
    if float(config["model"]["target_fidelity"]) != 1.0:
        raise ValueError("publication target fidelity must be fine (s=1)")
    if not bool(config["model"]["cost_aware"]):
        raise ValueError("multifidelity acquisition must be cost-aware")
    if float(config["design"]["operating_reynolds_number"]) <= 0.0:
        raise ValueError("the fixed BO operating Reynolds number must be positive")
    transforms = config["model"]["outcome_transforms"]
    if transforms["pressure_ratio_to_straight"] != "natural_log":
        raise ValueError("pressure-ratio observations must use a log transform")
    intensity_transform = transforms["flux_weighted_intensity_of_segregation"]
    if intensity_transform["transform"] != "log10" or float(intensity_transform["floor"]) <= 0.0:
        raise ValueError("segregation intensity requires a positive-floor log10 transform")

    for name, bounds in config["design"]["parameters"].items():
        if float(bounds["lower"]) >= float(bounds["upper"]):
            raise ValueError(f"invalid bounds for {name}")

    pilot = config["pilot"]
    if int(pilot["paired_designs"]) != int(
        config["initialization"]["paired_fine_anchors"]
    ):
        raise ValueError("pilot pairs must equal paired fine initialization anchors")
    rank_threshold = float(pilot["minimum_spearman_rank_correlation"])
    if not 0.0 < rank_threshold <= 1.0:
        raise ValueError("pilot Spearman threshold must lie in (0, 1]")


def materialize_design(parameters: dict[str, float], fidelity: str) -> dict:
    campaign = load_yaml(CAMPAIGN_PATH)
    reference = load_yaml(REFERENCE_PATH)
    validate_campaign(campaign)
    if fidelity not in campaign["fidelities"]:
        raise ValueError(f"unknown fidelity {fidelity!r}")

    bounds = campaign["design"]["parameters"]
    missing = set(bounds) - set(parameters)
    extra = set(parameters) - set(bounds)
    if missing or extra:
        raise ValueError(f"parameter mismatch: missing={sorted(missing)}, extra={sorted(extra)}")
    values = {name: float(parameters[name]) for name in bounds}
    for name, value in values.items():
        lower = float(bounds[name]["lower"])
        upper = float(bounds[name]["upper"])
        if not lower <= value <= upper:
            raise ValueError(f"{name}={value} lies outside [{lower}, {upper}]")

    P = float(reference["geometry"]["clear_pitch_P"])
    geometry = dict(reference["geometry"])
    geometry.update(
        {
            "main_span_H": values["H_over_P"] * P,
            "diagonal_width_w": values["w_over_P"] * P,
            "single_layer_depth_d": 0.5 * values["D_over_P"] * P,
            "vertical_segment_width_b": values["b_over_P"] * P,
            "diagonal_end_inset_over_w": values["diagonal_end_inset_over_w"],
            "crossing_phase_over_b": values["crossing_phase_over_b"],
            "number_of_units": int(campaign["design"]["fixed_number_of_units_during_bo"]),
        }
    )
    fidelity_cfg = campaign["fidelities"][fidelity]
    minimum_feature = min(
        geometry["diagonal_width_w"],
        geometry["vertical_segment_width_b"],
        geometry["single_layer_depth_d"],
    ) * 1.0e-3
    cell_size = float(fidelity_cfg["nominal_cell_size_m"])
    required_cells = int(campaign["constraints"][f"minimum_feature_cells_{fidelity}"])
    if minimum_feature / cell_size < required_cells:
        raise ValueError(
            f"design has only {minimum_feature / cell_size:.2f} {fidelity} cells "
            f"across its smallest feature; {required_cells} are required"
        )

    return {
        "schema_version": 1,
        "scientific_framing": campaign["scientific_framing"],
        "bo_parameters": values,
        "fidelity": {
            "name": fidelity,
            "coordinate": float(fidelity_cfg["coordinate"]),
            "cell_size_m": cell_size,
            "nominal_relative_cost": float(fidelity_cfg["nominal_relative_cost"]),
        },
        "geometry": geometry,
        "output": reference["output"],
    }


def reference_parameters() -> dict[str, float]:
    reference = load_yaml(REFERENCE_PATH)["geometry"]
    P = float(reference["clear_pitch_P"])
    return {
        "H_over_P": float(reference["main_span_H"]) / P,
        "w_over_P": float(reference["diagonal_width_w"]) / P,
        "D_over_P": 2.0 * float(reference["single_layer_depth_d"]) / P,
        "b_over_P": float(reference["vertical_segment_width_b"]) / P,
        "diagonal_end_inset_over_w": float(reference["diagonal_end_inset_over_w"]),
        "crossing_phase_over_b": float(reference["crossing_phase_over_b"]),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("validate-config")
    reference = subparsers.add_parser("reference")
    reference.add_argument("--fidelity", choices=("coarse", "fine"), required=True)
    reference.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    campaign = load_yaml(CAMPAIGN_PATH)
    validate_campaign(campaign)
    if args.command == "validate-config":
        print("M10-inspired multifidelity campaign configuration is valid")
        return
    design = materialize_design(reference_parameters(), args.fidelity)
    serialized = yaml.safe_dump(design, sort_keys=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
        print(json.dumps({"output": str(args.output), "fidelity": args.fidelity}))
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
