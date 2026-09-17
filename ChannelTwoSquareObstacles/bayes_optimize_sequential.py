#!/usr/bin/env python3
"""Sequential single-objective BO for the TN-05 two-square channel benchmark."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch
import yaml

torch.set_default_dtype(torch.float64)

CASE_ROOT = Path(__file__).resolve().parent
HYDRO_ROOT = CASE_ROOT / "ChannelTwoSquareObstaclesHydro"
if str(HYDRO_ROOT) not in sys.path:
    sys.path.insert(0, str(HYDRO_ROOT))

from channel_two_square_obstacles_common import CFG_NAME as HYDRO_CFG_NAME  # noqa: E402
from channel_two_square_obstacles_common import compute_resolved_geometry  # noqa: E402


CFG_PATH = CASE_ROOT / "bayes_optimize_sequential.yaml"
with open(CFG_PATH, encoding="utf-8") as handle:
    CFG = yaml.safe_load(handle)

CAD_MODE = str(CFG.get("cad_mode", "constrained")).strip().lower()
if CAD_MODE not in {"constrained", "unconstrained"}:
    raise ValueError(
        f"Unsupported cad_mode '{CAD_MODE}' in {CFG_PATH}. "
        "Expected 'constrained' or 'unconstrained'."
    )

MODE_CFG = CFG["bo_parameters"][CAD_MODE]
BO_PARAM_NAMES = list(MODE_CFG.keys())
N_INIT = int(CFG["n_init"])
N_BO = int(CFG["n_bo"])
CORES = int(CFG["cores"])

RESULTS_DIR = CASE_ROOT / "results"
SNAKEFILE = CASE_ROOT / "Snakefile"
TEMPLATE_YAML = HYDRO_ROOT / HYDRO_CFG_NAME
MODEL_PATH = CASE_ROOT / "ChannelTwoSquareObstacles.pt"
OBJ_PDROP = "pdrop_pressure_drop_Pa"

with open(TEMPLATE_YAML, encoding="utf-8") as handle:
    TEMPLATE_GEO = yaml.safe_load(handle)

BOUNDS = torch.tensor(
    [
        [float(MODE_CFG[name]["lower"]) for name in BO_PARAM_NAMES],
        [float(MODE_CFG[name]["upper"]) for name in BO_PARAM_NAMES],
    ]
)


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def snakemake_bin() -> str:
    exe = shutil.which("snakemake")
    if exe is None:
        raise FileNotFoundError(
            "snakemake executable not found in PATH. "
            "Activate the environment that provides Snakemake and retry."
        )
    return exe


def _base_geometry_for_mode(bo_params: dict) -> dict:
    raw = dict(TEMPLATE_GEO)
    raw["cad_mode"] = CAD_MODE
    raw["a"] = float(bo_params["a"])
    if CAD_MODE == "constrained":
        raw["d_ratio"] = _clamp01(bo_params["d_ratio"])
    else:
        raw["d"] = float(bo_params["d"])
    return raw


def bo_to_geo(bo_params: dict) -> dict:
    return compute_resolved_geometry(_base_geometry_for_mode(bo_params), CAD_MODE)


def geo_to_bo(geo_params: dict) -> dict:
    a = float(geo_params["a"])
    d = float(geo_params["d"])
    if CAD_MODE == "constrained":
        spacing_margin = float(TEMPLATE_GEO["spacing_margin_H"])
        d_max = float(TEMPLATE_GEO["d_max_design"])
        denom = max(d_max - a - spacing_margin, 1.0e-12)
        d_ratio = _clamp01((d - a - spacing_margin) / denom)
        return {"a": a, "d_ratio": d_ratio}
    return {"a": a, "d": d}


def plot_objective_history(
    sample_ids: list[int],
    pdrop_values: list[float],
    phases: list[str],
    title: str,
) -> None:
    if not sample_ids:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    x = np.asarray(sample_ids, dtype=int)
    y = np.asarray(pdrop_values, dtype=float)
    colors = np.asarray(
        ["steelblue" if phase == "sobol" else "darkorange" for phase in phases],
        dtype=object,
    )
    best_so_far = np.maximum.accumulate(y)

    fig, ax = plt.subplots(figsize=(7, 5))
    sobol_mask = colors == "steelblue"
    bo_mask = colors == "darkorange"
    if sobol_mask.any():
        ax.scatter(x[sobol_mask], y[sobol_mask], c="steelblue", s=45, alpha=0.8)
    if bo_mask.any():
        ax.scatter(x[bo_mask], y[bo_mask], c="darkorange", s=45, alpha=0.8)

    ax.plot(x, best_so_far, color="black", lw=1.5, zorder=4)
    ax.scatter(x[-1], y[-1], c="red", s=100, marker="D", zorder=5)

    ax.set_xlabel("Sample index")
    ax.set_ylabel(r"Pressure drop  $J_\mathrm{dp}$  (Pa)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=8, label="Sobol init"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="darkorange", markersize=8, label="BO suggested"),
            Line2D([0], [0], color="black", lw=1.5, label="Best so far"),
            Line2D([0], [0], marker="D", color="w", markerfacecolor="red", markersize=8, label="Current"),
        ],
        loc="best",
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "pressure_drop_history.png", dpi=150)
    plt.close(fig)
    print(f"[bo] plot saved -> {RESULTS_DIR / 'pressure_drop_history.png'}")


def aggregate_all() -> None:
    rows = []
    fieldnames = []
    seen = set()
    for sample_dir in sorted(RESULTS_DIR.iterdir()) if RESULTS_DIR.exists() else []:
        if not (sample_dir.is_dir() and sample_dir.name.isdigit()):
            continue
        obj_csv = sample_dir / "objectives.csv"
        if not obj_csv.exists():
            continue
        with open(obj_csv, newline="") as handle:
            reader = csv.DictReader(handle)
            for field in reader.fieldnames or []:
                if field not in seen:
                    fieldnames.append(field)
                    seen.add(field)
            rows.extend(reader)
    if not rows:
        return
    out = RESULTS_DIR / "all_samples.csv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[bo] all_samples.csv -> {len(rows)} sample(s)  ({out})")


def next_sample_id() -> str:
    ids = (
        [int(sample_dir.name) for sample_dir in RESULTS_DIR.iterdir() if sample_dir.is_dir() and sample_dir.name.isdigit()]
        if RESULTS_DIR.exists()
        else []
    )
    return f"{(max(ids) + 1 if ids else 0):05d}"


def write_sample_yaml(sample_dir: Path, bo_params: dict) -> dict:
    geo = bo_to_geo(bo_params)
    raw = _base_geometry_for_mode(bo_params)
    raw["d"] = float(geo["d"])
    raw["cad_mode"] = CAD_MODE

    sample_dir.mkdir(parents=True, exist_ok=True)
    with open(sample_dir / HYDRO_CFG_NAME, "w", encoding="utf-8") as handle:
        yaml.safe_dump(raw, handle, default_flow_style=False, sort_keys=False)
    return geo


def _annotation_fields(bo_params: dict, geo_params: dict, phase: str) -> dict:
    return {
        "phase": phase,
        "cad_mode": CAD_MODE,
        **{f"bo_{key}": bo_params[key] for key in BO_PARAM_NAMES},
        **{f"geo_{key}": geo_params[key] for key in geo_params},
    }


def _annotate_objectives_csv(sample_dir: Path, bo_params: dict, geo_params: dict, phase: str) -> None:
    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        return

    with open(obj_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    extra = _annotation_fields(bo_params, geo_params, phase)
    extra["feasible"] = 1
    for key in extra:
        if key not in fieldnames:
            fieldnames.append(key)
    for row in rows:
        row.update(extra)

    with open(obj_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_failed_objectives(
    sample_dir: Path,
    bo_params: dict,
    geo_params: dict,
    phase: str,
    failure_stage: str,
    failure_reason: str,
) -> None:
    row = {
        "sample_id": sample_dir.name,
        "results_dir": str(sample_dir.resolve()),
        "phase": phase,
        "cad_mode": CAD_MODE,
        "feasible": 0,
        "failure_stage": failure_stage,
        "failure_reason": failure_reason,
        **{f"bo_{key}": bo_params[key] for key in BO_PARAM_NAMES},
        **{f"geo_{key}": geo_params[key] for key in geo_params},
        OBJ_PDROP: "",
    }
    out = sample_dir / "objectives.csv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    print(f"[bo] failed objectives.csv written -> {out}")


def read_objective(sample_dir: Path) -> float | None:
    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        return None
    with open(obj_csv, newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                return float(row[OBJ_PDROP])
            except (KeyError, TypeError, ValueError):
                return None
    return None


def fit_model(X: torch.Tensor, Y: torch.Tensor, warm_start: dict | None = None):
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.mlls import ExactMarginalLogLikelihood

    model = SingleTaskGP(
        X,
        Y,
        input_transform=Normalize(d=X.shape[-1]),
        outcome_transform=Standardize(m=Y.shape[-1]),
    )
    if warm_start is not None:
        try:
            with torch.no_grad():
                for name, param in model.named_parameters():
                    if name in warm_start and warm_start[name].shape == param.shape:
                        param.data.copy_(warm_start[name])
        except Exception as exc:
            print(
                f"[bo] WARNING: warm start failed ({exc}) - using default init",
                file=sys.stderr,
            )
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))
    return model


def save_model(model) -> None:
    checkpoint = {name: value.detach().clone() for name, value in model.named_parameters()}
    torch.save(checkpoint, MODEL_PATH)
    print(f"[bo] model checkpoint saved -> {MODEL_PATH}")


def load_model_hyperparams() -> dict | None:
    if not MODEL_PATH.exists():
        return None
    try:
        checkpoint = torch.load(MODEL_PATH, weights_only=True)
        print(f"[bo] model checkpoint loaded <- {MODEL_PATH}")
        return checkpoint
    except Exception as exc:
        print(f"[bo] WARNING: could not load model checkpoint ({exc}) - fitting from scratch", file=sys.stderr)
        return None


def next_candidate(model, X: torch.Tensor, Y: torch.Tensor) -> torch.Tensor:
    try:
        from botorch.acquisition.analytic import LogExpectedImprovement as AcqClass
        acq_name = "LogExpectedImprovement"
    except ImportError:
        from botorch.acquisition.analytic import ExpectedImprovement as AcqClass
        acq_name = "ExpectedImprovement"
    from botorch.optim import optimize_acqf

    print(f"[bo] acquisition: {acq_name}")
    acqf = AcqClass(model=model, best_f=Y.max().item())
    candidate, _ = optimize_acqf(
        acqf,
        bounds=BOUNDS,
        q=1,
        num_restarts=10,
        raw_samples=256,
    )
    return candidate.squeeze(0)


def evaluate(bo_params: dict, phase: str) -> tuple[torch.Tensor | None, torch.Tensor | None, int | None, str | None]:
    geo_params = bo_to_geo(bo_params)
    sid = next_sample_id()
    sample_dir = RESULTS_DIR / sid
    write_sample_yaml(sample_dir, bo_params)

    print(f"\n[bo] sample {sid}: phase={phase}, bo={bo_params}")
    print(f"[bo] sample {sid}: geo={{'a': {geo_params['a']:.4f}, 'd': {geo_params['d']:.4f}}}")

    result = subprocess.run(
        [
            snakemake_bin(),
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(sample_dir),
            "--cores",
            str(CORES),
            "--config",
            f"results_dir={sample_dir}",
        ]
    )
    if result.returncode != 0:
        failure_stage = "geometry_overlap" if not bool(geo_params["feasible"]) else "workflow_failed"
        failure_reason = (
            "d <= a in unconstrained TN-05 geometry"
            if not bool(geo_params["feasible"])
            else "Snakemake pipeline returned non-zero exit status"
        )
        _write_failed_objectives(sample_dir, bo_params, geo_params, phase, failure_stage, failure_reason)
        return None, None, int(sid), phase

    objective = read_objective(sample_dir)
    if objective is None:
        _write_failed_objectives(
            sample_dir,
            bo_params,
            geo_params,
            phase,
            "missing_objective",
            "pressureDrop.csv/objectives.csv did not contain a valid pressure-drop value",
        )
        return None, None, int(sid), phase

    _annotate_objectives_csv(sample_dir, bo_params, geo_params, phase)
    print(f"[bo] {sid}: J_dp={objective:.4g} Pa")

    x_t = torch.tensor([[float(bo_params[key]) for key in BO_PARAM_NAMES]])
    y_t = torch.tensor([[objective]])
    return x_t, y_t, int(sid), phase


def collect_existing():
    if not RESULTS_DIR.exists():
        return None, None, [], [], 0

    xs, ys, sample_ids, phases = [], [], [], []
    n_recorded = 0
    for sample_dir in sorted(RESULTS_DIR.iterdir()):
        if not (sample_dir.is_dir() and sample_dir.name.isdigit()):
            continue
        obj_csv = sample_dir / "objectives.csv"
        if not obj_csv.exists():
            continue

        with open(obj_csv, newline="") as handle:
            for row in csv.DictReader(handle):
                n_recorded += 1
                try:
                    feasible = int(str(row.get("feasible", "1")).strip() or "1")
                except ValueError:
                    feasible = 1
                if not feasible:
                    continue

                try:
                    if all(row.get(f"bo_{key}", "") != "" for key in BO_PARAM_NAMES):
                        bo = {key: float(row[f"bo_{key}"]) for key in BO_PARAM_NAMES}
                    else:
                        geo = {
                            "a": float(row["geo_a"]),
                            "d": float(row["geo_d"]),
                        }
                        bo = geo_to_bo(geo)
                    objective = float(row[OBJ_PDROP])
                except (KeyError, TypeError, ValueError):
                    continue

                xs.append([float(bo[key]) for key in BO_PARAM_NAMES])
                ys.append([objective])
                sample_ids.append(int(row["sample_id"]))
                phases.append(str(row.get("phase", "bo")).strip().lower() or "bo")

    if not xs:
        return None, None, [], [], n_recorded
    return torch.tensor(xs), torch.tensor(ys), sample_ids, phases, n_recorded


def main() -> None:
    from botorch.utils.sampling import draw_sobol_samples

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    X_obs, Y_obs, sample_ids_obs, phases_obs, n_recorded = collect_existing()
    if n_recorded > 0:
        print(
            f"[bo] resuming: {n_recorded} recorded sample(s) found "
            f"({0 if Y_obs is None else Y_obs.shape[0]} feasible objective observations)"
        )
        aggregate_all()
        if Y_obs is not None:
            plot_objective_history(
                sample_ids_obs,
                Y_obs[:, 0].tolist(),
                phases_obs,
                f"Pressure-drop history - resumed [{len(sample_ids_obs)} feasible sample(s)]",
            )

    n_init_done = min(n_recorded, N_INIT)
    n_init_needed = max(0, N_INIT - n_recorded)
    if n_init_needed > 0:
        print(f"[bo] === Sobol initialisation: {n_init_needed} remaining of {N_INIT} ===")
        sobol_X = draw_sobol_samples(bounds=BOUNDS, n=n_init_needed, q=1).squeeze(1)

        for x in sobol_X:
            bo_params = {key: float(x[idx]) for idx, key in enumerate(BO_PARAM_NAMES)}
            x_t, y_t, sid, phase = evaluate(bo_params, phase="sobol")
            n_recorded += 1
            if x_t is not None and y_t is not None and sid is not None and phase is not None:
                X_obs = x_t if X_obs is None else torch.cat([X_obs, x_t])
                Y_obs = y_t if Y_obs is None else torch.cat([Y_obs, y_t])
                sample_ids_obs.append(sid)
                phases_obs.append(phase)

            aggregate_all()
            if Y_obs is not None:
                plot_objective_history(
                    sample_ids_obs,
                    Y_obs[:, 0].tolist(),
                    phases_obs,
                    f"Pressure-drop history - init {min(n_recorded, N_INIT)}/{N_INIT}",
                )

    if X_obs is None or X_obs.shape[0] < 2:
        sys.exit("[bo] not enough feasible samples to fit a GP - aborting")

    n_bo_existing = max(0, n_recorded - N_INIT)
    print(
        f"\n[bo] === Sequential BO: launching {N_BO} new iteration(s) "
        f"(existing BO samples: {n_bo_existing}) ==="
    )

    warm_start = load_model_hyperparams()
    for i in range(N_BO):
        bo_iter = n_bo_existing + i + 1
        print(
            f"\n[bo] --- BO iteration {bo_iter} "
            f"(launch {i + 1}/{N_BO}; feasible observations so far: {X_obs.shape[0]}) ---"
        )

        model = fit_model(X_obs, Y_obs, warm_start=warm_start)
        save_model(model)
        warm_start = {name: value.detach().clone() for name, value in model.named_parameters()}
        x_next = next_candidate(model, X_obs, Y_obs)
        bo_params = {key: float(x_next[idx]) for idx, key in enumerate(BO_PARAM_NAMES)}

        x_t, y_t, sid, phase = evaluate(bo_params, phase="bo")
        n_recorded += 1
        if x_t is not None and y_t is not None and sid is not None and phase is not None:
            X_obs = torch.cat([X_obs, x_t])
            Y_obs = torch.cat([Y_obs, y_t])
            sample_ids_obs.append(sid)
            phases_obs.append(phase)

        aggregate_all()
        plot_objective_history(
            sample_ids_obs,
            Y_obs[:, 0].tolist(),
            phases_obs,
            f"Pressure-drop history - BO iter {bo_iter} [{len(sample_ids_obs)} feasible sample(s)]",
        )

    best_index = int(torch.argmax(Y_obs[:, 0]).item())
    best_bo = {key: float(X_obs[best_index, idx]) for idx, key in enumerate(BO_PARAM_NAMES)}
    best_geo = bo_to_geo(best_bo)
    best_value = float(Y_obs[best_index, 0].item())

    print("\n[bo] Best design so far:")
    if CAD_MODE == "constrained":
        print(
            f"  a={best_geo['a']:.4f}  d={best_geo['d']:.4f}  "
            f"d_ratio={best_bo['d_ratio']:.4f}  J_dp={best_value:.4g} Pa"
        )
    else:
        print(
            f"  a={best_geo['a']:.4f}  d={best_geo['d']:.4f}  "
            f"J_dp={best_value:.4g} Pa"
        )

    aggregate_all()
    print("\n[bo] === Generating visualizations ===")
    vis_ret = subprocess.run(
        [
            sys.executable,
            str(CASE_ROOT / "visualize_results.py"),
            "--results-dir",
            str(RESULTS_DIR),
            "--output-dir",
            str(RESULTS_DIR / "visualizations"),
        ]
    )
    if vis_ret.returncode != 0:
        print("[bo] WARNING: visualization failed (non-fatal)", file=sys.stderr)


if __name__ == "__main__":
    main()
