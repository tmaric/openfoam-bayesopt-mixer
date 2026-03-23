#!/usr/bin/env python3
"""Bayesian optimisation loop for the SAR lamination ladder mixer.

Design variables (normalised, H = 1):
    w_s   – subchannel half-gap after split     [0.20, 0.44]
    t_s   – splitter thickness, split section   [0.04, 0.15]
    t_m   – splitter thickness, merge section   [0.02, 0.10]
    delta – top deflector vertical bias         [0.00, 0.15]

Objectives:
    J_dp  – pressure drop (Pa)          minimise   column: pdrop_pressure_drop_Pa
    J_mix – intensity of segregation    minimise   column: mixing_intensity_of_segregation
    (displayed as mixing quality = 1 - J_mix, higher = better mixed)

The loop supports resumption: completed samples in --results-dir are loaded
first, then BO picks up from where it left off.

After every new sample the Pareto front is redrawn and saved as
    results/pareto_front.png

Usage (from SplitAndRecombineMixer/):
    python bayes_optimize.py [--n-init N] [--n-bo N] [--cores N]
                             [--results-dir PATH] [--no-plot]
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import torch
import yaml

torch.set_default_dtype(torch.float64)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CASE_ROOT = Path(__file__).resolve().parent
SNAKEFILE = CASE_ROOT / "Snakefile"
TEMPLATE_YAML = CASE_ROOT / "SplitAndRecombineHydro" / "sar_mixer_cad.yaml"

# ---------------------------------------------------------------------------
# Design space
# ---------------------------------------------------------------------------
PARAM_NAMES = ["w_s", "t_s", "t_m", "delta"]

BOUNDS = torch.tensor(
    [
        [0.20, 0.04, 0.02, 0.00],  # lower bounds
        [0.44, 0.15, 0.10, 0.15],  # upper bounds
    ]
)

# ---------------------------------------------------------------------------
# Objective column names in objectives.csv
# ---------------------------------------------------------------------------
OBJ_PDROP = "pdrop_pressure_drop_Pa"
OBJ_MIX = "mixing_intensity_of_segregation"


# ---------------------------------------------------------------------------
# Pareto front plot
# ---------------------------------------------------------------------------

def plot_pareto_front(
    results_root: Path,
    X_obs: torch.Tensor,
    Y_obs: torch.Tensor,
    n_init: int,
    label: str = "",
) -> None:
    """
    Save results/pareto_front.png.

    x-axis: pressure drop J_dp (Pa)          — lower left is better
    y-axis: mixing quality = 1 - I_s (–)    — upper left is better

    Points are coloured by phase:
        steel-blue  — Sobol initialisation samples
        orange      — BO-suggested samples
    Pareto-optimal points are additionally ringed in black.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from botorch.utils.multi_objective.pareto import is_pareto_efficient
    except ImportError:
        print("[bo] matplotlib or botorch not available — skipping plot")
        return

    j_dp = Y_obs[:, 0].numpy()
    j_mix = Y_obs[:, 1].numpy()
    mixing_quality = 1.0 - j_mix  # higher = better mixed

    n = len(j_dp)
    point_colors = [
        "steelblue" if i < n_init else "darkorange" for i in range(n)
    ]

    # is_pareto_efficient expects maximisation; negate both minimised objectives
    pareto_mask = is_pareto_efficient(-Y_obs).numpy()

    fig, ax = plt.subplots(figsize=(7, 5))

    # ---- scatter all points ------------------------------------------------
    non_pareto = ~pareto_mask
    if non_pareto.any():
        ax.scatter(
            j_dp[non_pareto], mixing_quality[non_pareto],
            c=[point_colors[i] for i in range(n) if non_pareto[i]],
            alpha=0.5, s=40,
        )
    if pareto_mask.any():
        ax.scatter(
            j_dp[pareto_mask], mixing_quality[pareto_mask],
            c=[point_colors[i] for i in range(n) if pareto_mask[i]],
            s=90, edgecolors="black", linewidths=1.4, zorder=5,
        )

    # ---- Pareto staircase --------------------------------------------------
    if pareto_mask.sum() > 1:
        px = j_dp[pareto_mask]
        py = mixing_quality[pareto_mask]
        order = px.argsort()
        px, py = px[order], py[order]
        # extend staircase to the right (worst pressure drop) at bottom
        ax.step(px, py, where="post", color="black", lw=1.5, alpha=0.8)
        ax.scatter(px, py, color="none", edgecolors="black", s=90, zorder=6)

    # ---- legend ------------------------------------------------------------
    n_bo_done = max(0, n - n_init)
    handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="steelblue", markersize=9,
               label=f"Sobol init  ({min(n_init, n)})"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="darkorange", markersize=9,
               label=f"BO suggested  ({n_bo_done})"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="gray", markeredgecolor="black",
               markeredgewidth=1.4, markersize=9,
               label=f"Pareto-optimal  ({int(pareto_mask.sum())})"),
        Line2D([0], [0], color="black", lw=1.5, label="Pareto front"),
    ]
    ax.legend(handles=handles, fontsize=9, loc="lower left")

    ax.set_xlabel(r"Pressure drop  $J_\mathrm{dp}$  (Pa)", fontsize=11)
    ax.set_ylabel(r"Mixing quality  $1 - I_s$  (–)", fontsize=11)
    title = f"SAR mixer — Pareto front  [{n} sample(s)"
    if label:
        title += f",  {label}"
    title += "]"
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)

    out = results_root / "pareto_front.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[bo] plot saved: {out}")


# ---------------------------------------------------------------------------
# Sample management
# ---------------------------------------------------------------------------

def next_sample_id(results_root: Path) -> str:
    """Return the next zero-padded 5-digit sample ID."""
    existing = (
        sorted(int(d.name) for d in results_root.iterdir() if d.is_dir() and d.name.isdigit())
        if results_root.exists()
        else []
    )
    nxt = existing[-1] + 1 if existing else 0
    return f"{nxt:05d}"


def write_sample_yaml(sample_dir: Path, params: dict) -> None:
    """Write sar_mixer_cad.yaml into sample_dir, merging params into the template."""
    with open(TEMPLATE_YAML) as fh:
        base = yaml.safe_load(fh)
    base.update(params)
    sample_dir.mkdir(parents=True, exist_ok=True)
    with open(sample_dir / "sar_mixer_cad.yaml", "w") as fh:
        yaml.dump(base, fh, default_flow_style=False, sort_keys=False)


def run_snakemake(sample_dir: Path, cores: int) -> bool:
    """
    Run Snakemake for one sample.

    --snakefile  keeps the DAG definition in the repo
    --directory  isolates .snakemake/ metadata so parallel runs don't collide
    --config     tells the Snakefile where to stage results
    """
    cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--directory", str(sample_dir),
        "--cores", str(cores),
        "--config", f"results_dir={sample_dir}",
    ]
    print(f"[bo] snakemake {sample_dir.name}: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    return result.returncode == 0


def read_objectives(sample_dir: Path):
    """Return (J_dp, J_mix) floats from objectives.csv, or (None, None)."""
    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        return None, None
    with open(obj_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                return float(row[OBJ_PDROP]), float(row[OBJ_MIX])
            except (KeyError, ValueError):
                return None, None
    return None, None


def collect_existing(results_root: Path):
    """
    Load all completed samples from results_root.
    Returns (X, Y) tensors of shape (n, 4) and (n, 2), or (None, None).
    """
    xs, ys = [], []
    if not results_root.exists():
        return None, None
    for d in sorted(results_root.iterdir()):
        if not (d.is_dir() and d.name.isdigit()):
            continue
        params_yaml = d / "sar_mixer_cad.yaml"
        if not params_yaml.exists():
            continue
        j_dp, j_mix = read_objectives(d)
        if j_dp is None:
            continue
        with open(params_yaml) as fh:
            p = yaml.safe_load(fh)
        xs.append([float(p[k]) for k in PARAM_NAMES])
        ys.append([j_dp, j_mix])
    if not xs:
        return None, None
    return torch.tensor(xs), torch.tensor(ys)


def run_sample(results_root: Path, params: dict, cores: int):
    """
    Assign ID, write YAML, run Snakemake.
    Returns (x_tensor (1,4), y_tensor (1,2)) or (None, None).
    """
    sid = next_sample_id(results_root)
    sample_dir = results_root / sid
    print(f"\n[bo] sample {sid}: {params}")
    write_sample_yaml(sample_dir, params)

    if not run_snakemake(sample_dir, cores):
        print(f"[bo] WARNING: sample {sid} failed", file=sys.stderr)
        return None, None

    j_dp, j_mix = read_objectives(sample_dir)
    if j_dp is None:
        print(f"[bo] WARNING: no objectives for {sid}", file=sys.stderr)
        return None, None

    print(f"[bo] sample {sid}: J_dp={j_dp:.4g} Pa,  J_mix(I_s)={j_mix:.4g}"
          f"  →  mixing quality={1.0 - j_mix:.4g}")
    x = torch.tensor([[float(params[k]) for k in PARAM_NAMES]])
    y = torch.tensor([[j_dp, j_mix]])
    return x, y


# ---------------------------------------------------------------------------
# GP model and acquisition
# ---------------------------------------------------------------------------

def fit_model(X: torch.Tensor, Y: torch.Tensor):
    """Fit a batched SingleTaskGP (one output per objective, negated for maximisation)."""
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.mlls import ExactMarginalLogLikelihood

    Y_neg = -Y  # BOTorch maximises; we minimise both objectives
    model = SingleTaskGP(
        X,
        Y_neg,
        input_transform=Normalize(d=X.shape[-1]),
        outcome_transform=Standardize(m=Y_neg.shape[-1]),
    )
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    return model


def suggest_next(model, X_obs: torch.Tensor, Y_obs: torch.Tensor) -> torch.Tensor:
    """Return the next candidate point via qNEHVI."""
    from botorch.acquisition.multi_objective import (
        qNoisyExpectedHypervolumeImprovement,
    )
    from botorch.optim import optimize_acqf

    Y_neg = -Y_obs
    # Reference point: slightly below the worst observed value in each objective
    Y_range = (Y_neg.max(dim=0).values - Y_neg.min(dim=0).values).clamp(min=1e-6)
    ref_point = Y_neg.min(dim=0).values - 0.1 * Y_range

    acqf = qNoisyExpectedHypervolumeImprovement(
        model=model,
        ref_point=ref_point,
        X_baseline=X_obs,
        prune_baseline=True,
    )
    candidate, _ = optimize_acqf(
        acqf,
        bounds=BOUNDS,
        q=1,
        num_restarts=10,
        raw_samples=256,
    )
    return candidate.squeeze(0)


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def aggregate_all(results_root: Path) -> None:
    """Write all_objectives.csv aggregating every per-sample objectives.csv."""
    rows = []
    fieldnames = None
    for d in sorted(results_root.iterdir()):
        if not (d.is_dir() and d.name.isdigit()):
            continue
        obj_csv = d / "objectives.csv"
        if not obj_csv.exists():
            continue
        with open(obj_csv, newline="") as fh:
            reader = csv.DictReader(fh)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            for row in reader:
                rows.append(row)
    if not rows:
        return
    out = results_root / "all_objectives.csv"
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[bo] wrote {out} ({len(rows)} samples)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="BO loop for the SAR lamination mixer")
    ap.add_argument(
        "--n-init", type=int, default=8,
        help="Number of initial Sobol samples (default: 8)",
    )
    ap.add_argument(
        "--n-bo", type=int, default=20,
        help="Number of BO iterations after initialisation (default: 20)",
    )
    ap.add_argument(
        "--cores", type=int, default=4,
        help="CPU cores per Snakemake run (default: 4)",
    )
    ap.add_argument(
        "--results-dir", type=Path, default=CASE_ROOT / "results",
        help="Root directory for per-sample result subdirectories",
    )
    ap.add_argument(
        "--no-plot", action="store_true",
        help="Disable Pareto front PNG output",
    )
    args = ap.parse_args()

    results_root: Path = args.results_dir.resolve()
    results_root.mkdir(parents=True, exist_ok=True)
    do_plot = not args.no_plot

    # ---- resume from completed samples -------------------------------------
    X_obs, Y_obs = collect_existing(results_root)
    n_existing = 0 if X_obs is None else X_obs.shape[0]
    print(f"[bo] found {n_existing} completed sample(s) in {results_root}")

    # n_init tracks the boundary between Sobol and BO points in X_obs/Y_obs
    n_init_effective = min(n_existing, args.n_init)

    if n_existing > 0 and do_plot:
        plot_pareto_front(results_root, X_obs, Y_obs, n_init_effective, label="resumed")

    # ---- initial Sobol phase -----------------------------------------------
    n_init_needed = max(0, args.n_init - n_existing)
    if n_init_needed > 0:
        from botorch.utils.sampling import draw_sobol_samples

        print(f"[bo] drawing {n_init_needed} initial Sobol sample(s)")
        sobol_X = draw_sobol_samples(bounds=BOUNDS, n=n_init_needed, q=1).squeeze(1)

        for x in sobol_X:
            params = {k: float(x[j]) for j, k in enumerate(PARAM_NAMES)}
            x_t, y_t = run_sample(results_root, params, args.cores)
            if x_t is None:
                continue
            X_obs = x_t if X_obs is None else torch.cat([X_obs, x_t])
            Y_obs = y_t if Y_obs is None else torch.cat([Y_obs, y_t])
            n_init_effective += 1
            if do_plot:
                plot_pareto_front(
                    results_root, X_obs, Y_obs, n_init_effective,
                    label=f"init {n_init_effective}/{args.n_init}",
                )

    if X_obs is None or X_obs.shape[0] < 2:
        print("[bo] not enough successful samples to fit GP — exiting", file=sys.stderr)
        sys.exit(1)

    # ---- BO loop -----------------------------------------------------------
    for iteration in range(args.n_bo):
        n_total = X_obs.shape[0]
        print(f"\n[bo] === BO iteration {iteration + 1}/{args.n_bo} "
              f"(total observed: {n_total}) ===")

        model = fit_model(X_obs, Y_obs)
        x_next = suggest_next(model, X_obs, Y_obs)
        params = {k: float(x_next[j]) for j, k in enumerate(PARAM_NAMES)}

        x_t, y_t = run_sample(results_root, params, args.cores)
        if x_t is None:
            continue
        X_obs = torch.cat([X_obs, x_t])
        Y_obs = torch.cat([Y_obs, y_t])

        if do_plot:
            plot_pareto_front(
                results_root, X_obs, Y_obs, n_init_effective,
                label=f"BO iter {iteration + 1}/{args.n_bo}",
            )

    # ---- aggregate ---------------------------------------------------------
    aggregate_all(results_root)

    # ---- report Pareto front -----------------------------------------------
    from botorch.utils.multi_objective.pareto import is_pareto_efficient

    pareto_mask = is_pareto_efficient(-Y_obs)
    pareto_X = X_obs[pareto_mask]
    pareto_Y = Y_obs[pareto_mask]

    print("\n[bo] Pareto-optimal designs:")
    print(f"  {'w_s':>6}  {'t_s':>6}  {'t_m':>6}  {'delta':>6}"
          f"  {'J_dp [Pa]':>12}  {'mix quality':>12}")
    for x, y in zip(pareto_X.tolist(), pareto_Y.tolist()):
        print(f"  {x[0]:6.3f}  {x[1]:6.3f}  {x[2]:6.3f}  {x[3]:6.3f}"
              f"  {y[0]:12.4g}  {1.0 - y[1]:12.4f}")

    print("\n[bo] optimisation complete.")


if __name__ == "__main__":
    main()
