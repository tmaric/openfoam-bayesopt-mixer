#!/usr/bin/env python3
"""Simple sequential multi-objective Bayesian optimisation for the SAR mixer.

Objectives (both minimised):
    J_dp  – pressure drop (Pa)          pdrop_pressure_drop_Pa
    J_mix – intensity of segregation    mixing_intensity_of_segregation
    (plotted as mixing quality = 1 - J_mix, so higher = better mixed)

Design variables (Technical Note 02, normalised by H = 1 / L_cell = 4):
    w_s   [0.25, 0.45]   subchannel half-gap after split       (fraction of H)
    t_s   [0.02, 0.15]   splitter thickness – split section    (fraction of H)
    t_m   [0.02, 0.15]   splitter thickness – merge section    (fraction of H)
    L_s   [0.80, 1.80]   split section length  (0.20–0.45 * L_cell, L_cell=4)
    L_m   [0.80, 1.80]   merge section length  (0.20–0.45 * L_cell, L_cell=4)
    delta [0.00, 0.15]   top deflector vertical bias           (fraction of H)

Feasibility constraints (linear, from CAD geometry; see sar_mixer_cad.py G1-G6):
    C1 (≈G5): w_s - 0.5*t_s - delta >= 0.02
        Fluid gap between deflector peak and splitter surface must be resolvable.
        Includes the top-deflector delta shift.  Stricter than the CAD G5 alone.
    C2  (G1): L_s + L_m <= 3.60
        Interaction length L_c = L_cell - L_s - L_m >= 0.4 (10% of L_cell).
        CAD G1 allows L_c >= 0.01; this BO constraint is intentionally tighter.
    C3  (G4): t_s - t_m >= 0.011
        Split splitter must be strictly wider than merge splitter so that the
        Boolean cut at each cell boundary acts on real material (not a void).
        Threshold is 10% above CAD's _MESH_MIN (0.01) to prevent FP boundary
        failures when physical values t_s*H, t_m*H are compared.
    G2, G3a/b, G6 are guaranteed by the parameter bounds (see is_feasible()).

Algorithm:
    1. Sobol quasi-random initialisation  (N_INIT samples, feasibility-filtered)
    2. Sequential BO loop  (N_BO iterations):
         fit SingleTaskGP → optimise qNEHVI with C1,C2 → evaluate → update

Usage:
    cd SplitAndRecombineMixer
    python bayes_optimize_sequential.py
"""

import csv
import subprocess
import sys
from pathlib import Path

import torch
import yaml

torch.set_default_dtype(torch.float64)

# ---------------------------------------------------------------------------
# Configuration — edit bayes_optimize_sequential.yaml to change the run
# ---------------------------------------------------------------------------
CASE_ROOT   = Path(__file__).resolve().parent

_cfg_path = CASE_ROOT / "bayes_optimize_sequential.yaml"
with open(_cfg_path) as _f:
    _cfg = yaml.safe_load(_f)

N_INIT  = int(_cfg["n_init"])   # Sobol initialisation samples
N_BO    = int(_cfg["n_bo"])     # BO iterations
CORES   = int(_cfg["cores"])    # CPU cores per Snakemake call
RESULTS_DIR = CASE_ROOT / "results"
SNAKEFILE   = CASE_ROOT / "Snakefile"
TEMPLATE_YAML = CASE_ROOT / "SplitAndRecombineHydro" / "sar_mixer_cad.yaml"
MODEL_PATH    = CASE_ROOT / "SplitAndRecombineMixer.pt"

PARAM_NAMES = ["w_s", "t_s", "t_m", "L_s", "L_m", "delta"]
#               idx:    0      1      2      3      4      5

BOUNDS = torch.tensor([
    [0.25, 0.02, 0.02, 0.80, 0.80, 0.00],   # lower
    [0.45, 0.15, 0.15, 1.80, 1.80, 0.15],   # upper
])

L_CELL = 4.0   # fixed unit-cell length (normalised units, same as YAML)

# Feasibility constraints  (checked before every Snakemake call)
# C1: w_s - 0.5*t_s - delta >= 0.02   → deflector–splitter clearance
# C2: L_s + L_m             <= 3.60   → L_c = L_cell - L_s - L_m >= 0.4
CLEARANCE_MIN = 0.02
LC_MIN        = 0.4   # = L_CELL - 3.60

# BOTorch inequality_constraints format: (indices, coefficients, rhs)
# where  sum(coeff_i * x[idx_i]) >= rhs
TM_MARGIN = 0.011  # t_s must exceed t_m by at least this amount (normalised)
                   # Slightly above CAD's _MESH_MIN=0.01 to avoid FP boundary
                   # failures when t_s*H - t_m*H is evaluated in physical units.

INEQ_CONSTRAINTS = [
    # C1: 1*w_s  - 0.5*t_s  + 0*t_m  + 0*L_s  + 0*L_m  - 1*delta >= 0.02
    (torch.tensor([0, 1, 5]), torch.tensor([1.0, -0.5, -1.0]), CLEARANCE_MIN),
    # C2: 0*w_s  + 0*t_s   + 0*t_m  - 1*L_s  - 1*L_m  + 0*delta >= -(L_CELL - LC_MIN)
    (torch.tensor([3, 4]),    torch.tensor([-1.0, -1.0]),       -(L_CELL - LC_MIN)),
    # C3: 0*w_s  + 1*t_s   - 1*t_m  + 0*L_s  + 0*L_m  + 0*delta >= TM_MARGIN
    (torch.tensor([1, 2]),    torch.tensor([1.0, -1.0]),         TM_MARGIN),
]

OBJ_PDROP = "pdrop_pressure_drop_Pa"
OBJ_MIX   = "mixing_intensity_of_segregation"

# Penalty values returned (and written to objectives.csv) when a Snakemake run
# fails at any step.  The GP observes these as extremely bad outcomes and will
# steer future candidates away from the failed region.
#   PENALTY_PDROP — far above any physical pressure drop in the mixer (Pa)
#   PENALTY_MIX   — 1.0 is the worst possible intensity of segregation (fully unmixed)
PENALTY_PDROP = 1e3
PENALTY_MIX   = 1.0


def is_feasible(params: dict) -> bool:
    """Return True when all CAD geometry constraints are satisfied.

    Mirrors sar_mixer_cad.py::_check_geometry() in normalised units (H = 1).
    All G-constraints from _check_geometry() are covered:

      C1  ↔  G5 + top-deflector gap:
              w_s - 0.5*t_s - delta >= CLEARANCE_MIN (0.02)
              Stricter than G5 alone (G5 threshold 0.01) because delta shifts
              the top deflector toward the splitter.
      C2  ↔  G1:
              L_c = L_CELL - L_s - L_m >= LC_MIN (0.4)
              Intentionally conservative — CAD allows L_c >= 0.01.
      C3  ↔  G4:
              t_s - t_m >= TM_MARGIN (0.011)  (10% above _MESH_MIN = 0.01·H
              to guard against FP rounding when physical values t_s*H, t_m*H
              are subtracted and compared to _MESH_MIN)

    G2, G3a, G3b, G6 are guaranteed by the parameter bounds and are NOT active
    constraints here:
      G2  h_d = 0.5 - w_s >= 0.05  (w_s <= 0.45)              >> 0.01
      G3a t_s >= 0.02  (lower bound)                            > 0.01
      G3b t_m >= 0.02  (lower bound)                            > 0.01
      G6  2*w_s - delta >= 0.35  (w_s >= 0.25, delta <= 0.15)  >> 0.01
    """
    w_s   = params["w_s"]
    t_s   = params["t_s"]
    t_m   = params["t_m"]
    L_s   = params["L_s"]
    L_m   = params["L_m"]
    delta = params["delta"]
    c1 = w_s - 0.5 * t_s - delta >= CLEARANCE_MIN   # G5 + top-deflector clearance
    c2 = L_s + L_m <= L_CELL - LC_MIN               # G1: interaction length
    c3 = t_s - t_m >= TM_MARGIN                     # G4: split splitter wider than merge
    return c1 and c2 and c3


# ---------------------------------------------------------------------------
# Pareto front plot  (saved after every new sample)
# ---------------------------------------------------------------------------

def plot_pareto_front(X: torch.Tensor, Y: torch.Tensor, n_init: int, title: str) -> None:
    """Save RESULTS_DIR/pareto_front.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from botorch.utils.multi_objective.pareto import is_non_dominated

    j_dp = Y[:, 0].numpy()
    mq   = 1.0 - Y[:, 1].numpy()   # mixing quality: higher = better

    n = len(j_dp)
    colors = ["steelblue" if i < n_init else "darkorange" for i in range(n)]
    pareto = is_non_dominated(-Y).numpy()   # negate: minimisation → maximisation

    fig, ax = plt.subplots(figsize=(7, 5))

    # non-Pareto scatter
    if (~pareto).any():
        ax.scatter(j_dp[~pareto], mq[~pareto],
                   c=[colors[i] for i in range(n) if not pareto[i]],
                   alpha=0.5, s=40)

    # Pareto scatter (ringed)
    if pareto.any():
        ax.scatter(j_dp[pareto], mq[pareto],
                   c=[colors[i] for i in range(n) if pareto[i]],
                   s=90, edgecolors="black", linewidths=1.4, zorder=5)

    # Pareto staircase
    if pareto.sum() > 1:
        px, py = j_dp[pareto], mq[pareto]
        order  = px.argsort()
        ax.step(px[order], py[order], where="post",
                color="black", lw=1.5, alpha=0.8)

    ax.set_xlabel(r"Pressure drop  $J_\mathrm{dp}$  (Pa, log scale)", fontsize=11)
    ax.set_ylabel(r"Mixing quality  $1 - I_s$  (–)", fontsize=11)
    ax.set_xscale("log")
    ax.set_title(title, fontsize=11)
    ax.grid(True, which="both", alpha=0.3)

    handles = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor="steelblue",
               markersize=9, label=f"Sobol init ({min(n_init, n)})"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="darkorange",
               markersize=9, label=f"BO suggested ({max(0, n - n_init)})"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor="black", markeredgewidth=1.4, markersize=9,
               label=f"Pareto ({int(pareto.sum())})"),
        Line2D([0],[0], color="black", lw=1.5, label="Pareto front"),
    ]
    ax.legend(handles=handles, fontsize=9, loc="lower left")

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "pareto_front.png", dpi=150)
    plt.close(fig)
    print(f"[bo] plot saved → {RESULTS_DIR / 'pareto_front.png'}")


# ---------------------------------------------------------------------------
# Sample execution
# ---------------------------------------------------------------------------

def aggregate_all() -> None:
    """Rewrite results/all_samples.csv from all completed per-sample objectives.csv files.

    Each per-sample objectives.csv is a single flat row:
        sample_id, results_dir, geo_<param>, ..., pdrop_<col>, ..., mixing_<col>, ...

    Failed samples have fewer columns (penalty rows).  The union of all fieldnames
    is used; missing fields are written as empty strings so pd.read_csv / pd.concat
    can handle them gracefully (they become NaN in pandas).
    """
    rows, all_fieldnames = [], []
    seen_fields: set = set()
    for d in sorted(RESULTS_DIR.iterdir()):
        if not (d.is_dir() and d.name.isdigit()):
            continue
        obj_csv = d / "objectives.csv"
        if not obj_csv.exists():
            continue
        with open(obj_csv, newline="") as fh:
            reader = csv.DictReader(fh)
            for field in (reader.fieldnames or []):
                if field not in seen_fields:
                    all_fieldnames.append(field)
                    seen_fields.add(field)
            rows.extend(reader)
    if not rows:
        return
    out = RESULTS_DIR / "all_samples.csv"
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=all_fieldnames, extrasaction="ignore",
                                restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[bo] all_samples.csv → {len(rows)} sample(s)  ({out})")


def next_sample_id() -> str:
    ids = [int(d.name) for d in RESULTS_DIR.iterdir()
           if d.is_dir() and d.name.isdigit()] if RESULTS_DIR.exists() else []
    return f"{(max(ids) + 1 if ids else 0):05d}"


def _write_penalty_objectives(sample_dir: Path, params: dict) -> None:
    """Write a minimal objectives.csv with penalty values for a failed sample.

    This ensures that collect_existing() can reconstruct the penalty observation
    on resume, so the GP consistently avoids the failed region across restarts.
    The 'failed' column flags the row for post-hoc analysis.
    """
    row = {
        "sample_id":   sample_dir.name,
        "results_dir": str(sample_dir.resolve()),
        "failed":      "True",
        **{f"geo_{k}": params[k] for k in PARAM_NAMES},
        OBJ_PDROP: PENALTY_PDROP,
        OBJ_MIX:   PENALTY_MIX,
    }
    out = sample_dir / "objectives.csv"
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    print(f"[bo] penalty objectives.csv written → {out}")


def evaluate(params: dict) -> tuple[float, float] | tuple[None, None]:
    """Write YAML, run Snakemake, return (J_dp, J_mix) or (None, None).

    Returns (None, None) only when the parameter vector is infeasible (pre-filter).
    For Snakemake failures at any pipeline step, returns (PENALTY_PDROP, PENALTY_MIX)
    and persists those values to objectives.csv so the GP observes the failure and
    steers future candidates away from the failed region.
    """
    if not is_feasible(params):
        print(f"[bo] SKIP infeasible params: {params}", file=sys.stderr)
        return None, None

    sid        = next_sample_id()
    sample_dir = RESULTS_DIR / sid
    sample_dir.mkdir(parents=True, exist_ok=True)

    # write geometry YAML for this sample
    with open(TEMPLATE_YAML) as fh:
        geo = yaml.safe_load(fh)
    geo.update(params)
    with open(sample_dir / "sar_mixer_cad.yaml", "w") as fh:
        yaml.dump(geo, fh, default_flow_style=False, sort_keys=False)

    print(f"\n[bo] sample {sid}: {params}")

    ret = subprocess.run([
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--directory", str(sample_dir),
        "--cores",     str(CORES),
        "--config",    f"results_dir={sample_dir}",
    ])
    if ret.returncode != 0:
        print(f"[bo] WARNING: sample {sid} failed — penalising", file=sys.stderr)
        _write_penalty_objectives(sample_dir, params)
        return PENALTY_PDROP, PENALTY_MIX

    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        print(f"[bo] WARNING: objectives.csv missing for {sid} — penalising",
              file=sys.stderr)
        _write_penalty_objectives(sample_dir, params)
        return PENALTY_PDROP, PENALTY_MIX

    with open(obj_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            j_dp  = float(row[OBJ_PDROP])
            j_mix = float(row[OBJ_MIX])
            print(f"[bo] {sid}: J_dp={j_dp:.4g} Pa  "
                  f"mixing quality={1.0 - j_mix:.4f}")
            return j_dp, j_mix

    return None, None


# ---------------------------------------------------------------------------
# BO helpers
# ---------------------------------------------------------------------------

def fit_model(X: torch.Tensor, Y: torch.Tensor,
              warm_start: dict | None = None):
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.mlls import ExactMarginalLogLikelihood

    model = SingleTaskGP(X, -Y,                          # negate: BOTorch maximises
                         input_transform=Normalize(d=X.shape[-1]),
                         outcome_transform=Standardize(m=Y.shape[-1]))
    if warm_start is not None:
        try:
            with torch.no_grad():
                for name, param in model.named_parameters():
                    if name in warm_start and warm_start[name].shape == param.shape:
                        param.data.copy_(warm_start[name])
        except Exception as exc:
            print(f"[bo] WARNING: warm start failed ({exc}) — using default init",
                  file=sys.stderr)
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))
    return model


def save_model(model) -> None:
    """Persist GP hyperparameters to MODEL_PATH for warm-starting future runs."""
    checkpoint = {k: v.detach().clone() for k, v in model.named_parameters()}
    torch.save(checkpoint, MODEL_PATH)
    print(f"[bo] model checkpoint saved → {MODEL_PATH}")


def load_model_hyperparams() -> dict | None:
    """Load GP hyperparameters from MODEL_PATH, or return None if absent/corrupt."""
    if not MODEL_PATH.exists():
        return None
    try:
        ckpt = torch.load(MODEL_PATH, weights_only=True)
        print(f"[bo] model checkpoint loaded ← {MODEL_PATH}")
        return ckpt
    except Exception as exc:
        print(f"[bo] WARNING: could not load model checkpoint ({exc}) — fitting from scratch",
              file=sys.stderr)
        return None


def next_candidate(model, X: torch.Tensor, Y: torch.Tensor) -> torch.Tensor:
    from botorch.acquisition.multi_objective import qNoisyExpectedHypervolumeImprovement
    from botorch.optim import optimize_acqf

    Y_neg     = -Y
    Y_range   = (Y_neg.max(0).values - Y_neg.min(0).values).clamp(min=1e-6)
    ref_point = Y_neg.min(0).values - 0.1 * Y_range

    acqf = qNoisyExpectedHypervolumeImprovement(
        model=model, ref_point=ref_point, X_baseline=X, prune_baseline=True)

    cand, _ = optimize_acqf(acqf, bounds=BOUNDS, q=1,
                             num_restarts=10, raw_samples=256,
                             inequality_constraints=INEQ_CONSTRAINTS)
    return cand.squeeze(0)


# ---------------------------------------------------------------------------
# Resume helpers
# ---------------------------------------------------------------------------

def collect_existing() -> tuple:
    """Load all completed samples from RESULTS_DIR.

    Reads each sample's objectives.csv (written by postprocessing_agglomeration.py).
    The BO parameters are stored under geo_<name> columns; objectives under
    OBJ_PDROP and OBJ_MIX column names.

    Returns (X_obs, Y_obs, n) where X_obs is (n, 6) and Y_obs is (n, 2),
    or (None, None, 0) when no completed samples exist.
    """
    if not RESULTS_DIR.exists():
        return None, None, 0
    xs, ys = [], []
    for d in sorted(RESULTS_DIR.iterdir()):
        if not (d.is_dir() and d.name.isdigit()):
            continue
        obj_csv = d / "objectives.csv"
        if not obj_csv.exists():
            continue
        with open(obj_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    x     = [float(row[f"geo_{k}"]) for k in PARAM_NAMES]
                    j_dp  = float(row[OBJ_PDROP])
                    j_mix = float(row[OBJ_MIX])
                except (KeyError, ValueError):
                    continue
                xs.append(x)
                ys.append([j_dp, j_mix])
    if not xs:
        return None, None, 0
    return torch.tensor(xs), torch.tensor(ys), len(xs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ---- resume from completed samples ------------------------------------
    X_obs, Y_obs, n_existing = collect_existing()
    n_init_done = min(n_existing, N_INIT)
    if n_existing > 0:
        n_bo_existing = max(0, n_existing - N_INIT)
        print(f"[bo] resuming: {n_existing} completed sample(s) found "
              f"({n_init_done} init, {n_bo_existing} BO)")
        aggregate_all()
        plot_pareto_front(X_obs, Y_obs, n_init_done,
                          f"SAR mixer — resumed  [{n_existing} sample(s)]")

    # ---- Sobol initialisation (only remaining samples) --------------------
    n_init_needed = max(0, N_INIT - n_existing)
    if n_init_needed > 0:
        from botorch.utils.sampling import draw_sobol_samples

        print(f"[bo] === Sobol initialisation: {n_init_needed} remaining of {N_INIT} ===")
        # Over-sample by 5x to ensure enough feasible points after constraint filtering
        sobol_X = draw_sobol_samples(bounds=BOUNDS, n=N_INIT * 5, q=1).squeeze(1)
        sobol_X = sobol_X[
            torch.tensor([
                is_feasible({k: float(sobol_X[i, j]) for j, k in enumerate(PARAM_NAMES)})
                for i in range(sobol_X.shape[0])
            ])
        ]
        if sobol_X.shape[0] < n_init_needed:
            print(f"[bo] WARNING: only {sobol_X.shape[0]} feasible Sobol points "
                  f"(wanted {n_init_needed}) — proceeding with available", file=sys.stderr)
        sobol_X = sobol_X[:n_init_needed]

        for x in sobol_X:
            params = {k: float(x[j]) for j, k in enumerate(PARAM_NAMES)}
            j_dp, j_mix = evaluate(params)
            if j_dp is None:
                continue

            x_t = torch.tensor([[float(params[k]) for k in PARAM_NAMES]])
            y_t = torch.tensor([[j_dp, j_mix]])
            X_obs = x_t if X_obs is None else torch.cat([X_obs, x_t])
            Y_obs = y_t if Y_obs is None else torch.cat([Y_obs, y_t])
            n_init_done += 1

            aggregate_all()
            plot_pareto_front(X_obs, Y_obs, n_init_done,
                              f"SAR mixer — init {n_init_done}/{N_INIT}"
                              f"  [{X_obs.shape[0]} sample(s)]")

    if X_obs is None or X_obs.shape[0] < 2:
        sys.exit("[bo] not enough successful samples to fit a GP — aborting")

    # ---- BO loop (only remaining iterations) ------------------------------
    n_bo_done      = max(0, n_existing - N_INIT)
    n_bo_remaining = max(0, N_BO - n_bo_done)
    print(f"\n[bo] === Sequential BO: {n_bo_remaining} remaining iteration(s) "
          f"(of {N_BO}; {n_bo_done} already done) ===")

    warm_start = load_model_hyperparams()

    for i in range(n_bo_remaining):
        bo_iter = n_bo_done + i + 1
        print(f"\n[bo] --- BO iteration {bo_iter}/{N_BO}  "
              f"(observed so far: {X_obs.shape[0]}) ---")

        model = fit_model(X_obs, Y_obs, warm_start=warm_start)
        save_model(model)
        warm_start = {k: v.detach().clone() for k, v in model.named_parameters()}
        x_next = next_candidate(model, X_obs, Y_obs)
        params = {k: float(x_next[j]) for j, k in enumerate(PARAM_NAMES)}

        j_dp, j_mix = evaluate(params)
        if j_dp is None:
            continue

        x_t = torch.tensor([[float(params[k]) for k in PARAM_NAMES]])
        y_t = torch.tensor([[j_dp, j_mix]])
        X_obs = torch.cat([X_obs, x_t])
        Y_obs = torch.cat([Y_obs, y_t])

        aggregate_all()
        plot_pareto_front(X_obs, Y_obs, n_init_done,
                          f"SAR mixer — BO iter {bo_iter}/{N_BO}"
                          f"  [{X_obs.shape[0]} sample(s)]")

    # ---- final report -----------------------------------------------------
    from botorch.utils.multi_objective.pareto import is_non_dominated

    pareto = is_non_dominated(-Y_obs)
    print("\n[bo] Pareto-optimal designs:")
    print(f"  {'w_s':>6}  {'t_s':>6}  {'t_m':>6}  {'L_s':>6}  {'L_m':>6}  {'delta':>6}"
          f"  {'J_dp [Pa]':>12}  {'mix quality':>12}")
    for x, y in zip(X_obs[pareto].tolist(), Y_obs[pareto].tolist()):
        print(f"  {x[0]:6.3f}  {x[1]:6.3f}  {x[2]:6.3f}  {x[3]:6.3f}  {x[4]:6.3f}  {x[5]:6.3f}"
              f"  {y[0]:12.4g}  {1.0 - y[1]:12.4f}")

    print("\n[bo] done.")

    # ---- end-of-run visualisation -----------------------------------------
    aggregate_all()
    print("\n[bo] === Generating visualizations ===")
    vis_ret = subprocess.run([
        sys.executable,
        str(CASE_ROOT / "visualize_results.py"),
        "--results-dir", str(RESULTS_DIR),
        "--output-dir",  str(RESULTS_DIR / "visualizations"),
    ])
    if vis_ret.returncode != 0:
        print("[bo] WARNING: visualization failed (non-fatal)", file=sys.stderr)


if __name__ == "__main__":
    main()
