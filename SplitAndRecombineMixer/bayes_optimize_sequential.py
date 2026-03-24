#!/usr/bin/env python3
"""Sequential multi-objective Bayesian optimisation for the SAR mixer.

Objectives (both minimised):
    J_dp  - pressure drop (Pa)          pdrop_pressure_drop_Pa
    J_mix - intensity of segregation    mixing_intensity_of_segregation

The BO works in a mesh-aware latent parameterisation:
    w_s, t_s, L_c           are sampled directly.
    t_m_ratio, L_s_ratio,
    delta_ratio             are sampled on [0, 1] and mapped into CAD values
                            through interval-safe transforms.

This keeps the CAD generator oblivious to the BO policy while ensuring the
optimiser never proposes thin-feature combinations that violate the current
cfMesh resolution assumptions.
"""

import csv
import subprocess
import sys
from pathlib import Path

import torch
import yaml

torch.set_default_dtype(torch.float64)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CASE_ROOT = Path(__file__).resolve().parent

CFG_PATH = CASE_ROOT / "bayes_optimize_sequential.yaml"
with open(CFG_PATH) as _f:
    CFG = yaml.safe_load(_f)

N_INIT = int(CFG["n_init"])
N_BO = int(CFG["n_bo"])
CORES = int(CFG["cores"])

RESULTS_DIR = CASE_ROOT / "results"
SNAKEFILE = CASE_ROOT / "Snakefile"
TEMPLATE_YAML = CASE_ROOT / "SplitAndRecombineHydro" / "sar_mixer_cad.yaml"
MODEL_PATH = CASE_ROOT / "SplitAndRecombineMixer.pt"

with open(TEMPLATE_YAML) as _f:
    TEMPLATE_GEO = yaml.safe_load(_f)

BO_PARAM_NAMES = list(CFG["bo_parameters"].keys())
GEO_PARAM_NAMES = ["w_s", "t_s", "t_m", "L_s", "L_m", "delta"]

BOUNDS = torch.tensor([
    [float(CFG["bo_parameters"][name]["lower"]) for name in BO_PARAM_NAMES],
    [float(CFG["bo_parameters"][name]["upper"]) for name in BO_PARAM_NAMES],
])

BO_LOWER = {name: float(CFG["bo_parameters"][name]["lower"]) for name in BO_PARAM_NAMES}
BO_UPPER = {name: float(CFG["bo_parameters"][name]["upper"]) for name in BO_PARAM_NAMES}

L_CELL = float(TEMPLATE_GEO["L_cell"])

MESH_FINE_H = float(CFG["mesh_safety"]["fine_cell_size_H"])
T_M_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["t_m_min_cells"])
DELTA_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["delta_min_cells"])
SPLIT_GAP_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["split_gap_min_cells"])
TM_STEP_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["splitter_step_min_cells"])

T_M_MAX = float(CFG["cad_parameter_bounds"]["t_m"]["upper"])
L_S_MIN = float(CFG["cad_parameter_bounds"]["L_s"]["lower"])
L_S_MAX = float(CFG["cad_parameter_bounds"]["L_s"]["upper"])
L_M_MIN = float(CFG["cad_parameter_bounds"]["L_m"]["lower"])
L_M_MAX = float(CFG["cad_parameter_bounds"]["L_m"]["upper"])
DELTA_MAX = float(CFG["cad_parameter_bounds"]["delta"]["upper"])

LC_MIN = BO_LOWER["L_c"]
INEQ_CONSTRAINTS: list = []

OBJ_PDROP = "pdrop_pressure_drop_Pa"
OBJ_MIX = "mixing_intensity_of_segregation"

PENALTY_PDROP = float(CFG["penalties"]["pdrop"])
PENALTY_MIX = float(CFG["penalties"]["mix"])


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def _lerp(lo: float, hi: float, alpha: float) -> float:
    if hi <= lo:
        return lo
    a = _clamp01(alpha)
    return lo + a * (hi - lo)


def _ratio(value: float, lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    return _clamp01((value - lo) / (hi - lo))


def _within(value: float, lo: float, hi: float, tol: float = 1e-12) -> bool:
    return (lo - tol) <= value <= (hi + tol)


def _validate_config() -> None:
    if BO_LOWER["t_s"] < T_M_MIN + TM_STEP_MIN:
        raise ValueError(
            "Invalid BO config: t_s.lower must be >= t_m_min + splitter_step_min."
        )
    if BO_LOWER["w_s"] - 0.5 * BO_UPPER["t_s"] - SPLIT_GAP_MIN < DELTA_MIN:
        raise ValueError(
            "Invalid BO config: delta interval becomes empty at worst-case "
            "w_s/t_s bounds."
        )
    if BO_UPPER["L_c"] > L_CELL - L_S_MIN - L_M_MIN:
        raise ValueError(
            "Invalid BO config: L_c.upper exceeds the interval implied by "
            "L_s/L_m lower bounds."
        )
    if BO_LOWER["L_c"] < (L_CELL - L_S_MAX - L_M_MAX) - 1e-12:
        raise ValueError(
            "Invalid BO config: L_c.lower is too small for the configured "
            "L_s/L_m upper bounds."
        )


def bo_to_geo(bo_params: dict) -> dict:
    """Map BO parameters into actual CAD parameters."""
    w_s = float(bo_params["w_s"])
    t_s = float(bo_params["t_s"])
    t_m_ratio = float(bo_params["t_m_ratio"])
    L_c = float(bo_params["L_c"])
    L_s_ratio = float(bo_params["L_s_ratio"])
    delta_ratio = float(bo_params["delta_ratio"])

    t_m_lo = T_M_MIN
    t_m_hi = min(T_M_MAX, t_s - TM_STEP_MIN)
    if t_m_hi < t_m_lo - 1e-12:
        raise ValueError(
            f"Empty t_m interval for t_s={t_s:.6f}: [{t_m_lo:.6f}, {t_m_hi:.6f}]"
        )
    t_m = _lerp(t_m_lo, t_m_hi, t_m_ratio)

    L_s_lo = max(L_S_MIN, L_CELL - L_c - L_M_MAX)
    L_s_hi = min(L_S_MAX, L_CELL - L_c - L_M_MIN)
    if L_s_hi < L_s_lo - 1e-12:
        raise ValueError(
            f"Empty L_s interval for L_c={L_c:.6f}: [{L_s_lo:.6f}, {L_s_hi:.6f}]"
        )
    L_s = _lerp(L_s_lo, L_s_hi, L_s_ratio)
    L_m = L_CELL - L_c - L_s

    delta_lo = DELTA_MIN
    delta_hi = min(DELTA_MAX, w_s - 0.5 * t_s - SPLIT_GAP_MIN)
    if delta_hi < delta_lo - 1e-12:
        raise ValueError(
            "Empty delta interval for "
            f"w_s={w_s:.6f}, t_s={t_s:.6f}: [{delta_lo:.6f}, {delta_hi:.6f}]"
        )
    delta = _lerp(delta_lo, delta_hi, delta_ratio)

    geo = {
        "w_s": w_s,
        "t_s": t_s,
        "t_m": t_m,
        "L_s": L_s,
        "L_m": L_m,
        "delta": delta,
    }
    return geo


def geo_to_bo(geo_params: dict) -> dict:
    """Map actual CAD parameters back into the latent BO coordinates."""
    w_s = float(geo_params["w_s"])
    t_s = float(geo_params["t_s"])
    t_m = float(geo_params["t_m"])
    L_s = float(geo_params["L_s"])
    L_m = float(geo_params["L_m"])
    delta = float(geo_params["delta"])

    if not _within(w_s, BO_LOWER["w_s"], BO_UPPER["w_s"]):
        raise ValueError("w_s outside current BO bounds")
    if not _within(t_s, BO_LOWER["t_s"], BO_UPPER["t_s"]):
        raise ValueError("t_s outside current BO bounds")

    t_m_lo = T_M_MIN
    t_m_hi = min(T_M_MAX, t_s - TM_STEP_MIN)
    if not _within(t_m, t_m_lo, t_m_hi):
        raise ValueError("t_m outside admissible dependent interval")
    t_m_ratio = _ratio(t_m, t_m_lo, t_m_hi)

    L_c = L_CELL - L_s - L_m
    if not _within(L_c, BO_LOWER["L_c"], BO_UPPER["L_c"]):
        raise ValueError("L_c outside current BO bounds")

    L_s_lo = max(L_S_MIN, L_CELL - L_c - L_M_MAX)
    L_s_hi = min(L_S_MAX, L_CELL - L_c - L_M_MIN)
    if not _within(L_s, L_s_lo, L_s_hi):
        raise ValueError("L_s outside admissible dependent interval")
    L_s_ratio = _ratio(L_s, L_s_lo, L_s_hi)

    delta_lo = DELTA_MIN
    delta_hi = min(DELTA_MAX, w_s - 0.5 * t_s - SPLIT_GAP_MIN)
    if not _within(delta, delta_lo, delta_hi):
        raise ValueError("delta outside admissible dependent interval")
    delta_ratio = _ratio(delta, delta_lo, delta_hi)

    bo = {
        "w_s": w_s,
        "t_s": t_s,
        "t_m_ratio": t_m_ratio,
        "L_c": L_c,
        "L_s_ratio": L_s_ratio,
        "delta_ratio": delta_ratio,
    }
    return bo


def is_feasible(bo_params: dict) -> bool:
    """Return True when the transformed geometry stays mesh-safe."""
    try:
        geo = bo_to_geo(bo_params)
    except ValueError:
        return False

    w_s = geo["w_s"]
    t_s = geo["t_s"]
    t_m = geo["t_m"]
    L_s = geo["L_s"]
    L_m = geo["L_m"]
    delta = geo["delta"]

    c1 = w_s - 0.5 * t_s - delta >= SPLIT_GAP_MIN - 1e-12
    c2 = L_s + L_m <= L_CELL - LC_MIN + 1e-12
    c3 = t_s - t_m >= TM_STEP_MIN - 1e-12
    c4 = t_m >= T_M_MIN - 1e-12
    c5 = delta >= DELTA_MIN - 1e-12
    return c1 and c2 and c3 and c4 and c5


def _annotation_fields(bo_params: dict, geo_params: dict) -> dict:
    return {
        **{f"bo_{k}": bo_params[k] for k in BO_PARAM_NAMES},
        **{f"geo_{k}": geo_params[k] for k in GEO_PARAM_NAMES},
    }


def _annotate_objectives_csv(sample_dir: Path, bo_params: dict, geo_params: dict) -> None:
    """Add BO and geometry coordinates to a generated objectives.csv."""
    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        return

    with open(obj_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    extra = _annotation_fields(bo_params, geo_params)
    for key in extra:
        if key not in fieldnames:
            fieldnames.append(key)
    for row in rows:
        row.update(extra)

    with open(obj_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


_validate_config()


# ---------------------------------------------------------------------------
# Pareto front plot
# ---------------------------------------------------------------------------

def plot_pareto_front(X: torch.Tensor, Y: torch.Tensor, n_init: int, title: str) -> None:
    """Save RESULTS_DIR/pareto_front.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from botorch.utils.multi_objective.pareto import is_non_dominated

    j_dp = Y[:, 0].numpy()
    mq = 1.0 - Y[:, 1].numpy()

    n = len(j_dp)
    colors = ["steelblue" if i < n_init else "darkorange" for i in range(n)]
    pareto = is_non_dominated(-Y).numpy()

    fig, ax = plt.subplots(figsize=(7, 5))

    if (~pareto).any():
        ax.scatter(
            j_dp[~pareto],
            mq[~pareto],
            c=[colors[i] for i in range(n) if not pareto[i]],
            alpha=0.5,
            s=40,
        )

    if pareto.any():
        ax.scatter(
            j_dp[pareto],
            mq[pareto],
            c=[colors[i] for i in range(n) if pareto[i]],
            s=90,
            edgecolors="black",
            linewidths=1.4,
            zorder=5,
        )

    if pareto.sum() > 1:
        px, py = j_dp[pareto], mq[pareto]
        order = px.argsort()
        ax.step(px[order], py[order], where="post", color="black", lw=1.5, alpha=0.8)

    ax.set_xlabel(r"Pressure drop  $J_\mathrm{dp}$  (Pa, log scale)", fontsize=11)
    ax.set_ylabel(r"Mixing quality  $1 - I_s$  (-)", fontsize=11)
    ax.set_xscale("log")
    ax.set_title(title, fontsize=11)
    ax.grid(True, which="both", alpha=0.3)

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue",
               markersize=9, label=f"Sobol init ({min(n_init, n)})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="darkorange",
               markersize=9, label=f"BO suggested ({max(0, n - n_init)})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor="black", markeredgewidth=1.4, markersize=9,
               label=f"Pareto ({int(pareto.sum())})"),
        Line2D([0], [0], color="black", lw=1.5, label="Pareto front"),
    ]
    ax.legend(handles=handles, fontsize=9, loc="lower left")

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "pareto_front.png", dpi=150)
    plt.close(fig)
    print(f"[bo] plot saved -> {RESULTS_DIR / 'pareto_front.png'}")


# ---------------------------------------------------------------------------
# Sample execution
# ---------------------------------------------------------------------------

def aggregate_all() -> None:
    """Rewrite results/all_samples.csv from all completed objectives.csv files."""
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
        writer = csv.DictWriter(
            fh, fieldnames=all_fieldnames, extrasaction="ignore", restval=""
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"[bo] all_samples.csv -> {len(rows)} sample(s)  ({out})")


def next_sample_id() -> str:
    ids = (
        [int(d.name) for d in RESULTS_DIR.iterdir() if d.is_dir() and d.name.isdigit()]
        if RESULTS_DIR.exists()
        else []
    )
    return f"{(max(ids) + 1 if ids else 0):05d}"


def _write_penalty_objectives(sample_dir: Path, bo_params: dict, geo_params: dict) -> None:
    """Write a minimal objectives.csv with penalty values for a failed sample."""
    row = {
        "sample_id": sample_dir.name,
        "results_dir": str(sample_dir.resolve()),
        "failed": "True",
        **_annotation_fields(bo_params, geo_params),
        OBJ_PDROP: PENALTY_PDROP,
        OBJ_MIX: PENALTY_MIX,
    }
    out = sample_dir / "objectives.csv"
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    print(f"[bo] penalty objectives.csv written -> {out}")


def evaluate(bo_params: dict) -> tuple[float, float] | tuple[None, None]:
    """Write YAML, run Snakemake, return (J_dp, J_mix) or (None, None)."""
    if not is_feasible(bo_params):
        print(f"[bo] SKIP infeasible params: {bo_params}", file=sys.stderr)
        return None, None

    geo_params = bo_to_geo(bo_params)
    sid = next_sample_id()
    sample_dir = RESULTS_DIR / sid
    sample_dir.mkdir(parents=True, exist_ok=True)

    geo_yaml = dict(TEMPLATE_GEO)
    geo_yaml.update(geo_params)
    with open(sample_dir / "sar_mixer_cad.yaml", "w") as fh:
        yaml.dump(geo_yaml, fh, default_flow_style=False, sort_keys=False)

    print(f"\n[bo] sample {sid}: bo={bo_params}")
    print(f"[bo] sample {sid}: geo={geo_params}")

    ret = subprocess.run(
        [
            "snakemake",
            "--snakefile", str(SNAKEFILE),
            "--directory", str(sample_dir),
            "--cores", str(CORES),
            "--config", f"results_dir={sample_dir}",
        ]
    )
    if ret.returncode != 0:
        print(f"[bo] WARNING: sample {sid} failed - penalising", file=sys.stderr)
        _write_penalty_objectives(sample_dir, bo_params, geo_params)
        return PENALTY_PDROP, PENALTY_MIX

    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        print(f"[bo] WARNING: objectives.csv missing for {sid} - penalising",
              file=sys.stderr)
        _write_penalty_objectives(sample_dir, bo_params, geo_params)
        return PENALTY_PDROP, PENALTY_MIX

    _annotate_objectives_csv(sample_dir, bo_params, geo_params)

    with open(obj_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            j_dp = float(row[OBJ_PDROP])
            j_mix = float(row[OBJ_MIX])
            print(f"[bo] {sid}: J_dp={j_dp:.4g} Pa  mixing quality={1.0 - j_mix:.4f}")
            return j_dp, j_mix

    return None, None


# ---------------------------------------------------------------------------
# BO helpers
# ---------------------------------------------------------------------------

def fit_model(X: torch.Tensor, Y: torch.Tensor, warm_start: dict | None = None):
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.mlls import ExactMarginalLogLikelihood

    model = SingleTaskGP(
        X,
        -Y,
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
            print(f"[bo] WARNING: warm start failed ({exc}) - using default init",
                  file=sys.stderr)
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model))
    return model


def save_model(model) -> None:
    """Persist GP hyperparameters to MODEL_PATH for warm-starting future runs."""
    checkpoint = {k: v.detach().clone() for k, v in model.named_parameters()}
    torch.save(checkpoint, MODEL_PATH)
    print(f"[bo] model checkpoint saved -> {MODEL_PATH}")


def load_model_hyperparams() -> dict | None:
    """Load GP hyperparameters from MODEL_PATH, or return None if absent/corrupt."""
    if not MODEL_PATH.exists():
        return None
    try:
        ckpt = torch.load(MODEL_PATH, weights_only=True)
        print(f"[bo] model checkpoint loaded <- {MODEL_PATH}")
        return ckpt
    except Exception as exc:
        print(f"[bo] WARNING: could not load model checkpoint ({exc}) - fitting from scratch",
              file=sys.stderr)
        return None


def next_candidate(model, X: torch.Tensor, Y: torch.Tensor) -> torch.Tensor:
    from botorch.acquisition.multi_objective import qNoisyExpectedHypervolumeImprovement
    from botorch.optim import optimize_acqf

    Y_neg = -Y
    Y_range = (Y_neg.max(0).values - Y_neg.min(0).values).clamp(min=1e-6)
    ref_point = Y_neg.min(0).values - 0.1 * Y_range

    acqf = qNoisyExpectedHypervolumeImprovement(
        model=model, ref_point=ref_point, X_baseline=X, prune_baseline=True
    )

    opt_kwargs = {
        "bounds": BOUNDS,
        "q": 1,
        "num_restarts": 10,
        "raw_samples": 256,
    }
    if INEQ_CONSTRAINTS:
        opt_kwargs["inequality_constraints"] = INEQ_CONSTRAINTS
    cand, _ = optimize_acqf(acqf, **opt_kwargs)
    return cand.squeeze(0)


# ---------------------------------------------------------------------------
# Resume helpers
# ---------------------------------------------------------------------------

def collect_existing() -> tuple:
    """Load all completed samples from RESULTS_DIR.

    New runs store bo_<name> and geo_<name> columns. Older result folders may
    only contain geo_<name>; those are projected into the new latent variables
    when possible. Samples that fall outside the new mesh-safe admissible box
    are skipped.
    """
    if not RESULTS_DIR.exists():
        return None, None, 0

    xs, ys = [], []
    skipped_legacy = 0
    for d in sorted(RESULTS_DIR.iterdir()):
        if not (d.is_dir() and d.name.isdigit()):
            continue
        obj_csv = d / "objectives.csv"
        if not obj_csv.exists():
            continue

        with open(obj_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    if all(f"bo_{k}" in row and row[f"bo_{k}"] != "" for k in BO_PARAM_NAMES):
                        bo = {k: float(row[f"bo_{k}"]) for k in BO_PARAM_NAMES}
                    else:
                        geo = {k: float(row[f"geo_{k}"]) for k in GEO_PARAM_NAMES}
                        bo = geo_to_bo(geo)
                    j_dp = float(row[OBJ_PDROP])
                    j_mix = float(row[OBJ_MIX])
                except (KeyError, ValueError):
                    continue
                except Exception:
                    skipped_legacy += 1
                    continue

                if not is_feasible(bo):
                    skipped_legacy += 1
                    continue
                xs.append([float(bo[k]) for k in BO_PARAM_NAMES])
                ys.append([j_dp, j_mix])

    if skipped_legacy > 0:
        print(
            f"[bo] WARNING: skipped {skipped_legacy} legacy/incompatible sample(s) "
            "outside the current mesh-safe BO parameterisation",
            file=sys.stderr,
        )

    if not xs:
        return None, None, 0
    return torch.tensor(xs), torch.tensor(ys), len(xs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    X_obs, Y_obs, n_existing = collect_existing()
    n_init_done = min(n_existing, N_INIT)
    if n_existing > 0:
        n_bo_existing = max(0, n_existing - N_INIT)
        print(f"[bo] resuming: {n_existing} completed sample(s) found "
              f"({n_init_done} init, {n_bo_existing} BO)")
        aggregate_all()
        plot_pareto_front(X_obs, Y_obs, n_init_done,
                          f"SAR mixer - resumed  [{n_existing} sample(s)]")

    n_init_needed = max(0, N_INIT - n_existing)
    if n_init_needed > 0:
        from botorch.utils.sampling import draw_sobol_samples

        print(f"[bo] === Sobol initialisation: {n_init_needed} remaining of {N_INIT} ===")
        sobol_X = draw_sobol_samples(bounds=BOUNDS, n=n_init_needed, q=1).squeeze(1)

        for x in sobol_X:
            bo_params = {k: float(x[j]) for j, k in enumerate(BO_PARAM_NAMES)}
            j_dp, j_mix = evaluate(bo_params)
            if j_dp is None:
                continue

            x_t = torch.tensor([[float(bo_params[k]) for k in BO_PARAM_NAMES]])
            y_t = torch.tensor([[j_dp, j_mix]])
            X_obs = x_t if X_obs is None else torch.cat([X_obs, x_t])
            Y_obs = y_t if Y_obs is None else torch.cat([Y_obs, y_t])
            n_init_done += 1

            aggregate_all()
            plot_pareto_front(
                X_obs,
                Y_obs,
                n_init_done,
                f"SAR mixer - init {n_init_done}/{N_INIT}  [{X_obs.shape[0]} sample(s)]",
            )

    if X_obs is None or X_obs.shape[0] < 2:
        sys.exit("[bo] not enough successful samples to fit a GP - aborting")

    n_bo_done = max(0, n_existing - N_INIT)
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
        bo_params = {k: float(x_next[j]) for j, k in enumerate(BO_PARAM_NAMES)}

        j_dp, j_mix = evaluate(bo_params)
        if j_dp is None:
            continue

        x_t = torch.tensor([[float(bo_params[k]) for k in BO_PARAM_NAMES]])
        y_t = torch.tensor([[j_dp, j_mix]])
        X_obs = torch.cat([X_obs, x_t])
        Y_obs = torch.cat([Y_obs, y_t])

        aggregate_all()
        plot_pareto_front(
            X_obs,
            Y_obs,
            n_init_done,
            f"SAR mixer - BO iter {bo_iter}/{N_BO}  [{X_obs.shape[0]} sample(s)]",
        )

    from botorch.utils.multi_objective.pareto import is_non_dominated

    pareto = is_non_dominated(-Y_obs)
    print("\n[bo] Pareto-optimal designs:")
    print(f"  {'w_s':>6}  {'t_s':>6}  {'t_m':>6}  {'L_s':>6}  {'L_m':>6}  {'delta':>6}"
          f"  {'J_dp [Pa]':>12}  {'mix quality':>12}")
    for x, y in zip(X_obs[pareto].tolist(), Y_obs[pareto].tolist()):
        bo_params = {k: float(x[j]) for j, k in enumerate(BO_PARAM_NAMES)}
        geo = bo_to_geo(bo_params)
        print(
            f"  {geo['w_s']:6.3f}  {geo['t_s']:6.3f}  {geo['t_m']:6.3f}"
            f"  {geo['L_s']:6.3f}  {geo['L_m']:6.3f}  {geo['delta']:6.3f}"
            f"  {y[0]:12.4g}  {1.0 - y[1]:12.4f}"
        )

    print("\n[bo] done.")

    aggregate_all()
    print("\n[bo] === Generating visualizations ===")
    vis_ret = subprocess.run(
        [
            sys.executable,
            str(CASE_ROOT / "visualize_results.py"),
            "--results-dir", str(RESULTS_DIR),
            "--output-dir", str(RESULTS_DIR / "visualizations"),
        ]
    )
    if vis_ret.returncode != 0:
        print("[bo] WARNING: visualization failed (non-fatal)", file=sys.stderr)


if __name__ == "__main__":
    main()
