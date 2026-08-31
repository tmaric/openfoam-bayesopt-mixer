#!/usr/bin/env python3
"""Sequential multi-objective BO for the planar alternating-deflector mixer.

Objectives (both minimised):
    J_dp  - pressure ratio relative to the validated straight channel
    J_mix - flux-weighted segregation intensity

The raw kinematic and dimensional pressure drops, flow rate, pumping power and
literature-style mixing index are retained in every objective record.  The
default stage executes only the corrected 12-point feasibility screen.  Full
BO must be requested explicitly and is guarded by the predeclared screen.

The BO works in a mesh-aware physical parameterisation. The weak-wall cosine
amplitude is sampled directly and a_strong_ratio maps to an admissible strong
amplitude. This replaces the correlated w_s/delta coordinates used by the
pilot campaign. Splitter and section lengths retain interval-safe transforms.

This keeps the CAD generator oblivious to the BO policy while ensuring the
optimiser never proposes thin-feature combinations that violate the current
cfMesh resolution assumptions.
"""

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import torch
import yaml

import padm_runner

torch.set_default_dtype(torch.float64)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CASE_ROOT = Path(__file__).resolve().parent

CFG_PATH = CASE_ROOT / "bayes_optimize_sequential.yaml"
with open(CFG_PATH) as _f:
    CFG = yaml.safe_load(_f)

SCREENING_N_INIT = int(CFG["screening_n_init"])
N_INIT = int(CFG["n_init"])
N_BO = int(CFG["n_bo"])
# MPI ranks per CFD solve.  `cores` is the deprecated spelling, kept so the
# archived campaign configs still load unchanged.
NP = int(CFG.get("np", CFG.get("cores", 2)))
Q_BATCH = int(CFG.get("q", 1))
SOBOL_SEED = int(CFG.get("sobol_seed", 0))
TORCH_THREADS = int(CFG.get("torch_threads", 1))
torch.set_num_threads(TORCH_THREADS)

_results_path = Path(CFG.get("results_dir", "results"))
RESULTS_DIR = (_results_path if _results_path.is_absolute() else CASE_ROOT / _results_path).resolve()
SNAKEFILE = CASE_ROOT / "Snakefile"
# Upper bound on ranks.  At ~150k cells the useful ceiling is far below this;
# the limit exists to catch a typo, not to express a physical optimum.  Keep np
# EQUAL across the designs of one campaign: an MPI job runs at the pace of its
# slowest rank, so an np that varies between designs makes them incomparable.
MAX_NP = 16
# Resolved once by main() from --profile / $PADM_SNAKEMAKE_PROFILE.  Module level
# because evaluate() is reached from both the Sobol and the acquisition loop.
PROFILE_DIR: Path | None = None
# Catches the failure mode the per-design classifier structurally cannot: a
# broken mesher or solver fails a legitimate PHYSICS rule identically on every
# design, and would otherwise be recorded as "the whole space is infeasible".
FAILURE_STREAK = padm_runner.FailureStreak()
CAD_CONFIG_NAME = "alternating_deflector_cad.yaml"
TEMPLATE_YAML = CASE_ROOT / "FlowCase" / CAD_CONFIG_NAME
MODEL_PATH = RESULTS_DIR / str(CFG.get("model_file", "gp_checkpoint.pt"))
_baseline_path = Path(CFG["baseline_summary"])
BASELINE_SUMMARY = (
    _baseline_path if _baseline_path.is_absolute() else CASE_ROOT / _baseline_path
).resolve()
SCREENING_GATE_PATH = RESULTS_DIR / "screening_gate.json"

with open(TEMPLATE_YAML) as _f:
    TEMPLATE_GEO = yaml.safe_load(_f)

BO_PARAM_NAMES = list(CFG["bo_parameters"].keys())
GEO_PARAM_NAMES = [
    "w_s", "t_s", "t_m", "L_s", "L_m", "delta", "k",
    "a_weak", "a_strong",
]

BOUNDS = torch.tensor([
    [float(CFG["bo_parameters"][name]["lower"]) for name in BO_PARAM_NAMES],
    [float(CFG["bo_parameters"][name]["upper"]) for name in BO_PARAM_NAMES],
])

BO_LOWER = {name: float(CFG["bo_parameters"][name]["lower"]) for name in BO_PARAM_NAMES}
BO_UPPER = {name: float(CFG["bo_parameters"][name]["upper"]) for name in BO_PARAM_NAMES}

L_CELL = float(TEMPLATE_GEO["L_cell"])
L0 = float(TEMPLATE_GEO["L0"])
N_CELLS = int(TEMPLATE_GEO["N"])
TOTAL_L = 2.0 * L0 + N_CELLS * L_CELL

MESH_FINE_H = float(CFG["mesh_safety"]["fine_cell_size_H"])
T_M_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["t_m_min_cells"])
SPLIT_GAP_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["split_gap_min_cells"])
TM_STEP_MIN = MESH_FINE_H * float(CFG["mesh_safety"]["splitter_step_min_cells"])

T_M_MAX = float(CFG["cad_parameter_bounds"]["t_m"]["upper"])
A_STRONG_MIN = float(CFG["cad_parameter_bounds"]["a_strong"]["lower"])
A_STRONG_MAX = float(CFG["cad_parameter_bounds"]["a_strong"]["upper"])
STRONG_WEAK_CONTRAST_MIN = float(
    CFG["cad_parameter_bounds"]["strong_weak_contrast_min"]
)
L_S_MIN = float(CFG["cad_parameter_bounds"]["L_s"]["lower"])
L_S_MAX = float(CFG["cad_parameter_bounds"]["L_s"]["upper"])
L_M_MIN = float(CFG["cad_parameter_bounds"]["L_m"]["lower"])
L_M_MAX = float(CFG["cad_parameter_bounds"]["L_m"]["upper"])
LC_MIN = BO_LOWER["L_c"]
INEQ_CONSTRAINTS: list = []

OBJ_PDROP_RAW = "pdrop_pressure_drop_m2_s2"
OBJ_PDROP_LEGACY = "pdrop_pressure_drop_Pa"
OBJ_PRESSURE_PA = "metric_pressure_drop_Pa"
OBJ_PRESSURE_RATIO = "metric_pressure_ratio_to_straight"
OBJ_MIX = "mixing_flux_weighted_intensity_of_segregation"

ACQ_NUM_RESTARTS = int(CFG["acquisition"]["num_restarts"])
ACQ_RAW_SAMPLES = int(CFG["acquisition"]["raw_samples"])
REF_PRESSURE_PA = float(CFG["reference_point"]["pressure_drop_Pa"])
REF_MIX = float(CFG["reference_point"]["flux_intensity_of_segregation"])
SCREENING_MIN_MIXING_INDEX = float(CFG["screening_gate"]["minimum_mixing_index"])
SCREENING_MAX_FAILED_FRACTION = float(CFG["screening_gate"]["maximum_failed_fraction"])


def _clamp01(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def _lerp(lo: float, hi: float, alpha: float) -> float:
    if hi - lo <= 1e-12:
        return lo
    a = _clamp01(alpha)
    return lo + a * (hi - lo)


def _ratio(value: float, lo: float, hi: float) -> float:
    if hi - lo <= 1e-12:
        return 0.0
    return _clamp01((value - lo) / (hi - lo))


def _within(value: float, lo: float, hi: float, tol: float = 1e-12) -> bool:
    return (lo - tol) <= value <= (hi + tol)


def _validate_config() -> None:
    try:
        RESULTS_DIR.relative_to(CASE_ROOT)
    except ValueError as exc:
        raise ValueError("results_dir must remain inside the study directory") from exc
    try:
        BASELINE_SUMMARY.relative_to(CASE_ROOT)
    except ValueError as exc:
        raise ValueError("baseline_summary must remain inside the study directory") from exc
    if not (2 <= SCREENING_N_INIT <= N_INIT):
        raise ValueError("screening_n_init must be between 2 and n_init")
    if Q_BATCH != 1:
        raise ValueError("This research campaign requires strictly sequential BO (q = 1).")
    if NP < 1 or NP > MAX_NP:
        raise ValueError(f"np must be between 1 and {MAX_NP}.")
    if TORCH_THREADS != 1:
        raise ValueError("torch_threads must be 1 for the resource-limited campaign.")
    if N_CELLS < 1:
        raise ValueError("Invalid CAD template: N must be >= 1.")
    if BO_LOWER["t_s"] < T_M_MIN + TM_STEP_MIN:
        raise ValueError(
            "Invalid BO config: t_s.lower must be >= t_m_min + splitter_step_min."
        )
    if A_STRONG_MAX <= A_STRONG_MIN:
        raise ValueError("a_strong upper bound must exceed its lower bound.")
    if A_STRONG_MAX <= BO_UPPER["a_weak"] + STRONG_WEAK_CONTRAST_MIN:
        raise ValueError("a_strong bounds cannot satisfy the maximum weak amplitude.")
    if 1.0 - BO_UPPER["a_weak"] - A_STRONG_MAX - 4.0 * 0.01 < SPLIT_GAP_MIN:
        raise ValueError("worst-case strong/weak deflectors leave too little peak gap.")
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
    if ACQ_NUM_RESTARTS < 1 or ACQ_RAW_SAMPLES < ACQ_NUM_RESTARTS:
        raise ValueError("invalid acquisition restart/raw-sample configuration")
    if REF_PRESSURE_PA <= 0 or REF_MIX <= 0:
        raise ValueError("reference-point objectives must be positive")
    if not (0.0 <= SCREENING_MIN_MIXING_INDEX <= 1.0):
        raise ValueError("screening minimum_mixing_index must lie in [0, 1]")
    if not (0.0 <= SCREENING_MAX_FAILED_FRACTION < 1.0):
        raise ValueError("screening maximum_failed_fraction must lie in [0, 1)")
    transforms = CFG.get("objective_transforms", {})
    if any(value != "linear" for value in transforms.values()):
        raise ValueError(
            "objective transforms remain linear until corrected screening data "
            "justify and test an explicit model-space transform"
        )


def _interaction_midpoint_xhat(L_s: float, L_m: float, cell_idx: int) -> float:
    L_c = L_CELL - L_s - L_m
    x_mid = L0 + cell_idx * L_CELL + L_s + 0.5 * L_c
    return x_mid / TOTAL_L


def _cell_xhats(L_s: float, L_m: float) -> list[float]:
    return [_interaction_midpoint_xhat(L_s, L_m, cell_idx) for cell_idx in range(N_CELLS)]


def _realized_deltas(delta: float, k: float, L_s: float, L_m: float) -> list[float]:
    return [delta + k * xhat for xhat in _cell_xhats(L_s, L_m)]


def _realized_geometry_metrics(w_s: float, delta: float, k: float, L_s: float, L_m: float) -> dict:
    deltas = _realized_deltas(delta, k, L_s, L_m)
    delta_first = deltas[0]
    delta_last = deltas[-1]
    delta_min = min(deltas)
    delta_max = max(deltas)
    h_d = 0.5 - w_s
    return {
        "realized_deltas": deltas,
        "delta_first": delta_first,
        "delta_last": delta_last,
        "delta_min": delta_min,
        "delta_max": delta_max,
        "delta_span": delta_max - delta_min,
        "ramp_up": delta_last - delta_first,
        "peak_intrusion_max": h_d + delta_max,
        "centerline_margin_min": w_s - delta_max,
    }


def _format_realized_geometry_metrics(metrics: dict) -> str:
    return (
        f"delta_first={metrics['delta_first']:.4f}, "
        f"delta_last={metrics['delta_last']:.4f}, "
        f"delta_span={metrics['delta_span']:.4f}, "
        f"peak_intrusion_max={metrics['peak_intrusion_max']:.4f}, "
        f"centerline_margin_min={metrics['centerline_margin_min']:.4f}"
    )


def _is_geo_feasible(geo: dict) -> bool:
    """Return True when the realized geometry stays mesh-safe."""
    a_weak = float(geo["a_weak"])
    a_strong = float(geo["a_strong"])
    t_s = float(geo["t_s"])
    t_m = float(geo["t_m"])
    L_s = float(geo["L_s"])
    L_m = float(geo["L_m"])
    k = float(geo["k"])

    peak_gap = 1.0 - a_weak - a_strong - 4.0 * 0.01
    split_gap = 0.5 - a_weak - 0.5 * t_s
    checks = (
        _within(a_weak, BO_LOWER["a_weak"], BO_UPPER["a_weak"]),
        _within(a_strong, max(A_STRONG_MIN, a_weak + STRONG_WEAK_CONTRAST_MIN), A_STRONG_MAX),
        abs(k) <= 1e-12,
        peak_gap >= SPLIT_GAP_MIN - 1e-12,
        split_gap >= SPLIT_GAP_MIN - 1e-12,
        L_s + L_m <= L_CELL - LC_MIN + 1e-12,
        t_s - t_m >= TM_STEP_MIN - 1e-12,
        t_m >= T_M_MIN - 1e-12,
    )
    return all(checks)


def bo_to_geo(bo_params: dict) -> dict:
    """Map BO parameters into actual CAD parameters."""
    a_weak = float(bo_params["a_weak"])
    a_strong_ratio = float(bo_params["a_strong_ratio"])
    t_s = float(bo_params["t_s"])
    t_m_ratio = float(bo_params["t_m_ratio"])
    L_c = float(bo_params["L_c"])
    L_s_ratio = float(bo_params["L_s_ratio"])

    a_strong_lo = max(A_STRONG_MIN, a_weak + STRONG_WEAK_CONTRAST_MIN)
    if a_strong_lo > A_STRONG_MAX + 1e-12:
        raise ValueError(
            f"Empty strong-amplitude interval for a_weak={a_weak:.6f}"
        )
    a_strong = _lerp(a_strong_lo, A_STRONG_MAX, a_strong_ratio)

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

    w_s = 0.5 - a_weak
    delta = a_strong - a_weak
    k = 0.0

    geo = {
        "w_s": w_s,
        "t_s": t_s,
        "t_m": t_m,
        "L_s": L_s,
        "L_m": L_m,
        "delta": delta,
        "k": k,
        "a_weak": a_weak,
        "a_strong": a_strong,
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
    k = float(geo_params.get("k", 0.0))
    a_weak = float(geo_params.get("a_weak", 0.5 - w_s))
    a_strong = float(geo_params.get("a_strong", a_weak + delta))

    if abs(k) > 1e-12:
        raise ValueError("the verified campaign uses a constant strong amplitude (k=0)")
    if not _within(a_weak, BO_LOWER["a_weak"], BO_UPPER["a_weak"]):
        raise ValueError("a_weak outside current BO bounds")
    a_strong_lo = max(A_STRONG_MIN, a_weak + STRONG_WEAK_CONTRAST_MIN)
    if not _within(a_strong, a_strong_lo, A_STRONG_MAX):
        raise ValueError("a_strong outside admissible dependent interval")
    a_strong_ratio = _ratio(a_strong, a_strong_lo, A_STRONG_MAX)
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

    bo = {
        "a_weak": a_weak,
        "a_strong_ratio": a_strong_ratio,
        "t_s": t_s,
        "t_m_ratio": t_m_ratio,
        "L_c": L_c,
        "L_s_ratio": L_s_ratio,
    }
    return bo


def is_feasible(bo_params: dict) -> bool:
    """Return True when the transformed geometry stays mesh-safe."""
    try:
        geo = bo_to_geo(bo_params)
    except ValueError:
        return False
    return _is_geo_feasible(geo)


def _annotation_fields(bo_params: dict, geo_params: dict) -> dict:
    return {
        **{f"bo_{k}": bo_params[k] for k in BO_PARAM_NAMES},
        **{f"geo_{k}": geo_params[k] for k in GEO_PARAM_NAMES},
    }


def _straight_baseline() -> dict:
    """Return the validated straight-channel row required for normalization."""
    if not BASELINE_SUMMARY.exists():
        raise FileNotFoundError(
            "validated straight-channel baseline is missing. Run "
            "`python run_baselines.py --max-new-evaluations 1` until "
            f"{BASELINE_SUMMARY} exists before launching corrected BO."
        )
    with open(BASELINE_SUMMARY, newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"straight", "symmetric_deflector", "strong_alternating"}
    completed = {row.get("baseline") for row in rows}
    missing = required.difference(completed)
    if missing:
        raise ValueError(
            "complete all corrected baselines before BO; missing "
            + ", ".join(sorted(missing))
        )
    for row in rows:
        if row.get("baseline") != "straight":
            continue
        if str(row.get("analytical_check_passed", "")).lower() != "true":
            raise ValueError(
                "straight-channel baseline has not passed its analytical pressure check"
            )
        pressure_pa = float(row["pressure_drop_Pa"])
        if pressure_pa <= 0.0:
            raise ValueError("straight-channel pressure drop must be positive")
        return {**row, "pressure_drop_Pa": pressure_pa}
    raise ValueError("baseline summary contains no validated straight-channel row")


def _add_pressure_normalization(row: dict) -> dict:
    """Add dimensional pressure and straight-channel-normalized pressure."""
    normalized = dict(row)
    baseline = _straight_baseline()
    pressure_pa_raw = normalized.get(OBJ_PRESSURE_PA, "")
    if pressure_pa_raw in ("", None):
        # The workflow writes dimensional pressure for all corrected samples.
        # This fallback only supports manually agglomerated corrected rows.
        raw = normalized.get(OBJ_PDROP_RAW, "")
        if raw in ("", None):
            legacy = normalized.get(OBJ_PDROP_LEGACY, "")
            raw = legacy
        if raw not in ("", None):
            research = yaml.safe_load((CASE_ROOT / "research_config.yaml").read_text())
            density = float(research["operating_point"]["fluid_density_kg_m3"])
            normalized[OBJ_PRESSURE_PA] = density * float(raw)
    if normalized.get(OBJ_PRESSURE_PA, "") not in ("", None):
        normalized[OBJ_PRESSURE_RATIO] = (
            float(normalized[OBJ_PRESSURE_PA]) / baseline["pressure_drop_Pa"]
        )
    return normalized


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
    updated_rows = []
    for row in rows:
        row.update(extra)
        updated_rows.append(_add_pressure_normalization(row))
    for row in updated_rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with open(obj_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(updated_rows)


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

    pressure_ratio = Y[:, 0].numpy()
    # Literature-compatible relative-standard-deviation mixing index. The BO
    # still minimises I_s; this monotone transform is used only for reporting.
    mq = 1.0 - Y[:, 1].clamp(min=0.0).sqrt().numpy()

    n = len(pressure_ratio)
    colors = ["steelblue" if i < n_init else "darkorange" for i in range(n)]
    pareto = is_non_dominated(-Y).numpy()

    fig, ax = plt.subplots(figsize=(7, 5))

    if (~pareto).any():
        ax.scatter(
            pressure_ratio[~pareto],
            mq[~pareto],
            c=[colors[i] for i in range(n) if not pareto[i]],
            alpha=0.5,
            s=40,
        )

    if pareto.any():
        ax.scatter(
            pressure_ratio[pareto],
            mq[pareto],
            c=[colors[i] for i in range(n) if pareto[i]],
            s=90,
            edgecolors="black",
            linewidths=1.4,
            zorder=5,
        )

    if pareto.sum() > 1:
        px, py = pressure_ratio[pareto], mq[pareto]
        order = px.argsort()
        ax.step(px[order], py[order], where="post", color="black", lw=1.5, alpha=0.8)

    ax.set_xlabel(
        r"Pressure ratio  $\Delta p/\Delta p_\mathrm{straight}$  (-, log scale)",
        fontsize=11,
    )
    ax.set_ylabel(r"Flux-weighted mixing index  $1 - \sqrt{I_s^\phi}$  (-)", fontsize=11)
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
    expected_fields = [
        "sample_id",
        "results_dir",
        *[f"bo_{name}" for name in BO_PARAM_NAMES],
        *[f"geo_{name}" for name in GEO_PARAM_NAMES],
    ]
    rows = []
    all_fieldnames = list(expected_fields)
    seen_fields: set = set(expected_fields)
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
            for row in reader:
                row = _normalize_objectives_row(row)
                for field in row:
                    if field not in seen_fields:
                        all_fieldnames.append(field)
                        seen_fields.add(field)
                rows.append(row)
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


def _write_failure_objectives(sample_dir: Path, bo_params: dict, geo_params: dict) -> None:
    """Record a failed evaluation without introducing fictitious GP targets."""
    row = {
        "sample_id": sample_dir.name,
        "results_dir": sample_dir.name,
        "failed": "True",
        **_annotation_fields(bo_params, geo_params),
        OBJ_PDROP_RAW: "",
        OBJ_PRESSURE_PA: "",
        OBJ_PRESSURE_RATIO: "",
        OBJ_MIX: "",
    }
    out = sample_dir / "objectives.csv"
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)
    print(f"[bo] failure record written -> {out}")


def _extract_geo_from_row(row: dict) -> dict | None:
    geo = {}
    for name in GEO_PARAM_NAMES:
        raw = row.get(f"geo_{name}", "")
        if raw in ("", None):
            if name == "k":
                raw = 0.0
            else:
                return None
        try:
            geo[name] = float(raw)
        except (TypeError, ValueError):
            return None
    return geo


def _normalize_objectives_row(row: dict) -> dict:
    normalized = dict(row)
    if normalized.get("sample_id", "") not in ("", None):
        normalized["results_dir"] = normalized["sample_id"]
    if normalized.get(OBJ_PDROP_RAW, "") in ("", None):
        legacy_pdrop = normalized.get(OBJ_PDROP_LEGACY, "")
        if legacy_pdrop not in ("", None):
            normalized[OBJ_PDROP_RAW] = legacy_pdrop
    if normalized.get("geo_k", "") in ("", None):
        normalized["geo_k"] = "0.0"

    # Geometry is the physical source of truth. Always reproject it into the
    # current latent coordinates so archived samples remain consistent after a
    # constraint-preserving transform is tightened.
    geo = _extract_geo_from_row(normalized)
    if geo is not None:
        try:
            bo = geo_to_bo(geo)
        except Exception:
            bo = None
        if bo is not None:
            for name, value in bo.items():
                normalized[f"bo_{name}"] = f"{value:.17g}"

    if normalized.get(OBJ_PDROP_RAW, "") not in ("", None):
        normalized = _add_pressure_normalization(normalized)
    return normalized


def _bo_from_row(row: dict) -> dict:
    """Recover current BO coordinates from an objective row."""
    if all(f"bo_{name}" in row and row[f"bo_{name}"] != "" for name in BO_PARAM_NAMES):
        return {name: float(row[f"bo_{name}"]) for name in BO_PARAM_NAMES}

    geo = _extract_geo_from_row(row)
    if geo is None:
        raise KeyError("missing geometry columns")
    return geo_to_bo(geo)


def _candidate_key(bo: dict) -> tuple:
    """Hash a candidate in normalized coordinates for exact-repeat avoidance."""
    values = []
    for name in BO_PARAM_NAMES:
        width = BO_UPPER[name] - BO_LOWER[name]
        values.append(round((float(bo[name]) - BO_LOWER[name]) / width, 12))
    return tuple(values)


def evaluate(bo_params: dict) -> tuple[float, float] | tuple[None, None]:
    """Run one CFD case and return (pressure ratio, segregation intensity)."""
    try:
        geo_params = bo_to_geo(bo_params)
    except ValueError:
        print(f"[bo] SKIP infeasible params: {bo_params}", file=sys.stderr)
        return None, None
    if not _is_geo_feasible(geo_params):
        print(f"[bo] SKIP infeasible params: {bo_params}", file=sys.stderr)
        return None, None

    realized_metrics = _realized_geometry_metrics(
        geo_params["w_s"],
        geo_params["delta"],
        geo_params["k"],
        geo_params["L_s"],
        geo_params["L_m"],
    )
    sid = next_sample_id()
    sample_dir = RESULTS_DIR / sid
    sample_dir.mkdir(parents=True, exist_ok=True)

    geo_yaml = dict(TEMPLATE_GEO)
    geo_yaml.update(geo_params)
    with open(sample_dir / CAD_CONFIG_NAME, "w") as fh:
        yaml.dump(geo_yaml, fh, default_flow_style=False, sort_keys=False)

    print(f"\n[bo] sample {sid}: bo={bo_params}")
    print(f"[bo] sample {sid}: geo={geo_params}")
    print(f"[bo] sample {sid}: realized={_format_realized_geometry_metrics(realized_metrics)}")

    ret = padm_runner.run_design(sample_dir, NP, PROFILE_DIR)
    if ret.returncode != 0:
        kind, why = padm_runner.classify_failure(sample_dir)
        if kind == "launch":
            # NOT a design result.  Recording it as one would teach the GP that a
            # perfectly good region of the design space is infeasible, and nothing
            # in the campaign output would ever say otherwise.  Stop instead.
            raise padm_runner.LaunchFailure(
                f"sample {sid} could not be executed: {why}.  "
                "This is an environment failure, not a design failure, so it has "
                "NOT been recorded as one. Fix the environment and re-run; the "
                "sample directory is left in place and will be retried."
            )
        print(f"[bo] WARNING: sample {sid} failed ({why}) - excluded from GP",
              file=sys.stderr)
        # Raises if this is the Nth identical failure in a row.
        FAILURE_STREAK.record_failure(sample_dir)
        _write_failure_objectives(sample_dir, bo_params, geo_params)
        return None, None

    FAILURE_STREAK.record_success()

    obj_csv = sample_dir / "objectives.csv"
    if not obj_csv.exists():
        print(f"[bo] WARNING: objectives.csv missing for {sid} - excluded from GP",
              file=sys.stderr)
        _write_failure_objectives(sample_dir, bo_params, geo_params)
        return None, None

    _annotate_objectives_csv(sample_dir, bo_params, geo_params)

    with open(obj_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            row = _normalize_objectives_row(row)
            pressure_ratio = float(row[OBJ_PRESSURE_RATIO])
            pressure_pa = float(row[OBJ_PRESSURE_PA])
            j_mix = float(row[OBJ_MIX])
            print(
                f"[bo] {sid}: dp={pressure_pa:.4g} Pa, "
                f"dp/dp_straight={pressure_ratio:.4g}, "
                f"flux-weighted MI={1.0 - math.sqrt(max(0.0, j_mix)):.4f}"
            )
            return pressure_ratio, j_mix

    return None, None


# ---------------------------------------------------------------------------
# BO helpers
# ---------------------------------------------------------------------------

def fit_model(X: torch.Tensor, Y: torch.Tensor, warm_start: dict | None = None):
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.input import Normalize
    from botorch.models.transforms.outcome import Standardize
    from gpytorch.kernels import MaternKernel, ScaleKernel
    from gpytorch.mlls import ExactMarginalLogLikelihood

    output_batch = torch.Size([Y.shape[-1]])
    covar_module = ScaleKernel(
        MaternKernel(
            nu=2.5,
            ard_num_dims=X.shape[-1],
            batch_shape=output_batch,
        ),
        batch_shape=output_batch,
    )
    model = SingleTaskGP(
        X,
        -Y,
        input_transform=Normalize(d=X.shape[-1]),
        outcome_transform=Standardize(m=Y.shape[-1]),
        covar_module=covar_module,
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
    try:
        from botorch.acquisition.multi_objective import (
            qLogNoisyExpectedHypervolumeImprovement as qNEHVIClass,
        )
        acqf_name = "qLogNoisyExpectedHypervolumeImprovement"
    except ImportError:
        from botorch.acquisition.multi_objective import (
            qNoisyExpectedHypervolumeImprovement as qNEHVIClass,
        )
        acqf_name = "qNoisyExpectedHypervolumeImprovement"
    from botorch.optim import optimize_acqf

    baseline = _straight_baseline()
    ref_point_min = torch.tensor(
        [REF_PRESSURE_PA / baseline["pressure_drop_Pa"], REF_MIX]
    )
    ref_point = -ref_point_min

    print(f"[bo] acquisition: {acqf_name}")

    acqf = qNEHVIClass(
        model=model, ref_point=ref_point, X_baseline=X, prune_baseline=True
    )

    opt_kwargs = {
        "bounds": BOUNDS,
        "q": Q_BATCH,
        "num_restarts": ACQ_NUM_RESTARTS,
        "raw_samples": ACQ_RAW_SAMPLES,
    }
    if INEQ_CONSTRAINTS:
        opt_kwargs["inequality_constraints"] = INEQ_CONSTRAINTS
    cand, _ = optimize_acqf(acqf, **opt_kwargs)
    return cand.squeeze(0)


def fallback_novel_candidate(attempted_keys: set[tuple]) -> torch.Tensor:
    """Return a deterministic Sobol fallback not already attempted.

    Failed evaluations are deliberately absent from GP training. Without this
    guard, a deterministic acquisition optimizer can immediately select the
    same failed point again because the fitted model has not changed.
    """
    engine = torch.quasirandom.SobolEngine(
        dimension=len(BO_PARAM_NAMES), scramble=True, seed=SOBOL_SEED + 104729
    )
    batch_size = max(256, len(attempted_keys) + 1)
    for _ in range(16):
        unit = engine.draw(batch_size).to(dtype=BOUNDS.dtype)
        candidates = BOUNDS[0] + unit * (BOUNDS[1] - BOUNDS[0])
        for candidate in candidates:
            bo = {name: float(candidate[j]) for j, name in enumerate(BO_PARAM_NAMES)}
            if _candidate_key(bo) not in attempted_keys:
                print("[bo] acquisition repeated an attempted point; using novel Sobol fallback")
                return candidate
    raise RuntimeError("could not find a novel fallback candidate")


# ---------------------------------------------------------------------------
# Resume helpers
# ---------------------------------------------------------------------------

def collect_existing() -> tuple:
    """Load all completed samples from RESULTS_DIR.

    New runs store bo_<name> and geo_<name> columns. Physical geometry columns
    are always projected into the current latent transform, including for
    archived rows created by an earlier transform. Samples outside the current
    mesh-safe admissible box are skipped.
    """
    if not RESULTS_DIR.exists():
        return None, None, 0

    xs, ys = [], []
    seen_success: set[tuple] = set()
    skipped_legacy = 0
    for d in sorted(RESULTS_DIR.iterdir()):
        if not (d.is_dir() and d.name.isdigit()):
            continue
        obj_csv = d / "objectives.csv"
        if not obj_csv.exists():
            continue

        with open(obj_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                row = _normalize_objectives_row(row)
                try:
                    bo = _bo_from_row(row)
                    pressure_ratio = float(row[OBJ_PRESSURE_RATIO])
                    j_mix = float(row[OBJ_MIX])
                except (KeyError, ValueError):
                    continue
                except Exception:
                    skipped_legacy += 1
                    continue

                if not is_feasible(bo):
                    skipped_legacy += 1
                    continue
                key = _candidate_key(bo)
                if key in seen_success:
                    print(
                        f"[bo] WARNING: ignoring duplicate successful design in sample {d.name}",
                        file=sys.stderr,
                    )
                    continue
                seen_success.add(key)
                xs.append([float(bo[k]) for k in BO_PARAM_NAMES])
                ys.append([pressure_ratio, j_mix])

    if skipped_legacy > 0:
        print(
            f"[bo] WARNING: skipped {skipped_legacy} legacy/incompatible sample(s) "
            "outside the current mesh-safe BO parameterisation",
            file=sys.stderr,
        )

    if not xs:
        return None, None, 0
    return torch.tensor(xs), torch.tensor(ys), len(xs)


def collect_attempted_keys() -> set[tuple]:
    """Load successful and failed candidate locations from the campaign."""
    keys: set[tuple] = set()
    if not RESULTS_DIR.exists():
        return keys

    for sample_dir in sorted(RESULTS_DIR.iterdir()):
        if not (sample_dir.is_dir() and sample_dir.name.isdigit()):
            continue
        obj_csv = sample_dir / "objectives.csv"
        if not obj_csv.exists():
            continue
        with open(obj_csv, newline="") as handle:
            for row in csv.DictReader(handle):
                try:
                    bo = _bo_from_row(_normalize_objectives_row(row))
                except Exception:
                    continue
                if is_feasible(bo):
                    keys.add(_candidate_key(bo))
    return keys


def _attempt_counts() -> tuple[int, int]:
    """Return (attempted, failed) from persisted objective rows."""
    attempted = 0
    failed = 0
    if not RESULTS_DIR.exists():
        return attempted, failed
    for sample_dir in sorted(RESULTS_DIR.iterdir()):
        if not (sample_dir.is_dir() and sample_dir.name.isdigit()):
            continue
        path = sample_dir / "objectives.csv"
        if not path.exists():
            continue
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            continue
        attempted += 1
        row = rows[-1]
        if str(row.get("failed", "false")).lower() == "true" or row.get(OBJ_MIX, "") in ("", None):
            failed += 1
    return attempted, failed


def write_screening_gate(Y: torch.Tensor | None) -> dict:
    """Evaluate and persist the predeclared topology feasibility gate."""
    attempted, failed = _attempt_counts()
    successful = 0 if Y is None else int(Y.shape[0])
    best_mixing_index = (
        None
        if Y is None or successful == 0
        else float((1.0 - Y[:, 1].clamp(min=0.0).sqrt()).max())
    )
    failed_fraction = failed / attempted if attempted else 0.0
    enough_successes = successful >= SCREENING_N_INIT
    mixing_passed = (
        best_mixing_index is not None
        and best_mixing_index >= SCREENING_MIN_MIXING_INDEX
    )
    failure_rate_passed = failed_fraction <= SCREENING_MAX_FAILED_FRACTION
    passed = enough_successes and mixing_passed and failure_rate_passed
    report = {
        "schema_version": 1,
        "campaign": CFG["campaign"],
        "passed": passed,
        "successful_designs": successful,
        "required_successful_designs": SCREENING_N_INIT,
        "attempted_designs": attempted,
        "failed_designs": failed,
        "failed_fraction": failed_fraction,
        "maximum_failed_fraction": SCREENING_MAX_FAILED_FRACTION,
        "best_flux_weighted_mixing_index": best_mixing_index,
        "minimum_mixing_index": SCREENING_MIN_MIXING_INDEX,
        "decision": (
            "continue_to_full_sequential_bo"
            if passed
            else "adapt_topology_before_full_sequential_bo"
        ),
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SCREENING_GATE_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        f"[bo] screening gate: {'PASS' if passed else 'NO-GO'}; "
        f"best MI={best_mixing_index}, failed fraction={failed_fraction:.1%}"
    )
    print(f"[bo] screening report -> {SCREENING_GATE_PATH}")
    return report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-new-evaluations",
        type=int,
        default=None,
        help="cap new CFD evaluations in this invocation while preserving total targets",
    )
    parser.add_argument(
        "--stage",
        choices=("screening", "optimization"),
        default="screening",
        help="default stops at the 12-point feasibility gate; optimization requests the full campaign",
    )
    parser.add_argument(
        "--override-screening-gate",
        action="store_true",
        help="continue despite a recorded no-go decision (the override is printed and remains explicit)",
    )
    parser.add_argument(
        "--profile",
        default=None,
        help="Snakemake workflow profile directory selecting the execution backend "
             "(default $PADM_SNAKEMAKE_PROFILE, else profiles/local)",
    )
    return parser.parse_args()


def main() -> None:
    global PROFILE_DIR
    args = _parse_args()
    if args.max_new_evaluations is not None and args.max_new_evaluations < 1:
        raise ValueError("--max-new-evaluations must be positive")

    # Resolve and probe the backend BEFORE the first design is consumed: a
    # missing cfMesh or an unbuilt function-object library would otherwise
    # surface as a run of designs that all "failed".
    PROFILE_DIR = padm_runner.resolve_profile(args.profile)
    padm_runner.preflight(PROFILE_DIR)
    print(f"[bo] backend: profile {PROFILE_DIR.name}, np={NP}")

    # A corrected BO result is meaningful only relative to a validated
    # straight-channel baseline at the same operating point.
    _straight_baseline()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    remaining_this_run = args.max_new_evaluations
    target_init = SCREENING_N_INIT if args.stage == "screening" else N_INIT

    def budget_available() -> bool:
        return remaining_this_run is None or remaining_this_run > 0

    def consume_budget() -> None:
        nonlocal remaining_this_run
        if remaining_this_run is not None:
            remaining_this_run -= 1

    X_obs, Y_obs, n_existing = collect_existing()
    attempted_keys = collect_attempted_keys()
    if args.stage == "optimization" and n_existing < SCREENING_N_INIT:
        raise RuntimeError(
            f"complete the corrected {SCREENING_N_INIT}-design screening stage "
            "before requesting full optimization"
        )
    n_init_done = min(n_existing, target_init)
    if n_existing > 0:
        n_bo_existing = max(0, n_existing - N_INIT)
        print(f"[bo] resuming: {n_existing} completed sample(s) found "
              f"({n_init_done} init, {n_bo_existing} BO)")
        aggregate_all()
        plot_pareto_front(X_obs, Y_obs, min(n_existing, N_INIT),
                          f"PADM - resumed  [{n_existing} sample(s)]")

    if args.stage == "optimization" and n_existing >= SCREENING_N_INIT:
        gate = write_screening_gate(Y_obs)
        if not gate["passed"] and not args.override_screening_gate:
            raise RuntimeError(
                "corrected screening gate is NO-GO. Adapt the topology before the "
                "full BO campaign, or use --override-screening-gate to record an "
                "explicit methodological override."
            )
        if not gate["passed"]:
            print("[bo] WARNING: screening no-go explicitly overridden", file=sys.stderr)

    n_init_needed = max(0, target_init - n_init_done)
    if n_init_needed > 0:
        from botorch.utils.sampling import draw_sobol_samples

        print(f"[bo] === Sobol initialisation: {n_init_needed} remaining of {target_init} ===")
        # Use attempted (successful + failed) locations as the Sobol sequence
        # offset. A failed point is not training data, but it is never retried
        # indefinitely on the next bounded invocation.
        n_init_attempted = len(attempted_keys)
        n_proposals = target_init - n_init_done
        sobol_X = draw_sobol_samples(
            bounds=BOUNDS,
            n=n_init_attempted + n_proposals,
            q=1,
            seed=SOBOL_SEED,
        ).squeeze(1)[n_init_attempted:]

        for x in sobol_X:
            if not budget_available():
                break
            bo_params = {k: float(x[j]) for j, k in enumerate(BO_PARAM_NAMES)}
            j_dp, j_mix = evaluate(bo_params)
            attempted_keys.add(_candidate_key(bo_params))
            consume_budget()
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
                f"PADM - init {n_init_done}/{target_init}  [{X_obs.shape[0]} sample(s)]",
            )

    if n_init_done < target_init:
        print(
            f"[bo] paused after bounded sequential run: "
            f"{n_init_done}/{target_init} initial designs complete"
        )
        return

    if args.stage == "screening":
        write_screening_gate(Y_obs)
        print(
            "[bo] corrected feasibility screen complete; full optimization was "
            "not launched automatically"
        )
        return

    if X_obs is None or X_obs.shape[0] < 2:
        print("[bo] not enough successful samples to fit a GP")
        return

    n_bo_done = max(0, X_obs.shape[0] - N_INIT)
    n_bo_to_run = max(0, N_BO - n_bo_done)
    if remaining_this_run is not None:
        n_bo_to_run = min(n_bo_to_run, remaining_this_run)
    print(
        f"\n[bo] === Sequential BO: launching {n_bo_to_run} new iteration(s) "
        f"toward total target {N_BO} (existing BO samples: {n_bo_done}, q=1) ==="
    )

    warm_start = load_model_hyperparams()

    for i in range(n_bo_to_run):
        bo_iter = n_bo_done + i + 1
        print(
            f"\n[bo] --- BO iteration {bo_iter} "
            f"(launch {i + 1}/{n_bo_to_run}; observed so far: {X_obs.shape[0]}) ---"
        )

        model = fit_model(X_obs, Y_obs, warm_start=warm_start)
        save_model(model)
        warm_start = {k: v.detach().clone() for k, v in model.named_parameters()}
        x_next = next_candidate(model, X_obs, Y_obs)
        bo_params = {k: float(x_next[j]) for j, k in enumerate(BO_PARAM_NAMES)}
        if _candidate_key(bo_params) in attempted_keys:
            x_next = fallback_novel_candidate(attempted_keys)
            bo_params = {k: float(x_next[j]) for j, k in enumerate(BO_PARAM_NAMES)}

        j_dp, j_mix = evaluate(bo_params)
        attempted_keys.add(_candidate_key(bo_params))
        consume_budget()
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
            f"PADM - BO iter {bo_iter}  [{X_obs.shape[0]} sample(s)]",
        )

    from botorch.utils.multi_objective.pareto import is_non_dominated

    pareto = is_non_dominated(-Y_obs)
    print("\n[bo] Pareto-optimal designs:")
    print(
        f"  {'a_weak':>7}  {'a_strong':>8}  {'t_s':>6}  {'t_m':>6}  {'L_s':>6}  {'L_m':>6}"
        f"  {'dp/dp0':>12}  {'flux MI':>9}"
    )
    for x, y in zip(X_obs[pareto].tolist(), Y_obs[pareto].tolist()):
        bo_params = {k: float(x[j]) for j, k in enumerate(BO_PARAM_NAMES)}
        geo = bo_to_geo(bo_params)
        print(
            f"  {geo['a_weak']:7.3f}  {geo['a_strong']:8.3f}"
            f"  {geo['t_s']:6.3f}  {geo['t_m']:6.3f}"
            f"  {geo['L_s']:6.3f}  {geo['L_m']:6.3f}"
            f"  {y[0]:12.4g}  {1.0 - math.sqrt(max(0.0, y[1])):9.4f}"
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
