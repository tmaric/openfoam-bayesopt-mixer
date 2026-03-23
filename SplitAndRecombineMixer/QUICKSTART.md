# Quickstart — SplitAndRecombineMixer

## What is being optimised

The **Split-and-Recombine (SAR) laminar mixer** is a passive microfluidic
device.  A straight rectangular channel is interrupted by *N* identical unit
cells (default N = 5).  Each cell splits the flow into two sub-channels,
deflects each stream with a cosine-shaped ramp, then merges them back.
Repeated splitting and recombining stretches the interface between the two
inlet streams, increasing mixing by laminar advection.

The geometry is 2-D (a thin slab extruded in z).  Physical dimensions are set
by `scale = 1e-3` m/unit, so `H = 1` normalised unit = 1 mm channel height.

### Design variables

Six geometric parameters, all normalised by the channel height H = 1:

| Variable | Range        | Description |
|----------|-------------|---|
| `w_s`    | [0.25, 0.45] | Sub-channel half-gap after split (fraction of H) |
| `t_s`    | [0.02, 0.15] | Splitter plate thickness — split section |
| `t_m`    | [0.02, 0.15] | Splitter plate thickness — merge section |
| `L_s`    | [0.80, 1.80] | Split-section length (0.20–0.45 × L_cell) |
| `L_m`    | [0.80, 1.80] | Merge-section length (0.20–0.45 × L_cell) |
| `delta`  | [0.00, 0.15] | Top-deflector vertical bias (fraction of H) |

### Geometric feasibility constraints

Three constraints must hold before a geometry is meshed and simulated:

- **C1** `w_s − 0.5·t_s − delta ≥ 0.02`
  Minimum fluid gap between the deflector peak and the nearest splitter face.
  Prevents the channel from pinching off.

- **C2** `L_s + L_m ≤ 3.60`
  Keeps the central interaction region L_c = L_cell − L_s − L_m ≥ 0.4,
  where the two cosine deflectors face each other.

- **C3** `t_s − t_m ≥ 0.011`
  The split splitter must be strictly thicker than the merge splitter so that
  the Boolean STL cut at each cell boundary removes real material.

---

## Objective function

Two objectives are **minimised simultaneously**; the algorithm finds the
Pareto-optimal trade-off front between them.

### J1 — pressure drop  `pdrop_pressure_drop_Pa`

The area-averaged pressure difference between the inlet and outlet patches,
computed by OpenFOAM's `simpleFoam` steady-state flow solver:

```
J1 = p_inlet − p_outlet  [Pa]
```

Observed range across the 28 initial samples: **0.9 – 2.2 mPa**.
Higher `w_s` (wider sub-channels) and shorter sections generally reduce J1.

### J2 — intensity of segregation  `mixing_intensity_of_segregation`

The **Danckwerts intensity of segregation** I_s measures how far from uniform
the scalar concentration field T is at the outlet cross-section.  It is
defined as:

```
I_s = σ² / σ₀²
```

where σ² is the area-averaged variance of T at the outlet face and σ₀² is
the maximum possible variance for a binary mixture with inlet mean ā = 0.5:

```
σ₀² = ā(1 − ā) = 0.5 × 0.5 = 0.25
```

The variance σ² is computed about the area-averaged outlet mean (not the
flux-weighted mean), matching the `patchMixingQuality` function object
(`meanMode fromInletRatio`, `aMean 0.5`).

| I_s value | Physical meaning |
|-----------|---|
| 1.0 | Completely unmixed — step profile fully preserved at outlet |
| 0.0 | Perfectly mixed — uniform concentration at outlet |

**Column naming.**  The function object writes two related columns to
`mixing.csv`:

| CSV column | Formula | Used by BO? |
|---|---|---|
| `intensity_of_segregation` | I_s = σ²/σ₀² | **Yes** — minimised as J2 |
| `mixing_index_intensity` | 1 − I_s | No (for reference only) |

`postprocessing_agglomeration.py` prefixes all mixing columns with
`mixing_`, so the BO reads `mixing_intensity_of_segregation` = I_s directly.
The Pareto plot uses `1 − I_s` on its y-axis only for readability (higher =
better mixed), but the optimiser minimises I_s.

#### CSV output location and format

Each evaluated sample writes a `mixing.csv` at:

```
results/<sample_id>/SplitAndRecombineMixing/mixing.csv
```

One row is appended per solver iteration.  `postprocessing_agglomeration.py`
reads the **last (converged) row** and propagates `intensity_of_segregation`
into the per-sample `objectives.csv`, which the BO loop ingests.

Verbatim header and a representative converged row (sample `00000`,
iteration 98):

```
time,scalar_field,patch,weighting,mean_mode,mean_concentration,standard_deviation,coefficient_of_variation,mixing_coefficient,flux_weighted_mean_concentration,flux_weighted_standard_deviation,flux_weighted_coefficient_of_variation,flux_weighted_mixing_coefficient,intensity_of_segregation,flux_weighted_intensity_of_segregation,relative_standard_deviation,flux_weighted_relative_standard_deviation,mixing_index_rsd,mixing_index_intensity,max_absolute_relative_deviation,mean_absolute_relative_deviation,delta_x_min,delta_x_max,delta_x_mean,delta_x_average
98,T,outlet,phi,fromInletRatio,0.507547,0.432051,0.864102,0.864102,0.476483,0.428659,0.899632,0.899632,0.746673,0.734995,0.864102,0.857318,0.135898,0.253327,1,0.829418,6.74065e-06,1.42879e-05,1.23469e-05,1.23469e-05
```

Key columns at a glance:

| Column | Value | Note |
|---|---|---|
| `mean_concentration` | 0.5075 | Close to 0.5 — correct for a 50/50 split |
| `standard_deviation` | 0.4321 | σ at outlet |
| `intensity_of_segregation` | **0.7467** | **J2 — what the BO minimises** |
| `mixing_index_intensity` | 0.2533 | = 1 − I_s, for reference |

#### Expected values and what to hope for

For a 2-D laminar SAR mixer with N = 5 unit cells, molecular diffusion is
negligible (`scalarTransportFoam` with D → 0 by default).  Mixing is purely
advective: I_s decays only as fast as the geometry stretches the interface.

From the 28 Sobol-initialisation samples:

| | I_s | Mixing quality 1 − I_s |
|---|---|---|
| Best observed | ≈ 0.61 | ≈ 0.39 |
| Typical range | 0.63 – 0.90 | 0.10 – 0.37 |
| Worst observed | ≈ 0.90 | ≈ 0.10 |

**What to hope for from the BO loop:**
- **Target: I_s < 0.5** (mixing quality > 0.5).  This means more than half
  of the inlet concentration variance has been destroyed by the geometry
  alone — a strong result for purely advective mixing with only 5 unit cells.
- The BO targets the Pareto front between low pressure drop and low I_s,
  so expect it to find geometries with I_s in the 0.50 – 0.62 range that
  also improve on J1 compared to the Sobol samples.

---

## Bayesian optimisation

Run the sequential multi-objective BO loop:

```bash
cd SplitAndRecombineMixer
python bayes_optimize_sequential.py
```

The script:
1. **Sobol initialisation** (N_INIT = 8 samples, feasibility-filtered).
2. **Sequential BO loop** (N_BO = 20 iterations):
   fit `SingleTaskGP` → optimise `qNEHVI` acquisition with constraints →
   evaluate via Snakemake → update GP.
3. **Reports** Pareto-optimal designs and saves `results/pareto_front.png`.

### Resuming and extending runs

Each run is **automatically resumed** from all previous results:

- Sample directories are numbered `00000`, `00001`, … and new runs continue
  from the highest existing index.
- The Gaussian Process model is checkpointed to `SplitAndRecombineMixer.pt`
  after every BO iteration.  On restart the saved hyperparameters are used
  as a warm start, so the GP converges faster when data already exist.
- `results/all_samples.csv` is rebuilt from all per-sample `objectives.csv`
  files on every run, so visualisations always reflect the full dataset
  regardless of how many separate runs contributed to it.

---

## Prerequisites

Build the custom OpenFOAM function objects from the repository root:

```bash
./Allwmake
```

Source the OpenFOAM environment (once per shell session):

```bash
source $HOME/OpenFOAM/OpenFOAM-v2506/etc/bashrc
```

## Run with Snakemake

```bash
cd SplitAndRecombineMixer
snakemake -j N
```

`N` is the number of CPU cores available.  It is used in two ways:

- **Snakemake** uses it to schedule independent rules concurrently.
- **OpenFOAM** (`simpleFoam`, `scalarTransportFoam`) runs in parallel with
  `mpirun -np N`, so both the flow and scalar-transport solves use all N cores.

Outputs land in `results/`:

```
results/
├── SplitAndRecombineHydro/   — mesh, flow solution, pressureDrop.csv
├── SplitAndRecombineMixing/  — scalar transport solution, mixing.csv
└── objectives.csv            — agglomerated geometry + objectives (one row)
```

## Clean

```bash
snakemake -j 1 clean
```

Removes the entire `results/` directory.
