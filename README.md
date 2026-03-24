# openfoam-bayesopt-mixer

Bayesian optimization of passive laminar micromixers using a fully automated
CAD-to-CFD pipeline built on CADquery, cfMesh, OpenFOAM, Snakemake, and
BOTorch.

## Overview

The goal is to find Pareto-optimal mixer geometries that minimize two competing
objectives simultaneously:

- **Mixing defect** `J_mix` — variance of the passive scalar concentration at
  the outlet (lower = better mixed).
- **Pressure drop** `J_dp` — area-averaged pressure difference between inlet
  and outlet (lower = less pumping power).

Each design candidate is described by a small parameter vector `theta` (stored
in a YAML file). For every `theta` the pipeline

1. generates a 2-D CAD geometry (CADquery),
2. meshes the fluid domain (cfMesh `cartesian2DMesh`),
3. solves the steady laminar flow (OpenFOAM `simpleFoam`),
4. solves passive scalar transport (OpenFOAM `scalarTransportFoam`),
5. extracts `J_dp` and `J_mix` from the computed fields,
6. returns the objectives to the Bayesian optimizer (BOTorch).

## Repository layout

```
openfoam-bayesopt-mixer/
├── Allwmake                        # build custom OpenFOAM function objects
├── Allwclean
├── src/
│   └── functionObjects/
│       ├── pressureDrop/           # computes J_dp, writes pressureDrop.csv
│       └── patchMixingQuality/     # computes J_mix, writes mixing.csv
└── SplitAndRecombineMixer/                   # SAR lamination ladder mixer case
    ├── Allrun                                # run Hydro then Mixing sequentially
    ├── Allclean                              # clean both sub-cases
    ├── Snakefile                             # Snakemake workflow
    ├── QUICKSTART.md                         # quickstart for the Snakemake workflow
    ├── postprocessing_agglomeration.py       # merges YAML + CSVs into objectives.csv
    ├── SplitAndRecombineHydro/               # flow sub-case (template)
    │   ├── Allrun
    │   ├── Allclean
    │   ├── sar_mixer_cad.py                  # CADquery geometry script
    │   ├── sar_mixer_cad.yaml                # geometry parameters (edit this)
    │   └── system/, constant/, 0/
    └── SplitAndRecombineMixing/              # scalar transport sub-case (template)
        ├── Allrun                            # copies mesh+fields from Hydro, runs transport
        ├── Allclean
        └── system/, constant/, 0/
```

## The Split-and-Recombine (SAR) Lamination Ladder Mixer

The mixer is a 2-D rectangular channel of height `H` and total length
`L = 2 L0 + N L_cell` containing `N` identical unit cells.  Each unit cell
performs three actions:

```
 inlet ──► [ split ] ──► [ shuffle ] ──► [ recombine ] ──► next cell ──► outlet
```

1. **Split** — a thin horizontal splitter wall (`thickness t_s`) divides the
   channel into two subchannels of width `w_s` each.
2. **Shuffle** — cosine-shaped deflectors on the top and bottom walls
   (`height h_d = H/2 - w_s`) with a vertical offset `delta` redirect the
   sub-streams so they swap positions.
3. **Recombine** — a second thin splitter (`thickness t_m`) guides the
   sub-streams back into a single channel, now with laminated interfaces.

After `N` cells the number of laminae doubles with each cell, reducing striation
thickness and accelerating diffusion.

### Geometry parameters (`sar_mixer_cad.yaml`)

All lengths are in normalised units where `H = 1`.  The `scale` factor converts
to SI metres.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scale`   | `1e-3`  | 1 normalised unit = 1 mm |
| `H`       | `1.0`   | channel height |
| `L0`      | `2.0`   | inlet / outlet buffer length |
| `N`       | `5`     | number of unit cells |
| `L_cell`  | `4.0`   | unit-cell length |
| `L_s`     | `1.4`   | split section length |
| `L_m`     | `1.0`   | merge section length |
| `w_s`     | `0.38`  | subchannel half-gap after split (fraction of H) |
| `t_s`     | `0.10`  | splitter thickness — split section (fraction of H) |
| `t_m`     | `0.05`  | splitter thickness — merge section (fraction of H) |
| `delta`   | `0.08`  | top deflector vertical bias (fraction of H) |

Edit `SplitAndRecombineMixer/SplitAndRecombineHydro/sar_mixer_cad.yaml` to
change the geometry before running.

## Prerequisites

| Tool | Version tested |
|------|---------------|
| OpenFOAM | v2506 |
| CADquery | 2.x |
| cfMesh (cartesian2DMesh) | bundled with OpenFOAM-v2506 |
| Snakemake | 8.x |
| Python | 3.10+ |

Source the OpenFOAM environment before building or running:

```bash
source $HOME/OpenFOAM/OpenFOAM-v2506/etc/bashrc
```

## Build

Compile the custom function objects (`pressureDrop`, `patchMixingQuality`):

```bash
./Allwmake
```

This places `libbayesoptMixerFunctionObjects.so` under
`platforms/$WM_OPTIONS/lib/`.

## Running

### Option A — Allrun / Allclean scripts

Run both sub-cases sequentially from the mixer directory:

```bash
cd SplitAndRecombineMixer
./Allrun
```

`Allrun` calls `SplitAndRecombineHydro/Allrun` first, then
`SplitAndRecombineMixing/Allrun`.  The Mixing `Allrun` automatically copies
the mesh and latest flow fields from the Hydro case before running scalar
transport.

Clean all generated files:

```bash
./Allclean
```

Each sub-case also has its own `Allrun` and `Allclean` that can be run
independently.  Running `SplitAndRecombineMixing/Allrun` directly requires the
Hydro case to have been run first (the script checks for the Hydro mesh and
exits with a clear error if it is missing).

### Option B — Snakemake workflow

The Snakefile in `SplitAndRecombineMixer/` orchestrates the full pipeline as a
directed acyclic graph of rules, staging results under `results/` to keep
template cases untouched.

```bash
cd SplitAndRecombineMixer
snakemake -j 4          # run with 4 cores (used by mpirun for both solvers)
snakemake -j 1 clean    # remove results/
```

See `SplitAndRecombineMixer/QUICKSTART.md` for a concise reference.

The workflow proceeds as follows:

```
stage_template_cases
        │
        ├─► hydro_geometry      (CADquery: generate STL)      → log.sar_mixer_cad
        │        │
        │   hydro_mesh          (cartesian2DMesh)              → log.cartesian2DMesh
        │        │
        │   hydro_cell_volumes  (postProcess writeCellVolumes) → log.postProcess
        │        │
        │   hydro_decompose     (decomposePar, N = -j N)       → log.decomposePar
        │        │
        │   hydro_simpleFoam    (mpirun -np N simpleFoam)      → log.simpleFoam
        │        │  also writes: pressureDrop.csv
        │        │
        │   hydro_reconstruct   (reconstructPar -latestTime)   → log.reconstructPar
        │        │
        │   copy_hydro_to_mixing
        │        │
        │   mixing_set_expr_fields  (setExprFields)            → log.setExprFields
        │        │
        │   mixing_decompose    (decomposePar, N = -j N)       → log.decomposePar
        │        │
        │   mixing_scalar_transport (mpirun -np N scalarTransportFoam) → log.scalarTransportFoam
        │        │  also writes: mixing.csv
        │        │
        │   mixing_reconstruct  (reconstructPar -latestTime)   → log.reconstructPar
        │        │
        ├─► agglomerate         (postprocessing_agglomeration.py)
        │        │  writes: objectives.csv
        │
        └─► create_foam_files   (.foam files for ParaView)
```

Every rule captures its application output via `tee log.<appname>` inside the
results directory.  `pressureDrop.csv` and `mixing.csv` are declared Snakemake
outputs of the respective solver rules — if a function object fails to write its
CSV the rule fails immediately rather than silently producing empty results.

### Option C — Bayesian optimisation loop

Run the full multi-objective optimisation from `SplitAndRecombineMixer/`:

```bash
python bayes_optimize.py --n-init 8 --n-bo 20 --cores 4
```

| Flag | Default | Description |
|------|---------|-------------|
| `--n-init` | `8` | Sobol initial samples before GP is fit |
| `--n-bo` | `20` | BO iterations after initialisation |
| `--cores` | `4` | CPU cores passed to each Snakemake call |
| `--results-dir` | `results/` | Root directory for per-sample results |

The loop assigns sequential zero-padded IDs (`00000`, `00001`, …).  For each
candidate it:

1. writes `results/{id}/sar_mixer_cad.yaml` (geometry parameters),
2. calls Snakemake with `--snakefile`, `--directory results/{id}`, and
   `--config results_dir=results/{id}` so each run is fully isolated
   (`.snakemake/` metadata and all CFD outputs live under `results/{id}/`),
3. reads `results/{id}/objectives.csv` and updates the BOTorch GP model.

After all iterations `results/all_objectives.csv` aggregates every sample, and
the Pareto-optimal designs are printed to stdout.

The loop resumes automatically: if `results/` already contains completed
samples they are loaded and counted toward the Sobol initialisation budget
before any new Snakemake calls are made. Once the Sobol phase is complete,
each fresh invocation launches another batch of BO iterations and appends new
sample directories under `results/`.

#### Manual batch invocation

Each Snakemake run can also be triggered by hand for a specific parameter set:

```bash
# write sar_mixer_cad.yaml into the sample directory first, then:
snakemake \
  --snakefile SplitAndRecombineMixer/Snakefile \
  --directory results/00042 \
  --cores 4 \
  --config results_dir=results/00042
```

## Outputs

| File | Location | Contents |
|------|----------|----------|
| `pressureDrop.csv` | Hydro case dir | time, inlet/outlet average pressure, `J_dp` |
| `mixing.csv` | Mixing case dir | time, all mixing quality metrics |
| `objectives.csv` | `results/` (or `samples/{id}/`) | single-row: geometry params + all objectives |
| `log.sar_mixer_cad` | Hydro case dir | CADquery geometry script output |
| `log.cartesian2DMesh` | Hydro case dir | cfMesh output |
| `log.postProcess` | Hydro case dir | writeCellVolumes output |
| `log.decomposePar` | Hydro / Mixing case dir | domain decomposition output |
| `log.simpleFoam` | Hydro case dir | full simpleFoam stdout |
| `log.setExprFields` | Mixing case dir | scalar field initialisation output |
| `log.scalarTransportFoam` | Mixing case dir | full scalarTransportFoam stdout |
| `log.reconstructPar` | Hydro / Mixing case dir | reconstruction output |
| `*.foam` | Hydro / Mixing case dir | ParaView session files |
