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
└── SplitAndRecombineMixer/         # SAR lamination ladder mixer case
    ├── Allrun                      # run Hydro then Mixing sequentially
    ├── Allclean                    # clean both sub-cases
    ├── Snakefile                   # Snakemake workflow (geometry → mesh → flow → mixing)
    ├── SplitAndRecombineHydro/     # flow sub-case
    │   ├── Allrun
    │   ├── Allclean
    │   ├── sar_mixer_cad.py        # CADquery geometry script
    │   ├── sar_mixer_cad.yaml      # geometry parameters (edit this)
    │   └── system/, constant/, 0/
    └── SplitAndRecombineMixing/    # scalar transport sub-case
        ├── Allrun                  # copies mesh+fields from Hydro, runs transport
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
snakemake --cores 4          # run with 4 parallel MPI ranks
snakemake --cores 4 clean    # remove results/
```

The workflow proceeds as follows:

```
stage_template_cases
        │
        ├─► hydro_geometry      (CADquery: generate STL)
        │        │
        │   hydro_mesh          (cartesian2DMesh)
        │        │
        │   hydro_cell_volumes  (postProcess writeCellVolumes)
        │        │
        │   hydro_decompose     (decomposePar, N = --cores)
        │        │
        │   hydro_simpleFoam    (mpirun simpleFoam -parallel)
        │        │  outputs: pressureDrop.csv + log.simpleFoam
        │        │
        │   hydro_reconstruct   (reconstructPar -latestTime)
        │        │
        │   copy_hydro_to_mixing
        │        │
        │   mixing_set_expr_fields
        │        │
        │   mixing_decompose
        │        │
        │   mixing_scalar_transport  (mpirun scalarTransportFoam -parallel)
        │        │  outputs: mixing.csv + log.scalarTransportFoam
        │        │
        │   mixing_reconstruct
        │
        └─► create_foam_files   (.foam files for ParaView)
```

`pressureDrop.csv` and `mixing.csv` are declared Snakemake outputs of the
respective solver rules.  If a custom function object fails to write its CSV,
the rule fails immediately with the full solver log available in
`results/SplitAndRecombineHydro/log.simpleFoam` (or `log.scalarTransportFoam`).

## Outputs

| File | Location | Contents |
|------|----------|----------|
| `pressureDrop.csv` | Hydro case dir | time, inlet/outlet average pressure, `J_dp` |
| `mixing.csv` | Mixing case dir | time, outlet mixing defect `J_mix` |
| `log.simpleFoam` | Hydro case dir (Snakemake only) | full simpleFoam stdout |
| `log.scalarTransportFoam` | Mixing case dir (Snakemake only) | full transport solver stdout |
| `*.foam` | both case dirs | ParaView session files |
