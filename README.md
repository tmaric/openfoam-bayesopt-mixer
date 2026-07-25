# openfoam-bayesopt-mixer

Bayesian optimization of parameterized passive micromixers with a portable
CAD-to-CFD workflow built on CadQuery, cfMesh, OpenFOAM, Snakemake, and
BoTorch.

## Active study: Planar Alternating-Deflector Micromixer

The active study is `PlanarAlternatingDeflectorMixer/` (PADM). The former name
`SplitAndRecombineMixer` was retired because the implemented device is a
strictly planar channel with a centre baffle and alternating top/bottom wall
deflectors. It stretches and diffuses the inlet scalar interface, but it has no
three-dimensional branch exchange and the computed fields do not show the
layer multiplication claimed by a true split-and-recombine mixer.

The study minimizes two objectives:

- kinematic pressure drop, `J_dp = <p>_in - <p>_out`, in `m²/s²`;
- outlet intensity of segregation, `J_mix = I_s`, where zero is perfectly
  mixed and one retains the inlet variance.

For water, physical pressure drop in pascals is `rho * J_dp`. The original CSV
column suffix `_Pa` was dimensionally incorrect; new runs use `_m2_s2`, while
the Python tools continue to read legacy results.

## Repository layout

```text
openfoam-bayesopt-mixer/
├── Allwmake                         build the local function-object library
├── Allwclean
├── src/functionObjects/
│   ├── pressureDrop/                area-averaged kinematic pressure drop
│   └── patchMixingQuality/          outlet scalar statistics
└── PlanarAlternatingDeflectorMixer/
    ├── FlowCase/                    CadQuery + mesh + simpleFoam template
    ├── ScalarTransportCase/         scalarTransportFoam template
    ├── Snakefile                    isolated CAD-to-objectives workflow
    ├── bayes_optimize_sequential.py sequential multi-objective BO
    ├── bayes_optimize_sequential.yaml
    ├── QUICKSTART.md
    └── docs/                        Reveal.js study deck
```

`ChannelTwoSquareObstacles/` is a separate experimental benchmark and is not
part of the PADM study.

## Geometry and physics

The default channel is `H = 1 mm`, `L = 24 mm`, with five repeated cells. Each
cell contains:

1. a centre baffle of thickness `t_s` over the split-length segment;
2. two cosine wall deflectors over an interaction length `L_c`;
3. a thinner centre baffle of thickness `t_m` over the merge-length segment.

The additional bias alternates between the top and bottom deflector. Its
realized cell value is

```text
delta_i = delta + k * xhat_i
```

where `xhat_i` is the normalized streamwise midpoint of cell `i`.

The flow is steady and laminar (`simpleFoam`, `Re = 10`). Scalar transport uses
`DT = 1e-9 m²/s` and bounded first-order upwind convection. The latter is
robust but numerically diffusive, so the observed mixing should not be treated
as mesh-independent physical validation without a discretization study.

## Prerequisites

The case was developed with OpenFOAM v2506 and Python 3.10+. The Python
environment must provide CadQuery, Snakemake, PyYAML, PyTorch, BoTorch,
GPyTorch, NumPy, Matplotlib, and Pillow. ParaView batch mode is optional and is
used only for field images.

Source the `etc/bashrc` belonging to the OpenFOAM installation you want to use;
the repository does not assume where OpenFOAM is installed.

```bash
source /path/to/OpenFOAM/etc/bashrc
```

## Build and run

From the repository root:

```bash
./Allwmake
cd PlanarAlternatingDeflectorMixer
```

The build is written under the repository's `platforms/$WM_OPTIONS/lib/`, and
the Snakemake workflow loads the library from that relative location.

Run one design:

```bash
snakemake --cores 4
```

Run sequential multi-objective Bayesian optimization:

```bash
python3 bayes_optimize_sequential.py
```

The default campaign uses eight Sobol samples followed by twenty sequential
qLogNEHVI/qNEHVI suggestions. A rerun resumes completed samples and adds
another configured BO batch.

Clean generated data:

```bash
snakemake --cores 1 clean
```

## Documentation deck

The Reveal.js deck documents the actual geometry, parameter transforms,
objectives, workflow, results, and limitations:

```bash
cd PlanarAlternatingDeflectorMixer/docs
python3 -m http.server 8000
```

Then open `http://localhost:8000/`.

## Main outputs

```text
results/<sample_id>/
├── FlowCase/pressureDrop.csv
├── ScalarTransportCase/mixing.csv
├── objectives.csv
└── visualizations/<sample_id>_T.png   optional
```

`results/all_samples.csv` aggregates all completed samples and
`results/pareto_front.png` summarizes the current non-dominated set.
