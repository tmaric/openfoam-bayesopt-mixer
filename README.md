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
- flux-weighted outlet intensity of segregation, `J_mix = I_s,flux`, where
  zero is perfectly mixed and one retains the inlet variance.

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
    ├── research_sequence.py         gated, one-step research driver
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

The verified BO coordinates sample the weak-wall amplitude directly and map a
ratio into an admissible strong-wall amplitude:

```text
w_s = 0.5 - a_weak
delta = a_strong - a_weak
k = 0
```

The strong wall alternates between the top and bottom across the five cells.
The constant-amplitude choice removes the old downstream slope until evidence
supports adding it back.

The flow is steady and laminar (`simpleFoam`, `Re = 10`). Scalar transport uses
`DT = 1e-9 m²/s`, bounded second-order `limitedLinear 1` convection, and a
PBiCGStab/DILU scalar solve. Outlet objectives must satisfy a final-window
stability test. Results still require a mesh/scheme study before physical
publication claims.

## Prerequisites

The case was developed with OpenFOAM v2506 and Python 3.10+. The Python
environment must provide CadQuery, Snakemake, foamlib, Python VTK, PyYAML,
PyTorch, BoTorch, GPyTorch, NumPy, Matplotlib, and Pillow. Per-sample field
images are rendered directly in Python without ParaView, OpenGL, or a display
server.

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

Run one explicitly configured design with no more than two MPI ranks:

```bash
snakemake --cores 2 --config results_dir=results/manual_00
```

Advance the sequential multi-objective campaign by one evaluation:

```bash
python3 research_sequence.py status
python3 research_sequence.py next --max-new-evaluations 1
```

The corrected campaign first requires three matched baselines and a 12-point
feasibility screen. After a pass it targets 32 total Sobol designs followed by 80
strictly sequential (`q = 1`) qLogNEHVI suggestions. A rerun resumes toward
those totals; it does not add a fresh 80-point batch. OpenFOAM uses two MPI
ranks and Torch uses one thread. Failed cases are retained in the audit trail,
excluded from GP fitting, and never converted into fictitious penalties.

The corrected screen completed on 2026-07-26 with 12/12 numerically successful
designs but a scientific NO-GO: the best mixing index was 0.1698 versus the
required 0.60, at 34.060 Pa. Full BO is intentionally blocked until a new
topology with an interface-multiplication mechanism passes the same screen.

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
results/corrected_boundary_v3/<sample_id>/
├── FlowCase/pressureDrop.csv
├── ScalarTransportCase/mixing.csv
├── objectives.csv
└── visualizations/<sample_id>_T.png   optional
```

The campaign directory contains `all_samples.csv`, `pareto_front.png`, and the
GP checkpoint. Generated results are ignored by Git. See
`PlanarAlternatingDeflectorMixer/docs/RESEARCH_PLAN.md` for the staged
verification and publication gates.
