# Claude Code handover: planar alternating-deflector micromixer

Last updated: 2026-07-28

## Purpose

Continue the numerical verification and strictly sequential Bayesian-
optimization study in this repository. The active device is the **Planar
Alternating-Deflector Micromixer (PADM)**. It must not be described as a true
split-and-recombine (SAR) mixer: the model is two-dimensional and has no 3-D
branch exchange, layer permutation, or restacking.

The immediate research objective is to finish the remaining eleven corrected
feasibility-screen points without overloading the workstation. Full BO is
blocked until that screen passes. The longer-term objective is a defensible
passive-mixer comparison and a publication-quality numerical study. The
present results are not yet a publishable claim of a global optimum.

> **Current-campaign warning:** the historical `verified_flux_sequential_v2`
> discussion retained below is provenance only. Its objectives are invalid.
> The final “Superseding corrected-boundary update” is authoritative, and the
> only permitted campaign entry point is `research_sequence.py`.

## Non-negotiable instructions

- Run OpenFOAM and the optimization from **Ubuntu WSL**, not directly from
  PowerShell.
- Keep the BO strictly sequential: `q = 1`.
- Use no more than two OpenFOAM MPI ranks and one Torch/BLAS thread.
- Prefer bounded invocations of one to four new evaluations. One evaluation at
  a time is the safest default.
- Preserve failed evaluations in the audit trail and exclude them from GP
  fitting. Do not replace failures with arbitrary objective penalties.
- Do not modify, run, stage, or commit `ChannelTwoSquareObstacles/`; it is not
  part of this study.
- Preserve unrelated user changes in the working tree. Do not reset, restore,
  overwrite, or stage them without an explicit request.
- Keep repository paths portable. Source files and workflows must use paths
  resolved relative to the repository or study directory, never machine-
  specific absolute repository paths.
- Commit messages should describe the scientific or software change and must
  not mention an assistant or its identity.

## Repository and Git state

- Repository: `openfoam-bayesopt-mixer`
- Active branch: `main`
- Current committed head: `79405cd Implement verified sequential BO campaign`
- Previous rename commit: `a1132ba Rename`
- Repository-local Git identity:
  - name: `Tomislav Maric`
  - email: `tomislav.maric@gmx.com`
- Generated campaign results are ignored by Git.

At handover time, the following pre-existing changes are unstaged and must be
treated as user-owned:

```text
docs/obsidian/01-Weekly/2026-W10-(2026-03-02-to-2026-03-06).md
docs/obsidian/01-Weekly/Milestones.md
docs/obsidian/01-Weekly/Weekly-Template.md
docs/obsidian/02-Technical-Notes/Benchmark-Problem-Specification.md
docs/obsidian/02-Technical-Notes/Technical-Note-01-Passive-Mixer-Benchmark.md
docs/obsidian/02-Technical-Notes/Technical-Note-02-2D-SAR-Lamination-Ladder-Mixer.md
docs/obsidian/03-Resources/passive_mixer_parameterization_sketch.py
ChannelTwoSquareObstacles/                       (untracked and unrelated)
```

This handover file is also a new working-tree file until the user explicitly
asks for it to be staged or committed.

## Active study layout

The active study is under `PlanarAlternatingDeflectorMixer/`:

```text
FlowCase/                              CAD, mesh, and simpleFoam template
ScalarTransportCase/                   passive-scalar template
Snakefile                              isolated CAD-to-objectives workflow
bayes_optimize_sequential.py           resumable sequential BO driver
bayes_optimize_sequential.yaml         campaign and parameter configuration
verify_scalar_convergence.py           final-window scalar stability gate
visualize_results.py                   campaign-level plots
render_latest_mixing_field.py          headless foamlib/VTK/Pillow rendering
QUICKSTART.md                          concise technical and run documentation
docs/index.html                        Reveal.js presentation
docs/RESEARCH_PLAN.md                  verification/publication roadmap
results/verified_flux_sequential_v2/   generated current campaign
```

The root `Allwmake` builds the two custom OpenFOAM function objects into the
repository-local `platforms/$WM_OPTIONS/lib/` directory:

- `pressureDrop`: area-averaged kinematic pressure loss;
- `patchMixingQuality`: outlet scalar statistics, including positive-flux-
  weighted quantities.

## Device, physics, and objectives

The geometry is a 2-D channel with a centre baffle and five repeated cells of
alternating strong/weak cosine wall deflectors. The stronger intrusion
alternates between the top and bottom wall. The nominal channel height is
`H = 1 mm`; all CAD dimensions are normalized by `H` and scaled to SI units.

The current flow/scalar setup is:

- steady laminar `simpleFoam` flow at `Re = 10`;
- passive scalar diffusivity `DT = 1e-9 m^2/s`;
- bounded second-order `limitedLinear 1` scalar convection;
- PBiCGStab/DILU scalar solver with equation relaxation `0.7`;
- explicit SIMPLE convergence, allowing up to 2000 iterations;
- scalar evolution to at most 600 iterations;
- final 50 scalar samples must have spans no greater than `1e-4` in both
  flux-weighted outlet mean and flux-weighted segregation intensity.

Both optimization objectives are minimized:

```text
J_dp  = <p>_inlet - <p>_outlet                  [m^2/s^2]
J_mix = flux-weighted intensity of segregation  [-]
```

OpenFOAM pressure is kinematic. Convert to pascals with
`DeltaP = rho * J_dp`. The literature-facing mixing index is
`1 - sqrt(J_mix)`, not the legacy `1 - J_mix` display.

## Six-dimensional BO parameterization

The checked-in campaign samples:

| BO coordinate | Bounds | Meaning/mapping |
|---|---:|---|
| `a_weak` | 0.00--0.12 | weak-wall intrusion; `w_s = 0.5 - a_weak` |
| `a_strong_ratio` | 0--1 | maps into the admissible strong-wall amplitude interval |
| `t_s` | 0.04--0.15 | split-section centre-baffle thickness |
| `t_m_ratio` | 0--1 | maps into a mesh-safe merge-baffle thickness interval |
| `L_c` | 0.40--2.40 | deflector interaction length |
| `L_s_ratio` | 0--1 | partitions remaining cell length between split and merge sections |

The strong amplitude is constrained by
`max(0.12, a_weak + 0.04) <= a_strong <= 0.35`. The realized CAD variables
include `w_s`, `t_s`, `t_m`, `L_s`, `L_m`, `delta`, `k`, `a_weak`, and
`a_strong`, where `delta = a_strong - a_weak` and `k = 0`. Interval transforms
enforce minimum mesh-cell counts across thin features and ensure
`L_s + L_c + L_m = 4` with `0.8 <= L_s,L_m <= 1.8`.

Do not reintroduce the old correlated `w_s/delta` parameterization or the
downstream amplitude slope `k` without new evidence and a documented campaign
version change.

## BO configuration

The active configuration is
`PlanarAlternatingDeflectorMixer/bayes_optimize_sequential.yaml`:

- campaign: `verified_flux_sequential_v2`;
- 32 successful scrambled-Sobol initial designs;
- Sobol seed: `20260725`;
- 80 subsequent sequential BO evaluations;
- independent-output exact GPs;
- normalized six-dimensional inputs and standardized outputs;
- ARD Matern-5/2 covariance;
- qLogNEHVI with `q = 1`;
- 32 acquisition restarts and 1024 raw samples;
- fixed minimization reference point `(0.02 m^2/s^2, 1.0)`;
- two Snakemake/OpenFOAM cores;
- one Torch thread.

The driver validates the resource limits and sequential setting at startup.
Campaign targets are totals. Re-running the driver resumes toward 32 + 80; it
does not add another 80 evaluations on each invocation.

## Current campaign state

Generated results are currently present locally under
`PlanarAlternatingDeflectorMixer/results/verified_flux_sequential_v2/`.
They are ignored by Git, so a fresh clone will not contain this state unless
the result directory is copied separately.

Current status: **3/32 successful Sobol designs, 0/80 BO designs**.

| Sample | Status | `J_dp` (`m^2/s^2`) | Flux `J_mix` | `1-sqrt(J_mix)` | Notes |
|---|---|---:|---:|---:|---|
| `00000` | failed | -- | -- | -- | early workflow failure; first Sobol point |
| `00001` | failed | -- | -- | -- | repeated early failure of the same point |
| `00002` | successful | 0.00141115 | 0.780094 | 0.1168 | first verified result |
| `00003` | successful | 0.00119273 | 0.889416 | 0.0569 | verified after scalar-solver correction |
| `00004` | failed | -- | -- | -- | flow converged; scalar missed the final-window stability gate |
| `00005` | successful | 0.00108524 | 0.871213 | 0.0666 | flow converged in 425 iterations |

Sample `00004` was numerically bounded and close to steady, but the final
unweighted mixing-intensity drift over roughly 35 iterations was about
`2.2e-4`, above the configured `1e-4` gate. The exact flux-weighted verifier
span was not preserved after Snakemake removed its failed output. Do not relax
the acceptance criterion on the evidence of this one sample. Review failure
clustering at the planned checkpoints.

The driver tracks attempted parameter vectors, so a failed design is not
immediately proposed again. Failed rows have blank targets and are excluded
from GP fitting.

## Portable WSL procedure

Use an interactive Ubuntu WSL terminal and substitute the local clone and
OpenFOAM installation locations:

```bash
source /path/to/OpenFOAM-v2506/etc/bashrc

cd /path/to/openfoam-bayesopt-mixer
./Allwmake

cd PlanarAlternatingDeflectorMixer
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

python -u \
    research_sequence.py next \
    --max-new-evaluations 1
```

The Python environment currently passes imports for PyYAML, CadQuery,
Snakemake, Torch, BoTorch, GPyTorch, NumPy, Matplotlib, and Pillow.

Repeat the bounded command to advance the campaign. Increasing
`--max-new-evaluations` to 2--4 still runs candidates sequentially. Omitting
the option runs all remaining evaluations and should be done only when the
machine can remain occupied for the full campaign.

Run the command in the foreground or in a persistent WSL terminal/tmux
session. Do not launch it through a short-lived `wsl.exe` wrapper whose stdout
pipe will be closed: Snakemake rules use `tee`, and closing the host-side pipe
can produce a false `SIGPIPE`/exit-141 meshing failure. For managed unattended
execution, use a persistent WSL session or a WSL `systemd-run --user` service.

## Useful inspection commands

From `PlanarAlternatingDeflectorMixer/`:

```bash
tail -n 8 results/verified_flux_sequential_v2/all_samples.csv

find results/verified_flux_sequential_v2 -maxdepth 2 \
    -name objectives.csv -print

pgrep -af 'bayes_optimize_sequential|snakemake|simpleFoam|scalarTransportFoam'
```

For one sample:

```bash
cat results/verified_flux_sequential_v2/00005/objectives.csv
grep 'SIMPLE solution converged' \
    results/verified_flux_sequential_v2/00005/FlowCase/log.simpleFoam
tail -n 50 \
    results/verified_flux_sequential_v2/00005/ScalarTransportCase/log.scalarTransportFoam
```

Main campaign outputs are:

```text
results/verified_flux_sequential_v2/all_samples.csv
results/verified_flux_sequential_v2/pareto_front.png
results/verified_flux_sequential_v2/gp_checkpoint.pt   (once BO fitting begins)
results/verified_flux_sequential_v2/<sample>/objectives.csv
results/verified_flux_sequential_v2/<sample>/FlowCase/pressureDrop.csv
results/verified_flux_sequential_v2/<sample>/ScalarTransportCase/mixing.csv
```

## Known issues and cautions

### Headless field rendering

The former ParaView renderer aborted under headless WSL when it tried to open
an X render window. It has been replaced by a display-independent Python
pipeline: foamlib selects and validates the latest field, Python VTK reads the
reconstructed OpenFOAM mesh, and Pillow rasterizes the top faces and colorbar.
The rule remains non-fatal because a missing PNG is not a failed physical
sample and must not invalidate completed BO objectives.

### Scalar stability failures

One candidate (`00004`) failed only the final-window scalar stability test.
If this repeats, first quantify how many candidates fail and by how much.
Consider an adaptive scalar continuation beyond 600 iterations before
contracting the design space or relaxing tolerances. Any such change alters
the numerical protocol and must be documented and regression-tested.

### Results are local and ignored

Do not assume Git contains the live campaign. Before changing branches,
machines, or clones, archive/copy `results/verified_flux_sequential_v2/`,
including failed sample directories, `all_samples.csv`, and the GP checkpoint.

## Research roadmap and remaining work

### Phase 0: workflow pilot

Operationally complete. Three distinct successful Sobol designs now exist.
The workflow uses corrected second-order scalar transport, explicit flow
convergence, flux-weighted objectives, and exclusion of failed cases.

### Phase 1: space-filling initialization -- in progress

- Complete the remaining 29 successful Sobol designs.
- Review after 8, 16, and 32 successes:
  - meshing/solver/stability failure rates and geometric clustering;
  - duplicate or near-duplicate proposals;
  - mass balance and scalar boundedness;
  - objective distributions;
  - validity of the fixed reference point;
  - concentration of nondominated points at parameter bounds.
- If failures cluster, consider feasibility modelling or documented bound
  contraction. Never use fictitious objective penalties.

### Phase 2: sequential multi-objective BO -- pending

- Run up to 80 successful qLogNEHVI proposals with `q = 1`.
- Every 10 successful BO evaluations, examine Pareto front, dominated
  hypervolume, GP length scales/residuals, proposal novelty, repeatability, and
  bound concentration.
- Fix any early-stopping rule before interpreting the final candidates. The
  current plan suggests stopping only after 15 consecutive practically
  negligible hypervolume gains and low posterior uncertainty along the front.

### Phase 3: numerical verification -- pending

- Select at least five Pareto designs spanning low pressure loss to high
  mixing, plus straight-channel and symmetric-deflector baselines.
- Run at least three systematically refined meshes per selected design.
- Compare `limitedLinear 1` with another bounded high-resolution scheme.
- Quantify mass conservation, scalar boundedness, repeatability, and
  mesh/refinement uncertainty (for example GCI where justified).
- Reassess Pareto dominance using uncertainty intervals.

### Phase 4: scientific comparison and literature survey -- pending

- Start from the user-supplied review:
  “Raza et al. (2020), A Review of Passive Micromixers with a Comparative
  Analysis.” The PDF was supplied outside the repository, so ask for its local
  path if it is not attached in the next session.
- Compare at matched hydraulic diameter, channel length or residence time,
  Reynolds number, diffusivity/Schmidt/Peclet number, inlet condition, and
  mixing-index definition.
- Include at minimum a straight channel, symmetric PADM, and one reproducible
  planar passive-mixer benchmark.
- Use mixing index versus pressure drop or pumping power, not mixing alone.
- Literature values must not be pooled unless their metric definitions and
  operating conditions can be converted defensibly.

### Phase 5: publication gate -- pending

The current best observation is not a demonstrated global optimum and is not
ready for a Chemical Engineering Science global-optimum claim. Required gates
include campaign completion or predeclared defensible stopping, mesh/scheme
independence, quantified uncertainty, fair baselines, a physical mechanism
supported by fields, reproducible successful/failed histories, and preferably
3-D confirmation or experiments for broader passive-mixer claims.

Without these gates, the defensible contribution is a reproducible 2-D
multi-objective numerical methodology and PADM design study, not a universal
passive-mixer optimum or a SAR-mixer result.

## Documentation and slides

The Reveal.js deck documents the renamed device, geometry, parameterization,
objectives, workflow, legacy limitations, revised BO, and research plan:

```bash
cd PlanarAlternatingDeflectorMixer/docs
./serve.sh
```

Then open `http://localhost:8000/`. All repository asset paths in the deck are
relative. The deck uses pinned Reveal.js/KaTeX assets from jsDelivr.

Primary reading order for continuation:

1. `README.md`
2. `PlanarAlternatingDeflectorMixer/QUICKSTART.md`
3. `PlanarAlternatingDeflectorMixer/docs/RESEARCH_PLAN.md`
4. `PlanarAlternatingDeflectorMixer/bayes_optimize_sequential.yaml`
5. `PlanarAlternatingDeflectorMixer/bayes_optimize_sequential.py`
6. `PlanarAlternatingDeflectorMixer/Snakefile`
7. `PlanarAlternatingDeflectorMixer/FlowCase/alternating_deflector_cad.py`

## Recommended first actions for Claude Code

1. Run `git status --short` and preserve every pre-existing user-owned change.
2. Confirm no campaign process is active with `pgrep` before launching one.
3. Source OpenFOAM and run the next bounded WSL evaluation with
   `--max-new-evaluations 1`.
4. Confirm the new sample either has valid `objectives.csv` or a recorded
   failure row with blank objectives.
5. At 8 successful Sobol samples, stop and perform the first Phase-1 audit
   rather than blindly continuing.
6. Keep the user informed of sample IDs, success/failure reason, objectives,
   successful-design count, CPU/rank limits, and whether any process remains
   active.

## Superseding corrected-boundary update (2026-07-26)

This section supersedes every instruction above that treats
`verified_flux_sequential_v2` as physically valid.

The CAD formerly classified boundaries by x-normal. Internal upstream-facing
deflector faces therefore became inlets and downstream-facing faces became
outlets. Velocity and scalar were imposed on internal obstacles, and both
objectives sampled internal faces. All legacy and v2 pressure/mixing values are
invalid and must not seed any corrected GP.

Implemented repairs:

- inlet/outlet classification now requires the complete CAD face at `x=0` or
  `x=L`;
- the CAD emits `geometry_manifest.json` with physical patch bounds and areas;
- `validate_mesh_patches.py` reconstructs ASCII polyMesh boundary geometry;
- full `checkMesh`, `validate_flow_balance.py`, and flow convergence checks run
  before objective aggregation;
- the deflector endpoint floor is 0.02 H, its tessellation is mesh-aware, and
  centre-baffle transitions span several fine cells;
- scalar transport advances in resumable 600-iteration chunks to at most 2400
  and is accepted only on final-50 objective stability;
- result rows include Pa, flow rate, pumping power, mass balance, and mixing
  index; BO pressure is normalized by the corrected straight channel;
- full BO is gated behind twelve corrected Sobol designs and a mixing-index
  threshold of 0.60.

Corrected matched baselines at `Re=10`, `Sc=1000`:

| Baseline | Pressure (Pa) | Pressure ratio | Mixing index |
|---|---:|---:|---:|
| Straight | 2.8737 | 1.000 | 0.1003 |
| Symmetric deflectors | 10.8116 | 3.762 | 0.0901 |
| Strong alternating | 30.5091 | 10.617 | 0.1450 |

The straight CFD pressure differs from the analytical value by 0.039%. The
strong reference exceeds the 20 Pa budget and is far below the 0.60 gate.

Use only this bounded sequence from WSL:

```bash
source /path/to/OpenFOAM-v2506/etc/bashrc
cd /path/to/openfoam-bayesopt-mixer/PlanarAlternatingDeflectorMixer
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python research_sequence.py status
python research_sequence.py next --max-new-evaluations 1
```

As of 2026-07-26, corrected baselines and all twelve screening designs are
complete. Every design passed mesh, flow, scalar-convergence, and rendering
checks; there were no failures. The formal gate is NO-GO: the best mixing
index was 0.1698 at 34.0596 Pa (`DeltaP/DeltaP0 = 11.852`), while the target
was 0.60. The best design below 20 Pa reached only 0.1126. Full optimization
must remain blocked; proceed to topology adaptation. The tracked numerical
summary is `PlanarAlternatingDeflectorMixer/research/corrected_screening_summary.yaml`.

Do not request `research_sequence.py optimization` until the generated
`results/corrected_boundary_v3/screening_gate.json` passes. On a no-go, adapt
the topology—preferably to a genuine 3-D transverse/SAR mechanism—before
spending the 32+80 budget.

## Superseding M10 research direction (2026-07-27)

The active follow-on study is `TwoLayerSerpentineCrossingMixer/`, a genuinely
three-dimensional reconstruction and six-parameter extension of the Hossain
et al. M10 mixer. Do not use PADM observations to train its surrogate. Its BO
contract is explicit coarse/fine multifidelity, strictly sequential (`q=1`),
with at most four OpenFOAM ranks and one Torch thread.

The accepted fine cfMesh Cartesian mesh is 99.967% hexahedral for the
six-unit review geometry. Passive scalar `T` uses
`Gauss linearUpwind gradT` with
`gradT cellLimited pointCellsLeastSquares 1`; the former bounded-scheme values
are retained only for sensitivity comparison.

The separate six-unit, 5.05 mm Raza-review protocol was rerun at Re=1, 10, and
20 with diffusivity `1e-10 m2/s`, the 14 um fine mesh, 2,400 scalar
pseudo-iterations, and four ranks:

| Re | Pressure (Pa) | Area MI | Flux MI | Gate status |
|---:|---:|---:|---:|---|
| 1 | 15.258 | 0.956333 | 0.957138 | pressure pass; mixing fail |
| 10 | 160.829 | 0.917259 | 0.943334 | no tabulated review target |
| 20 | 354.132 | 0.891928 | 0.907025 | pass |

All numerical checks pass. Relative pressure errors at the tabulated points
are -6.39% and -9.20%, within 10%. The area-MI errors are +0.04133 and -0.00907
against an absolute tolerance of 0.03. The follow-up audit established
that the area metric matches the review definition and plausible 80 um lead or
outlet shifts change coarse MI by at most 0.00333. A 30/20/14 um sequence gives
area MI 0.913232/0.942197/0.956333, so the remaining blocker is the
inter-layer aperture topology or its mesh representation.

The CAD manifest now measures the six-unit open interface: 1.12350 mm2 over
the seven vertical segments and 0.557822 mm2 over the six X crossings, for
1.681322 mm2 total. These are full projected overlaps. Verify against an
original model, fabrication mask, or author data before altering them.

The per-design exact-series straight pressure reference is implemented in
`run_case.py`; it is 43.5885 Pa for the four-unit reference at Re=10. Do not
launch the 30-case initialization until the study is explicitly classified as
an exact M10 reproduction or an M10-inspired reconstruction and a small paired
coarse/fine pilot demonstrates useful objective rank correlation.

A controlled fine Re=1 pressure-correction sensitivity reused the exact
330,440-cell mesh and increased `nNonOrthogonalCorrectors` from one to two.
Pressure remained `15.257880624 Pa`; relative L2 differences were `8.72e-7`
for `U` and `8.66e-7` for `phi`, and patch mass balance was `1.95e-9`.
Retain one corrector and do not spend scalar-solver time on this sensitivity.
Pressure coupling is closed as an explanation of the low-Re mixing gap; see
`research/pressure_nonorthogonality_sensitivity.md`.

Read `TwoLayerSerpentineCrossingMixer/research/reproduction_status.md`,
`research/second_order_scalar_validation.md`, and
`research/mesh_qualification.md`, plus `research/review_protocol_audit.md` and
`research/pressure_nonorthogonality_sensitivity.md`, before continuing.
Runtime CFD products below
`results/` are intentionally ignored; the tracked templates, configuration,
tests, and documentation are the reproducibility source of truth.

## OpenFOAM-v2606 snappyHexMesh update (2026-07-28)

This section supersedes the M10 mesh-generator instructions immediately above.

OpenFOAM-v2606 was installed and built from the official source distributions
under `$HOME/OpenFOAM/OpenFOAM-v2606` with its matching
`$HOME/OpenFOAM/ThirdParty-v2606`. `foamSystemCheck` passed, the complete
four-job build finished without logged errors, and the repository-local
function-object library was rebuilt against v2606. Use:

```bash
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd openfoam-bayesopt-mixer
./Allwmake
```

The M10-inspired runner no longer uses cfMesh. Its portable production chain
is `blockMesh`, `surfaceFeatureExtract`, four-rank `snappyHexMesh`,
`reconstructParMesh`, `createPatch`, and independent `checkMesh`. Generated
processor mesh directories are removed after reconstruction. Solver cases
still use four ranks internally, while the BO/campaign remains strictly
sequential (`q=1`). No checkout or OpenFOAM installation path is embedded in
the case templates.

Reference meshes pass at 24 um coarse and 13 um fine. The four-rank nine-unit
reference contains 84,528 cells/93.562% strict hex coarse and 579,096
cells/96.352% strict hex fine. Both contain zero tetrahedra. Oblique fitted
walls create bounded prism/polyhedral cut cells; the fine mesh has one
tet-wedge cell, which is recorded separately and is not a tetrahedron.

The full six-anchor mesh-only preflight passes all 12 coarse/fine cases under
the fixed policy. Evidence is intentionally ignored under
`TwoLayerSerpentineCrossingMixer/results/mesh_qualification/` and is summarized
in `research/mesh_qualification.md`. `multifidelity_pilot.py preflight`
reapplies the current policy to retained meshes, so stale validation JSON is
not trusted blindly.

The earlier cfMesh reproduction and six-unit review results remain historical
scheme/literature evidence, not snappyHexMesh convergence evidence. A complete
four-rank coarse Re=20 smoke test has passed on the new mesh: `493.282 Pa`,
flux MI `0.990867`, mass-balance error `2.96e-9`, final-50 intensity span
`2e-10`, and scalar bounds `[-2.38e-7, 1.0000007]`. This verifies the full
flow/scalar chain and repository-local function objects but is not a fine-grid
convergence result. Continue one pilot observation at a time, only when no
other MPI launcher is active, with:

```bash
cd TwoLayerSerpentineCrossingMixer
python multifidelity_pilot.py status
python multifidelity_pilot.py next --max-new-evaluations 1
```
