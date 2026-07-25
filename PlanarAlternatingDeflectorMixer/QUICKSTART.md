# Quickstart — Planar Alternating-Deflector Micromixer

## Scope

This study models a two-dimensional passive micromixer with a centre baffle
and five repeated pairs of cosine wall deflectors. The stronger deflector
alternates between the top and bottom wall. It bends and stretches one inlet
scalar interface; it is not a three-dimensional split-and-recombine network.

All CAD dimensions are normalized by `H = 1` and converted to SI units by
`scale = 1e-3`.

## Verified six-dimensional parameterization

| BO coordinate | Bounds | CAD mapping |
|---|---:|---|
| `a_weak` | 0–0.12 | weak-wall intrusion; `w_s = 0.5 - a_weak` |
| `a_strong_ratio` | 0–1 | maps into `max(0.12, a_weak+0.04) <= a_strong <= 0.35` |
| `t_s` | 0.04–0.15 | centre-baffle thickness |
| `t_m_ratio` | 0–1 | maps into the mesh-safe merge-baffle interval |
| `L_c` | 0.40–2.40 | deflector interaction length |
| `L_s_ratio` | 0–1 | partitions the remaining cell length into `L_s` and `L_m` |

The realized CAD vector stored in each sample is

```text
(w_s, t_s, t_m, L_s, L_m, delta, k, a_weak, a_strong)
```

where `delta = a_strong - a_weak`, `k = 0`, and
`L_s + L_c + L_m = L_cell = 4`. Fixing `k = 0` removes the old downstream
amplitude trend until evidence supports adding it back.

The transform also enforces two fine mesh cells across `t_m`, one fine cell
across the splitter-thickness step, four fine cells across the split-side gap,
and `0.8 <= L_s,L_m <= 1.8`.

## Objectives

Both objectives are minimized:

```text
J_dp  = <p>_inlet - <p>_outlet                    [m²/s²]
J_mix = flux_weighted_intensity_of_segregation    [-]
```

OpenFOAM uses kinematic pressure, so physical pressure loss is
`DeltaP = rho * J_dp`. For water, a value of `0.0014 m²/s²` is approximately
`1.4 Pa`.

`J_mix` weights every outlet face by positive scalar flux. For communication
with passive-mixer literature, the plots report
`mixing index = 1 - sqrt(J_mix)`. This is the usual relative-standard-deviation
form; it is not the legacy `1 - J_mix` display.

## Sequential BO campaign

The checked-in campaign is deliberately resumable and resource limited:

```text
32 feasible scrambled-Sobol designs
                 ↓
ARD Matérn-5/2 SingleTaskGP on the negated two-objective response
Normalize inputs + standardize outputs
                 ↓
qLogNEHVI, fixed physical reference point
32 restarts, 1024 raw samples, q = 1
                 ↓
80 sequential BO evaluations
```

`q = 1` is validated at startup. Each CFD evaluation uses two MPI ranks and
Torch uses one thread. A failed CFD case is logged with blank targets and is
excluded from GP fitting; it never becomes an artificial penalty observation.

## Run from any clone location

Source an OpenFOAM v2506 environment and use the Python environment containing
the packages in `requirements.txt`. No installation or repository path is
stored in the workflow.

```bash
source /path/to/OpenFOAM-v2506/etc/bashrc
cd /path/to/openfoam-bayesopt-mixer
./Allwmake
cd PlanarAlternatingDeflectorMixer
python bayes_optimize_sequential.py --max-new-evaluations 1
```

Repeat the last command to advance the same campaign by one sequential
evaluation. Omit `--max-new-evaluations` only when the machine may remain
occupied until the configured total of 32 initial plus 80 BO evaluations is
reached. Results and the GP checkpoint are kept below
`results/verified_flux_sequential_v2/` and are not version controlled.

For a single explicitly parameterized workflow, provide the required
Snakemake configuration and retain `--cores 2`.

## Verification status

Two distinct Sobol geometries have completed the revised scalar calculation:

| Pilot | `J_dp` (m²/s²) | Flux `I_s` | `1-sqrt(I_s)` | Scalar stability |
|---|---:|---:|---:|---:|
| first | 0.00141115 | 0.780094 | 0.1168 | final-50 span `3.8e-5` |
| second | 0.00119273 | 0.889416 | 0.0569 | final-50 span `2.0e-6` |

The second geometry exposed a GAMG coarse-grid failure. The checked-in scalar
solver now uses PBiCGStab/DILU with equation relaxation 0.7; rerunning the exact
case then completed all 600 iterations. This is a workflow pilot, not evidence
that the Pareto front has converged.

## Legacy 28-sample campaign

The original eight-Sobol plus twenty-BO result set remains under `results/` for
provenance. It used a different seven-dimensional parameterization,
first-order upwind scalar convection, an unweighted mixing objective, only 200
flow iterations, and penalty values for failures. Its front and reported
`1-I_s` mixing values must not be pooled with or used to validate the revised
campaign.

See `docs/index.html` for the Reveal.js presentation and
`docs/RESEARCH_PLAN.md` for the validation gates required before publication.
