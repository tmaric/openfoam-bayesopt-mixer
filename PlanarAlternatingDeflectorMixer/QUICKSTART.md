# Quickstart — Planar Alternating-Deflector Micromixer

## Scope

This study models a two-dimensional passive micromixer with a centre baffle
and five repeated pairs of cosine wall deflectors. It is not a three-dimensional
split-and-recombine network.

## Corrected-boundary prerequisite

The former campaigns are invalid because the CAD generator classified every
x-normal obstacle face as an inlet or outlet. The corrected generator requires
the complete inlet face to lie at `x=0` and the complete outlet face to lie at
`x=L`. Every generated mesh is independently checked for patch location, area,
mesh quality, flow rate, and mass balance before objectives are accepted.

The archived `results/` and `results/verified_flux_sequential_v2/` values must
not be pooled with the corrected campaign.

## Six-dimensional parameterization

| BO coordinate | Bounds | CAD mapping |
|---|---:|---|
| `a_weak` | 0–0.12 | weak-wall intrusion; `w_s = 0.5 - a_weak` |
| `a_strong_ratio` | 0–1 | maps into `max(0.12, a_weak+0.04) <= a_strong <= 0.35` |
| `t_s` | 0.04–0.15 | centre-baffle thickness |
| `t_m_ratio` | 0–1 | maps into the mesh-safe merge-baffle interval |
| `L_c` | 0.40–2.40 | deflector interaction length |
| `L_s_ratio` | 0–1 | partitions the remaining cell length into `L_s` and `L_m` |

The CAD now uses a 0.02 H endpoint floor so every deflector end spans more
than one fine cell. Its mesh-safety algebra includes the floor on both walls.
The merge-to-split thickness transition is resolved over four CAD minimum
feature lengths, and the cosine tessellation is commensurate with the CFD mesh.

## Objectives and retained metrics

Both BO objectives are minimized:

```text
J_dp  = DeltaP / DeltaP_straight                  [-]
J_mix = flux_weighted_intensity_of_segregation    [-]
```

Every objective row also retains kinematic and dimensional pressure drop,
flow rate, pumping power, mass-balance error, and the literature-style mixing
index `1-sqrt(J_mix)`.

## Corrected research sequence

The default command advances at most one CFD evaluation and never runs cases
in parallel:

```bash
python research_sequence.py status
python research_sequence.py next --max-new-evaluations 1
```

The sequence is gated:

1. validated straight, symmetric-deflector, and strong-alternating baselines;
2. twelve corrected scrambled-Sobol screening designs;
3. a no-go unless the best mixing index reaches 0.60 and failures stay below
   25%;
4. only after a pass, explicit `python research_sequence.py optimization` to
   expand to 32 Sobol designs and then 80 strictly sequential BO evaluations.

The corrected baseline values at `Re=10`, `Sc=1000` are:

| Baseline | Pressure (Pa) | `DeltaP/DeltaP0` | Mixing index |
|---|---:|---:|---:|
| Straight | 2.874 | 1.000 | 0.1003 |
| Symmetric deflectors | 10.812 | 3.762 | 0.0901 |
| Strong alternating | 30.509 | 10.617 | 0.1450 |

The straight pressure result differs from the fully developed parallel-plate
solution by only 0.039%. The strong alternating reference already exceeds the
predeclared 20 Pa budget and remains far below the 0.60 mixing gate.

Corrected screen result (2026-07-26): all twelve designs completed and passed
every numerical validation, with no failures. The formal decision is NO-GO.
The best mixing index was only 0.1698 at 34.060 Pa (pressure ratio 11.852),
and the best design below 20 Pa reached only 0.1126 at 16.469 Pa. Both are far
below the predeclared 0.60 continuation threshold. The full 32+80 BO stage
must not be run on this topology.

## Sequential BO configuration

```text
12-design feasibility gate
                 ↓
32 total feasible scrambled-Sobol designs
                 ↓
ARD Matérn-5/2 SingleTaskGP on the negated two-objective response
Normalize inputs + standardize outputs
                 ↓
qLogNEHVI, 32 restarts, 1024 raw samples, q = 1
                 ↓
80 sequential BO evaluations
```

Scalar transport advances in sequential 600-iteration chunks up to 2400 and
stops at the first chunk where the final 50 outlet measurements satisfy the
fixed stability tolerance. Algebraic residual convergence alone is not used
as a substitute for objective stability.

## Run from any clone location

Everything the study needs -- OpenFOAM v2512, cfMesh (`cartesian2DMesh`, which
is **not** part of a stock OpenFOAM install), CadQuery, BoTorch, Snakemake -- is
in one Apptainer image, so a clone plus the image is a complete environment.

```bash
cd /path/to/openfoam-bayesopt-mixer
./apptainer/build.sh                 # or --remote to build on the cluster

# the study's OpenFOAM function objects must be built IN this environment:
# they are dlopen'ed by the container's OpenFOAM and produce both objectives
apptainer exec --bind "$PWD" apptainer/padm.sif bash -c "./Allwclean && ./Allwmake"

cd PlanarAlternatingDeflectorMixer
apptainer exec --bind "$PWD/.." ../apptainer/padm.sif \
    python3 research_sequence.py next --max-new-evaluations 1 \
        --profile profiles/local
```

Without the container, source an OpenFOAM `etc/bashrc` that has cfMesh built
into its `FOAM_USER_APPBIN`, run `./Allwmake`, then invoke
`research_sequence.py` directly; the thread-pinning the image sets for you is
`export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

### Choosing where the work runs

The backend is a profile, and nothing else changes:

| `--profile` | where |
|---|---|
| `profiles/local` | laptop, workstation, or a cluster **login node** |
| `profiles/local2` | as above but one design at a time (shared box; or when wall time is part of the result) |
| `profiles/slurm` | cluster **compute nodes**, one sbatch per design |

`--np` sets the MPI ranks per CFD solve (default 2). Keep it EQUAL across the
designs of one campaign: an MPI job runs at the pace of its slowest rank, so a
varying `np` makes designs incomparable.

See `CLUSTER.md` for the SLURM path, the bind list, and the one diagnostic that
outranks the rest (`ClockTime/ExecutionTime` must be ~1, never ~np).

All stored paths are clone-relative. Corrected results are written below
`results/corrected_boundary_v3_baselines/` and
`results/corrected_boundary_v3/`; runtime results are not version controlled.

See `docs/index.html` for the Reveal.js presentation and
`docs/RESEARCH_PLAN.md` for the publication decision gates.
