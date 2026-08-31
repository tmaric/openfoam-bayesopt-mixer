# Planar Alternating-Deflector Micromixer — Bayesian optimization study

A two-dimensional passive micromixer optimised with multi-objective Bayesian
optimization over a CAD → mesh → CFD → objectives pipeline. Everything the study
needs is in one Apptainer image, and the execution backend is a Snakemake
workflow profile, so the same study runs on a laptop and on a SLURM cluster.

```
CadQuery geometry → cfMesh (cartesian2DMesh) → simpleFoam → scalarTransportFoam → objectives.csv
                                     ↑                                                    │
                                     └────────── BoTorch qLogNEHVI proposes θ ←────────────┘
```

| | |
|---|---|
| **Design space** | 6 dimensions, mesh-safe (`bayes_optimize_sequential.yaml`) |
| **Objectives** | minimise Δp/Δp<sub>straight</sub> and the flux-weighted intensity of segregation |
| **Operating point** | Re = 10, Sc = 1000, H = 1 mm, L = 24 mm, water |
| **Model** | ARD Matérn-5/2 `SingleTaskGP` per objective, q = 1 |
| **Acquisition** | `qLogNoisyExpectedHypervolumeImprovement` |
| **Status** | corrected screen complete — **NO-GO**, see [Results](#results) |

---

## 1. Get the environment

The image carries OpenFOAM v2512, **cfMesh** (`cartesian2DMesh`, which no
OpenFOAM release ships), CadQuery, BoTorch, foamlib, VTK and Snakemake. Nothing
else has to be installed — not on a laptop, not on a cluster.

```bash
cd ..                      # repository root
./apptainer/build.sh       # or: ./apptainer/build.sh --remote
```

`--remote` builds on Lichtenberg and copies the `.sif` back; use it when the
local Apptainer cannot build (see `../apptainer/README.md`).

Then build the study's OpenFOAM function objects **inside the image**. They are
`dlopen`ed by the container's OpenFOAM, and `pressureDrop` and
`patchMixingQuality` produce *both* BO objectives — without them every design
yields an empty `objectives.csv`:

```bash
apptainer exec --bind "$PWD" apptainer/padm.sif bash -c "./Allwclean && ./Allwmake"
```

> `Allwclean` first is not optional when switching environments. `wmake` keys its
> object files on `$WM_OPTIONS`, which is `linux64GccDPInt32Opt` for *every*
> OpenFOAM version, so a stale tree from a different version looks current and
> is not.

## 2. Run the study

Everything below runs inside the image. The `--profile` flag — and nothing else
— chooses where the CFD happens.

```bash
cd PlanarAlternatingDeflectorMixer

apptainer exec --bind "$PWD/.." ../apptainer/padm.sif \
    python3 research_sequence.py status

apptainer exec --bind "$PWD/.." ../apptainer/padm.sif \
    python3 research_sequence.py next --max-new-evaluations 1 --profile profiles/local
```

`research_sequence.py` is the gated driver. It advances the campaign by **one
bounded step** and stops:

| stage | what runs |
|---|---|
| baselines missing | `run_baselines.py` — straight, symmetric-deflector, strong-alternating |
| baselines done | `bayes_optimize_sequential.py --stage screening` — the 12-design gate |
| screen passed | `research_sequence.py optimization` — must be asked for explicitly |

### Choosing the backend

| `--profile` | where the CFD runs |
|---|---|
| `profiles/local` | laptop, workstation, **or a cluster login node** |
| `profiles/local2` | as above, one design at a time |
| `profiles/slurm` | cluster compute nodes, one sbatch job per design |

`--np` sets the MPI ranks per solve (default 2) and is the single source of truth
for both the launcher and `numberOfSubdomains`. Keep it **equal across the
designs of one campaign**: an MPI job runs at the pace of its slowest rank, so a
varying `np` makes designs incomparable.

On the cluster, submit the orchestrator rather than running the driver on a
login node — login nodes reap long-lived processes:

```bash
sbatch --parsable --export=ALL,N=1 run-bo.sbatch >> .padm_jobs
```

See [CLUSTER.md](CLUSTER.md) for the bind list, the SLURM specifics, and the one
diagnostic that outranks the rest.

### One explicit design, without the BO driver

```bash
apptainer exec --bind "$PWD/.." ../apptainer/padm.sif \
    snakemake --workflow-profile profiles/local --config results_dir=results/manual_00
```

### Rebuild the figures and the animation

```bash
apptainer exec --bind "$PWD/.." ../apptainer/padm.sif \
    snakemake visualize --workflow-profile profiles/local \
        --config results_dir=results/corrected_boundary_v3
```

Writes `pareto_animation.gif` (always) and `pareto_animation.mp4` (ffmpeg) into
`<results_dir>/visualizations/`. This is the animation embedded in the slide deck.

## 3. Running natively, without the container

Supported and verified to give the same numbers. Source an OpenFOAM v2512
`etc/bashrc`, provide the Python packages from `requirements.txt`, and build
cfMesh from the **ESI integration fork** — not the SourceForge repository:

```bash
git clone --depth 1 --branch develop \
    https://develop.openfoam.com/Community/integration-cfmesh.git ~/OpenFOAM/integration-cfmesh
source ~/OpenFOAM/OpenFOAM-v2512/etc/bashrc
cd ~/OpenFOAM/integration-cfmesh && ./Allwmake -j "$(nproc)"
```

A SourceForge-built `cartesian2DMesh` aborted here with `free(): invalid pointer`
during surface smoothing on geometry it had previously meshed in 7 s. The ESI
fork is the maintained integration and is what the image builds.

## 4. Verified equivalence across backends

The `straight` baseline at np = 2, one design, four ways:

| where | profile | Δp [Pa] | mixing index | wall/CPU |
|---|---|---|---|---|
| laptop, native | `profiles/local` | 2.873688034 | 0.10028449 | 1.004 |
| laptop, in image | `profiles/local` | 2.873688034 | 0.10028894 | 1.006 |
| cluster login node | `profiles/local` | 2.873688034 | 0.10028338 | 1.027 |
| cluster compute node | `profiles/slurm` | 2.873688034 | 0.10028171 | 1.056 |
| *archived v2506 reference* | — | *2.874* | *0.1003* | — |

Δp is identical to every digit; the mixing index agrees to five significant
figures (the scalar solve's iterative tolerance). **`ClockTime/ExecutionTime` ≈ 1
everywhere** — the check that catches MPI ranks silently sharing one core, which
produces no error message of its own.

## 5. A diverged design is a result; a launch failure is not

A design that will not mesh or will not converge is recorded as a failure and
excluded from the GP. That is correct. A missing binary, an unbound path or a
rejected `sbatch` is **not** a result: recorded as one it teaches the GP that a
perfectly good region of the design space is infeasible, and nothing in the
output ever says otherwise.

`padm_runner.py` keeps them apart — every rule records its name in
`<sample>/.workflow/failed_rule`, `preflight()` probes the environment before the
first design is touched, and `FailureStreak` aborts when three consecutive
designs fail in the *same* rule (physics failures scatter; a broken tool fails
identically every time).

<a id="results"></a>
## 6. Results — the corrected screen is a NO-GO

Twelve scrambled-Sobol designs, all twelve numerically successful, zero failures:

| | value | gate |
|---|---|---|
| best mixing index | **0.1698** at 34.06 Pa | ≥ 0.60 |
| best under the 20 Pa budget | 0.1126 at 16.47 Pa | — |
| failed fraction | 0.0 | ≤ 0.25 |
| **decision** | **`adapt_topology_before_full_sequential_bo`** | |

The full 32 + 80 campaign is intentionally blocked. This planar topology
stretches and diffuses the inlet interface but does not multiply it, and no
setting of its six parameters reaches a useful mixing index. The finding is the
result — the sibling `TwoLayerSerpentineCrossingMixer` study is the response.

`results/corrected_boundary_v3_baselines/` and `results/corrected_boundary_v3/`
hold the corrected data. The archived `results/` and
`results/verified_flux_sequential_v2/` predate the boundary repair and **must not
be pooled** with it.

## 7. Using this as a teaching study

Part 5 of the slide deck is a hands-on assignment: six exercises that alternate
**task → answer**, with the exact commands on the answer slide. Tasks 1–4 run
CFD (minutes each); tasks 5–6 run none at all, reusing the twelve designs the
campaign already paid for.

| # | Task | Learning outcome |
|---|---|---|
| 1 | Run the straight baseline; check it against 12νUL/H² | the workflow, the container, the profile switch |
| 2 | Find which YAML keys the CAD actually reads | read the source, not the field names |
| 3 | Halve the asymmetry; predict, then measure | parameters → geometry → objectives; the trade-off |
| 4 | Delete outputs; predict what Snakemake re-runs | the DAG is target-driven, not a file watcher |
| 5 | Fit a GP to the finished designs; rank the parameters | GP + ARD as a free sensitivity analysis |
| 6 | Maximise UCB for κ = 0, 2, 10 | exploration vs exploitation, on real data |

Every number on the answer slides was produced by running the exercise. The
reference solution for tasks 5 and 6 is `docs/exercises/surrogate.py`:

```bash
apptainer exec --bind "$PWD/.." ../apptainer/padm.sif python3 docs/exercises/surrogate.py
```

Serve the deck with `./docs/serve.sh` and open <http://localhost:8000/>.

## Repository map

| path | what |
|---|---|
| `Snakefile` | the per-design DAG: geometry → mesh → flow → scalar → objectives |
| `padm_runner.py` | shared Snakemake invocation, profile resolution, failure classification |
| `research_sequence.py` | gated one-step driver (`status` / `next` / `optimization`) |
| `bayes_optimize_sequential.py` | the BO loop: GP fit, acquisition, Pareto bookkeeping |
| `run_baselines.py` | the three matched baselines and the analytic straight-channel check |
| `config/config.yaml` | backend-independent defaults (`np`, launcher, preamble) |
| `profiles/` | `local`, `local2`, `slurm` |
| `FlowCase/`, `ScalarTransportCase/` | OpenFOAM case templates |
| `FlowCase/alternating_deflector_cad.py` | the CadQuery geometry generator |
| `validate_*.py`, `verify_*.py` | independent mesh, flow and convergence checks |
| `docs/index.html` | the slide deck, incl. the hands-on assignment (`./docs/serve.sh`) |
| `docs/exercises/` | reference solutions for the assignment |
| [QUICKSTART.md](QUICKSTART.md) | parameterization, objectives, campaign ladder |
| [CLUSTER.md](CLUSTER.md) | SLURM, binds, and the measured failure modes |
| [docs/RESEARCH_PLAN.md](docs/RESEARCH_PLAN.md) | phase gates and publication criteria |
