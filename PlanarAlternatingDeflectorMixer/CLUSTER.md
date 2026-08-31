# Running the PADM campaign on a cluster

The same study runs on a laptop, on a cluster login node, and on SLURM compute
nodes. **The only thing that changes is `--profile`.** There are no
per-campaign `.sbatch` scripts and no `if cluster:` branches in the workflow.

| where | profile | how the work runs |
|---|---|---|
| laptop / workstation | `profiles/local` | everything in one shell, `mpirun` |
| shared box, one design at a time | `profiles/local2` | as above, `jobs: 1` |
| cluster login node (smoke test) | `profiles/local` | as above |
| cluster compute nodes | `profiles/slurm` | one sbatch per design, DAG runs in the image on the node |

Everything else — OpenFOAM v2512, cfMesh, CadQuery, BoTorch, Snakemake — lives
in one Apptainer image, so nothing is installed natively on the cluster.

---

## Why the profiles look the way they do

**`np` has exactly one source of truth.** The rank count lives in the launcher
*and* in `system/decomposeParDict`'s `numberOfSubdomains`. If they disagree,
`decomposePar` writes one number of `processor*` directories and the solver is
launched with another. The decompose rules therefore write `numberOfSubdomains`
from `config["np"]` at run time, so the committed dictionary is never read. The
method is **scotch**: it partitions the cell connectivity graph, so there are no
per-direction coefficients whose product must equal `np`, and any `np` works.

**One sbatch job per DESIGN, running the design's whole DAG in the image.**
The driver submits it and waits; `profiles/slurm/sbatch.yaml` describes the
allocation and its mere presence is what switches `padm_runner` from running
snakemake inline to submitting it. The DAG, the rules and the groups are
byte-identical to the laptop run — which is the point: a discrepancy between the
two is then a real result, not a difference in the driver.

*Why not snakemake's own SLURM executor, one job per rule group?* It was tried
first and the grouping worked (the generated job was named
`group_case_pre_hydro_cell_volumes_hydro_decompose_hydro_geometry_hydro_mesh_stage_template_cases`).
It is nonetheless incompatible with running the driver inside the image, which
is what "Apptainer everywhere" means here. **MEASURED:** the executor writes its
own `sys.executable` into every generated job script — inside the container that
is `/opt/venv/bin/python3` — and slurmd then runs that script on the compute
node *outside* any container:

```
/var/spool/slurmd/job54434597/slurm_script: line 4:
    /opt/venv/bin/python3: No such file or directory
```

The path does not exist out there, and the plugin exposes no hook to override
the interpreter. Making it work would mean installing snakemake, its SLURM
plugin, torch and botorch natively on the cluster — the native install this
design exists to avoid.

**The groups still matter.** They are what keeps the sixteen rules from becoming
sixteen serialised process launches, and they are what a future per-rule backend
would submit. The two solves stay **ungrouped**: an MPI solve grouped with
serial steps is submitted as one allocation that does not reliably carry the
solve's rank count, and the launcher then fails with *"More processors requested
than permitted"* or *"not enough slots"*.

**The design job is one task with `np` CPUs on one node** — `-n 1
--cpus-per-task=np`, not `-n np`. That is exactly what lets the container's own
`mpirun -np {np}` see `np` slots in the job's cpuset with no `srun` nesting
anywhere.

**The allocation is sized by the FATTEST rule**, because they all share it: a
scotch `decomposePar` holds the cell connectivity graph on top of the whole mesh
in one process, and `reconstructPar` merges `np` processor meshes plus every
written time. Serial rules do **not** get cheaper as `np` rises.

**MPI is container-internal, so the launcher stays `mpirun`.** The usual cluster
trap — `mpirun` inside Snakemake's one-task `srun` jobstep sees a single slot,
and the `srun --overlap` fix then inherits a one-CPU mask and silently collapses
OpenMPI's busy-wait progress engine — does not arise here, because SLURM never
launches the ranks. The chain is:

```
sbatch job (host) → apptainer exec padm.sif → mpirun -np {np} simpleFoam -parallel
```

`mpi_launcher` is still a config string, so the hybrid form is a one-line
profile change if a site needs it.

### The one check that outranks everything else

OpenFOAM prints `ExecutionTime` (CPU) and `ClockTime` (wall) every step.

> **`ClockTime / ExecutionTime` must be ≈ 1, not ≈ `np`.**

`≈ np` means all the ranks are sharing one core. **There is no error message**
for this — the run is simply slow and the CPU time looks plausible, because
spinning *is* CPU. Check it on the first parallel job in any new environment,
before trusting a single timing or scaling anything up.

---

## Bootstrap on Lichtenberg

```bash
git clone git@github.com:tmaric/openfoam-bayesopt-mixer.git \
    /work/scratch/tm83tomy/openfoam-bayesopt-mixer
cd /work/scratch/tm83tomy/openfoam-bayesopt-mixer
```

Build the image (or `rsync` one built elsewhere to `apptainer/padm.sif`):

```bash
./apptainer/build.sh
```

Build the study's OpenFOAM function-object library **inside the image** — it is
`dlopen`ed by the container's OpenFOAM, and `pressureDrop` / `patchMixingQuality`
produce both BO objectives, so without it every design yields an empty
`objectives.csv`:

```bash
apptainer exec --bind /work/scratch/tm83tomy apptainer/padm.sif \
    bash -c "./Allwclean && ./Allwmake"
```

> **Run `./Allwclean` first when switching environments.** `wmake` leaves its
> object files and dependency lists in `src/functionObjects/Make/$WM_OPTIONS/`,
> and `$WM_OPTIONS` is `linux64GccDPInt32Opt` for *every* OpenFOAM version — so
> artifacts from a different OpenFOAM look current and are not.  MEASURED: a
> `.dep` carried over from a newer OpenFOAM demanded
> `src/OpenFOAM/lnInclude/FixedList.txx`, a file v2512 does not have, and the
> build stopped with `No rule to make target`.

Smoke-test on the **login node** first — this is the cheapest proof the image is
self-sufficient there, and it needs no scheduler at all:

```bash
cd PlanarAlternatingDeflectorMixer
apptainer exec --bind /work/scratch/tm83tomy ../apptainer/padm.sif \
    python3 research_sequence.py next --max-new-evaluations 1 --profile profiles/local
```

## Submitting the campaign

```bash
cd PlanarAlternatingDeflectorMixer
sbatch --parsable --export=ALL,N=1 run-bo.sbatch >> .padm_jobs
```

`run-bo.sbatch` is a **one-core** orchestrator: it runs the BO driver inside the
image, and the driver submits one job per design and waits. It is a batch job
rather than a login-node process because login nodes reap long-lived processes;
a multi-hour driver started there is killed mid-campaign. Resuming is free —
re-submit, and both the BO driver and Snakemake skip what already exists.

Every design job id it submits is appended to `<results>/.padm_jobs`.

Variables it honours: `N` (new evaluations), `PROFILE` (default
`profiles/slurm`), `COMMAND` (`next` / `status` / `optimization`), `SIF`.

### The binds, and why each is there

There are **two different bind lists**, and copying one into the other breaks
things.

**The orchestrator** (`run-bo.sbatch`) runs the BO driver inside the image, and
the driver shells out to `sbatch`/`squeue`/`sacct`, so the SLURM client has to
work in there:

| bind | why |
|---|---|
| `/work/scratch/$USER` | Apptainer does not bind `/work`; without it the case directory is simply absent inside the container |
| `/opt/slurm` | the real client binaries and libs. `/shared/bin/sbatch` is only a site wrapper |
| `/shared/etc/slurm` | **`/opt/slurm/current/etc` is a symlink to it**, so it is *outside* the bind above. Without this every slurm command dies with `cannot stat file .../slurm.conf` |
| `/run/munge` | `AuthType=auth/munge`; without the socket nothing authenticates |
| a synthesized `/etc/passwd`, `/etc/group` | see below |

**The design jobs** run on a compute node and never submit anything, so they
need only the data bind — `--bind /work/scratch/tm83tomy`, from
`binds:` in `sbatch.yaml`. OpenFOAM, cfMesh, Python and `mpirun` are all inside
the image, and MPI never leaves it. Do **not** copy the orchestrator's bind list
here: a design job has no use for the SLURM client, and binding the host
`/etc/passwd` would actively break it (see below).

They are launched with `apptainer exec --cleanenv`. **MEASURED:** without it the
host environment is inherited and `$HOME` is bind-mounted, so a host conda ahead
of `/opt/venv` on `PATH` shadows the image's Python — `import
snakemake_executor_plugin_slurm` failed inside a container that demonstrably had
it — and the host's exported `module` shell function floods every job with
`/opt/lmod/8.7.14/libexec/lmod: No such file or directory`.

#### The `/etc/passwd` trap, measured twice

* **Binding the host `/etc/passwd` breaks the SLURM client.** Site accounts are
  directory-based — `grep tm83tomy /etc/passwd` finds nothing while `getent`
  resolves it — so the bind *replaces* Apptainer's synthesized entry for the
  current user with a file that does not contain them:
  `squeue: error: Invalid user`.
* **Synthesizing only the current user is not enough.** `slurm.conf` names
  `SlurmUser=hrzlslurm`, and an unresolvable SlurmUser is fatal:
  `Invalid user for SlurmUser hrzlslurm` → `Unable to process configuration file`.

`run-bo.sbatch` therefore resolves `root`, `$USER` and the configured
`SlurmUser` **on the host**, where NSS works, writes them to a temporary file
and binds that. Far more robust than plumbing sssd sockets and libnss modules
into an Ubuntu image.

**The design jobs need the same treatment, for a different caller.** They never
touch SLURM, but snakemake resolves the current user on startup and dies the
moment it cannot:

```
File "/usr/lib/python3.12/getpass.py", line 169, in getuser
    return pwd.getpwuid(os.getuid())[0]
KeyError: 'getpwuid(): uid not found: 643395244'
```

`padm_runner._write_nss()` writes `root` and `$USER` into
`<sample>/.workflow/nss/` at submission time and binds them, so this does not
depend on Apptainer's own passwd synthesis being enabled at the site.

Verified inside the image with this recipe: `squeue` lists jobs, and
`sbatch --test-only` authenticates through munge and gets routed by the site's
Lua `job_submit` plugin.

#### What is NOT bound

The **apptainer binary itself**. The orchestrator never runs it — it writes
`apptainer exec` into the design job script, and the compute node runs that with
the host binary. Binding the host one does not work anyway: it needs
`libsubid.so.3`, which an Ubuntu image does not have.

## Cancelling jobs

`run-bo.sbatch` is submitted with `--parsable` so its id lands in `.padm_jobs`.
**Cancel only from that list, or by job name:**

```bash
squeue -u "$USER" -n padm-bo
scancel <id-from-.padm_jobs>
```

**Never `scancel -u $USER` on a shared account.** `special00004` is shared;
that command destroys other people's work, and `sacct` afterwards shows
`CANCELLED by 64+` — the truncated uid of the shared account itself, not an
administrator, which makes it look like a site action rather than your own.

## Verifying, in order

1. `snakemake --workflow-profile profiles/local -n` — does the DAG build?
2. One design locally at `np=2`, in the image.
3. One design on the **login node** with `profiles/local`. This is the cheapest
   proof the image is self-sufficient and needs no scheduler at all.
4. One design on **compute nodes** with `profiles/slurm`. Then, on the job:
   - `sacct -j <id> -o JobID,NCPUS,ReqCPUS,NNodes` — `np` CPUs on **one** node?
   - **`ClockTime/ExecutionTime` from `log.simpleFoam` — ≈ 1 and not ≈ `np`?**
5. Compare step 4's `objectives.csv` against step 3's field by field. Same `np`,
   same decomposition, same image ⇒ they must agree to solver tolerance. A
   difference is a real parallel defect, not a profile problem.
6. Only then scale up.

### Results of this verification, measured

The same design (the `straight` baseline, np=2) run three ways through the same
image:

| where | profile | Δp [Pa] | mixing index | wall/CPU |
|---|---|---|---|---|
| laptop (WSL) | `profiles/local` | 2.873688034 | 0.10028894 | 1.006 |
| cluster login node | `profiles/local` | 2.873688034 | 0.10028338 | 1.027 |
| cluster compute node | `profiles/slurm` | 2.873688034 | 0.10028171 | 1.056 |
| laptop, **native** (no container) | `profiles/local` | 2.873688034 | 0.10028449 | 1.004 |
| *archived v2506 reference* | — | *2.874* | *0.1003* | — |

**Δp is identical to every digit on all three.** The mixing index agrees to five
significant figures; the spread (1.7e-5 relative) is the scalar solve's iterative
tolerance, not a backend difference. Mass balance error, flow rate and the
kinematic pressure drop are bit-identical across all three. So laptop, login node
and compute node are the same computation, and the v2506 → v2512 container move
is sound.

Other checks:

| check | result |
|---|---|
| image self-test (`apptainer test`) | OpenFOAM v2512, `cartesian2DMesh`, and the full Python stack present |
| `./Allwclean && ./Allwmake` in the image | library built against the image's OpenFOAM, on both hosts |
| `sacct -j <design job>` | `NCPUS=2 ReqCPUS=2 NNodes=1 COMPLETED 0:0`, 8m34s |
| `nProcs` reported by the solver | 2 everywhere, matching `np` and the `processor*` count |
| launch-vs-physics classifier | exercised for real: both a missing `config/config.yaml` and a `getpwuid()` failure inside the design container were reported as **launch** failures, not recorded as design results |

A fourth run, **native on the laptop with no container at all** (`profiles/local`
in a sourced OpenFOAM-v2512 shell), gives Δp = 2.873688034 Pa and mixing index
0.10028449 at wall/CPU 1.004 — so the profile mechanism is not container-specific
either.

That native run needed a fix first, and it is worth recording. The laptop's
cfMesh, built from the **SourceForge** `cfmesh` repository, had started aborting
with `free(): invalid pointer` during surface smoothing on geometry the archived
log shows it meshing in 7 s — outside Snakemake too, so it was the binary and not
the workflow. Rebuilding from the **ESI integration fork** that the image uses,
`https://develop.openfoam.com/Community/integration-cfmesh` (branch `develop`),
against the same OpenFOAM-v2512, fixes it: the identical case now meshes in 7 s
and prints `End`. The old binaries are kept at
`$HOME/OpenFOAM/cfmesh-sourceforge-backup/`.

```bash
git clone --depth 1 --branch develop \
    https://develop.openfoam.com/Community/integration-cfmesh.git \
    "$HOME/OpenFOAM/integration-cfmesh"
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
cd "$HOME/OpenFOAM/integration-cfmesh" && ./Allwmake -j "$(nproc)"
```

It installs into `$FOAM_USER_APPBIN`, which is already on `PATH` once OpenFOAM
is sourced. Prefer the ESI fork over the SourceForge one for any current
OpenFOAM: it is the maintained integration and is what the image builds.

A serial pass proves nothing about the code paths that exist only in parallel.
Read the *result*, not the exit code: a deadlock leaves a job alive with a
truncated log and no non-zero return anywhere.

## A diverged design is a result; a launch failure is not

The driver keeps these apart, and the distinction is load-bearing. A design that
will not mesh or will not converge is recorded as a failure and excluded from
the GP — correct. A missing binary, an unbound path, a broken image or a
rejected `sbatch` is **not** a result: recorded as one, it teaches the GP that a
perfectly good region of the design space is infeasible, and nothing in the
output ever says otherwise.

Every rule records its name in `<sample>/.workflow/failed_rule` on a non-zero
exit. `padm_runner.classify_failure()` reads it and raises `LaunchFailure`
— aborting the campaign rather than consuming designs — when the failing rule
cannot fail for a physical reason, or when the marker is absent entirely (no
rule body ever ran, so the DAG could not be built or nothing could be
submitted). `padm_runner.preflight()` runs the same kind of check once, up
front, before the first design is touched.

There is one case the per-design classifier structurally cannot catch, and it is
worth knowing about. `hydro_mesh` is a *legitimate* physics rule — cfMesh really
does give up on sliver features — so a cfMesh binary that is simply **broken**
looks exactly like a design that cannot be meshed. This is not hypothetical: the
laptop's native cfMesh began aborting with `free(): invalid pointer` at the
smoothing stage on geometry it had meshed in 7 s before, which would have marked
every design infeasible. `padm_runner.FailureStreak` covers it — physics
failures scatter across the design space, so **three consecutive failures in the
same rule** are treated as a broken environment and abort the campaign.
