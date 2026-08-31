# Apptainer image for the PADM study

`padm.def` builds a self-contained environment: OpenFOAM v2512, cfMesh
(`cartesian2DMesh`), CadQuery, BoTorch/GPyTorch, foamlib, VTK, Snakemake and the
Snakemake SLURM executor. The image holds no study data — bind the repository in.

```bash
./apptainer/build.sh            # build here
./apptainer/build.sh --remote   # build on Lichtenberg, rsync the .sif back
```

## Which build path works where — MEASURED

| host | result |
|---|---|
| WSL laptop | **fails.** No `newuidmap`/`newgidmap` and no setuid starter (`~/apptainer` is a user-local install), so both `apptainer build` and `apptainer build --fakeroot` abort with *"newuidmap was not found in PATH ... required with fakeroot and unprivileged installation when user is in /etc/subuid"*. Fix with `sudo apt install uidmap`, or use `--remote`. |
| Lichtenberg login node | **works** unprivileged. `/usr/libexec/apptainer/bin/starter-suid` is setuid; the build reports *"User not listed in /etc/subuid, trying root-mapped namespace"* and proceeds. This is what `--remote` uses. |

`--remote` builds in `/work/scratch/tm83tomy/padm-image`, runs `apptainer test`,
and rsyncs the `.sif` back. The cluster copy is left in place, so the cluster
does not need the laptop's copy pushed back to it.

## After building

The study's OpenFOAM function objects are `dlopen`ed by the container's
OpenFOAM, and `pressureDrop` / `patchMixingQuality` produce **both** BO
objectives — without them every design yields an empty `objectives.csv`. Build
them inside the image:

```bash
apptainer exec --bind "$PWD" apptainer/padm.sif bash -c "./Allwclean && ./Allwmake"
```

> **Run `./Allwclean` first when switching environments.** `wmake` leaves its
> object files and dependency lists in `src/functionObjects/Make/$WM_OPTIONS/`,
> and `$WM_OPTIONS` is `linux64GccDPInt32Opt` for *every* OpenFOAM version — so
> artifacts from a different OpenFOAM look current and are not.  MEASURED: a
> `.dep` carried over from a newer OpenFOAM demanded
> `src/OpenFOAM/lnInclude/FixedList.txx`, a file v2512 does not have, and the
> build stopped with `No rule to make target`.

Both the native and the container build use the same `WM_OPTIONS`
(`linux64GccDPInt32Opt`), so this **overwrites** any host-built copy in
`platforms/`. That is expected: the image is the environment.

## Working inside the image

For anything interactive — and for the whole teaching assignment — enter it once
and stay there:

```bash
apptainer shell --bind "$PWD" apptainer/padm.sif
```

The image carries `sed`, `grep`, `awk`, `git`, `less` and `ffmpeg` alongside
OpenFOAM, cfMesh and the Python stack, so the only step that cannot happen
inside is building the image itself.

## Notes for anyone editing `padm.def`

Four things bit during the builds and are now guarded:

* **`%post` runs under `/bin/sh` (dash on Ubuntu), not bash.** `set -o pipefail`
  is a bashism and aborts the build outright. The code therefore avoids
  depending on any pipe's exit status — note `curl -o file && sh file` rather
  than `curl | bash`, which would feed an empty script to a shell that exits 0
  if the download failed.
* **OpenFOAM's `etc/bashrc` must be sourced under `set +eu`.** It reads unset
  variables *and* its internal helpers `return 1`, so with either option active
  the build dies inside the sourcing with hundreds of trace lines and no
  intelligible error.
* **cfMesh must not land in `FOAM_USER_APPBIN`.** That is
  `$HOME/OpenFOAM/$USER-v2512/...`, and `$HOME` is bind-mounted from the *host*
  at run time — anything installed there during the build simply disappears.
  It goes into the image's own `platforms/$WM_OPTIONS` instead, and the `%post`
  asserts `cartesian2DMesh` is reachable so a bad build fails at build time
  rather than three hours into a campaign.

## If you build cfMesh natively instead

Use the **ESI integration fork**, not the SourceForge repository:

```bash
git clone --depth 1 --branch develop \
    https://develop.openfoam.com/Community/integration-cfmesh.git
```

MEASURED: a SourceForge-built `cartesian2DMesh` on this laptop aborted with
`free(): invalid pointer` during surface smoothing on geometry it had previously
meshed in 7 s. Rebuilt from the ESI fork against the same OpenFOAM-v2512 it
meshes the identical case in 7 s and prints `End`. The image uses the ESI fork
for the same reason.

* **Never source OpenFOAM's `etc/bashrc` from `%environment`.** Apptainer
  sources the environment scripts while the user's command is still in `$@`, and
  OpenFOAM's `etc/config.sh/setup` loops over its positional parameters —
  `eval export`ing anything shaped like `name=value` and **sourcing anything
  that happens to be a file**. MEASURED: with the bashrc sourced at run time,
  `apptainer exec img bash -c "cmd1 && cmd2"` ran the script **twice**, the
  first time without `/opt/venv` on `PATH`; and `apptainer exec img python3
  script.py` printed `source: .../python3: invalid UTF-8 encoding` because
  OpenFOAM had tried to source the interpreter. `%post` now freezes the
  resulting variables into `/opt/openfoam-env.sh` as plain `export` lines, which
  `%environment` sources instead — no argument parsing on the runtime path, and
  a faster container start. There is a build-time assertion that the capture
  produced `WM_PROJECT_DIR`.

`requirements-container.txt` is pinned rather than ranged because CadQuery
*generates the geometry*: a different minor version can change the STL, the mesh
and the objectives with nothing in the output to say why.
