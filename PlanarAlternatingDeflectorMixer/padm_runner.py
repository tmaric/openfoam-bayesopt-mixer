#!/usr/bin/env python3
"""Shared Snakemake invocation for the PADM study.

Both drivers that evaluate a design -- ``bayes_optimize_sequential.py`` and
``run_baselines.py`` -- go through here, so that ``--profile`` means exactly the
same thing to each of them and so that the launch-versus-physics distinction is
made in one place.

THE DISTINCTION THAT MATTERS.  A design that will not mesh or will not converge
is a RESULT: it is recorded as a failure and excluded from the GP, which is
correct.  A missing binary, an unbound path, a broken container or a rejected
sbatch is NOT a result, and recording it as one silently poisons the campaign --
the GP learns that a perfectly good region of the design space is infeasible,
and nothing in the output says otherwise.  ``run_design`` raises
:class:`LaunchFailure` for the second kind so the driver aborts loudly instead
of consuming designs.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import time
from pathlib import Path

import yaml

CASE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = CASE_ROOT.parent
SNAKEFILE = CASE_ROOT / "Snakefile"

# Overridable by --profile or by $PADM_SNAKEMAKE_PROFILE, in that order.
DEFAULT_PROFILE = "profiles/local"

# Rules that cannot fail for a physical reason.  Staging only copies files;
# everything after `hydro_simpleFoam` runs on a case whose physics has already
# succeeded.  If one of these fails, the environment is broken.
_ENVIRONMENT_ONLY_RULES = frozenset(
    {
        "stage_template_cases",
        "hydro_cell_volumes",
        "hydro_decompose",
        "hydro_reconstruct",
        "mixing_decompose",
        "mixing_reconstruct",
        "mixing_post_process",
        "create_foam_files",
        "visualize_sample",
    }
)


# How many CONSECUTIVE designs may fail in the same rule before the campaign is
# treated as environmentally broken rather than as a run of infeasible designs.
#
# The per-design classifier below cannot catch everything: `hydro_mesh` is a
# legitimate physics rule -- cfMesh really does give up on sliver features -- so
# a cfMesh binary that is simply broken looks exactly like a design that cannot
# be meshed.  MEASURED: a local cfMesh install started aborting with
# `free(): invalid pointer` at the smoothing stage on geometry it had meshed in
# 7 s before.  Every design would have been recorded as infeasible, and the GP
# would have learned that the entire space is.  Physics failures scatter; a
# broken tool fails identically every time, and that is the signal.
MAX_CONSECUTIVE_SAME_RULE_FAILURES = 3


class LaunchFailure(RuntimeError):
    """The workflow could not be executed. This is not a design result."""


class FailureStreak:
    """Consecutive failures of one rule, across designs."""

    def __init__(self, limit: int = MAX_CONSECUTIVE_SAME_RULE_FAILURES):
        self.limit = limit
        self.rule: str | None = None
        self.count = 0

    def record_success(self) -> None:
        self.rule = None
        self.count = 0

    def record_failure(self, sample_dir: Path) -> None:
        """Raise LaunchFailure once the same rule has failed `limit` times running."""
        recorded = failed_rules(sample_dir)
        rule = recorded[0][0] if recorded else None
        if rule is not None and rule == self.rule:
            self.count += 1
        else:
            self.rule, self.count = rule, 1
        if rule is not None and self.count >= self.limit:
            raise LaunchFailure(
                f"{self.count} consecutive designs failed in the SAME rule "
                f"({rule}).  Physics failures scatter across the design space; an "
                f"identical failure every time is a broken tool or environment.  "
                f"Check {sample_dir / '.workflow' / 'failed_rule'} and the rule's "
                f"log in the sample directory before recording any of these as "
                f"design results."
            )


def resolve_profile(profile: str | None) -> Path:
    """Absolute path to a workflow profile directory.

    Absolute on purpose: the callers pass ``--directory <sample_dir>`` as well,
    and Snakemake would resolve a relative profile path against a root that is
    not the one the user typed it in.
    """
    raw = profile or os.environ.get("PADM_SNAKEMAKE_PROFILE") or DEFAULT_PROFILE
    path = Path(raw)
    if not path.is_absolute():
        path = CASE_ROOT / path
    path = path.resolve()
    if not (path / "config.yaml").is_file():
        raise LaunchFailure(
            f"no Snakemake profile at {path} (expected {path / 'config.yaml'}); "
            f"available: {sorted(p.name for p in (CASE_ROOT / 'profiles').iterdir())}"
        )
    return path


def _profile_config(profile_dir: Path) -> dict:
    """The profile's own `config:` overrides, as a dict.

    Snakemake takes them as a list of ``key=value`` strings; we only need to
    read them back for the preflight check.
    """
    data = yaml.safe_load((profile_dir / "config.yaml").read_text()) or {}
    entries = data.get("config") or []
    overrides = {}
    for entry in entries:
        key, _, value = str(entry).partition("=")
        overrides[key.strip()] = value.strip()
    return {"_raw": data, **overrides}


def preflight(profile_dir: Path) -> None:
    """Check the execution environment BEFORE consuming any design.

    Cheap, and it converts the most common setup mistakes from "twelve designs
    recorded as infeasible" into one clear message.
    """
    spec = sbatch_spec(profile_dir)

    missing = [tool for tool in ("snakemake",) if shutil.which(tool) is None]
    if missing:
        raise LaunchFailure(f"not on PATH: {', '.join(missing)}")

    if spec is not None:
        # The rules run inside the image on a compute node, so we cannot probe
        # the binaries directly -- check the things whose absence guarantees
        # failure instead.
        # Submission mode.  Deliberately NOT checking for `apptainer` on PATH:
        # it is the COMPUTE NODE that runs `apptainer exec`, from the job script
        # this driver writes.  The driver itself usually runs INSIDE the image,
        # where the host apptainer is neither present nor loadable.
        if shutil.which("sbatch") is None:
            raise LaunchFailure(
                "this profile submits jobs but `sbatch` is not on PATH.  When "
                "the driver runs inside the container the SLURM client must be "
                "bound in -- see run-bo.sbatch and CLUSTER.md."
            )
        image = _resolve_image(spec)
        if not image.is_file():
            raise LaunchFailure(
                f"the container image named by {profile_dir / 'sbatch.yaml'} does "
                f"not exist: {image}.  Build it with apptainer/build.sh."
            )
        for bind in spec.get("binds") or []:
            if not Path(bind).exists():
                raise LaunchFailure(f"bind path does not exist: {bind}")
    else:
        # Local executor: the rules run right here, so probe for real.
        needed = [
            "cartesian2DMesh",
            "simpleFoam",
            "scalarTransportFoam",
            "setExprFields",
            "foamDictionary",
            "decomposePar",
            "reconstructPar",
            "mpirun",
        ]
        absent = [tool for tool in needed if shutil.which(tool) is None]
        if absent:
            raise LaunchFailure(
                "OpenFOAM tooling missing from PATH: "
                + ", ".join(absent)
                + ".  Source an OpenFOAM etc/bashrc, or run inside the Apptainer "
                "image (apptainer exec apptainer/padm.sif ...).  "
                "cartesian2DMesh comes from cfMesh, which is NOT part of a stock "
                "OpenFOAM install."
            )
        wm_options = os.environ.get("WM_OPTIONS")
        if not wm_options:
            raise LaunchFailure("WM_OPTIONS is unset; the OpenFOAM environment is not sourced")
        lib = REPO_ROOT / "platforms" / wm_options / "lib" / "libbayesoptMixerFunctionObjects.so"
        if not lib.is_file():
            raise LaunchFailure(
                f"the study's function-object library is missing: {lib}.  "
                "Build it in THIS environment with ./Allwmake from the repository "
                "root -- pressureDrop and patchMixingQuality produce both objectives, "
                "so without it every design yields an empty objectives.csv."
            )


def failed_rules(sample_dir: Path) -> list[tuple[str, int]]:
    """(rule, exit code) pairs recorded by the Snakefile's EXIT trap."""
    marker = Path(sample_dir) / ".workflow" / "failed_rule"
    if not marker.is_file():
        return []
    out = []
    for line in marker.read_text().splitlines():
        rule, _, code = line.partition("\t")
        if rule:
            out.append((rule.strip(), int(code) if code.strip().isdigit() else -1))
    return out


def classify_failure(sample_dir: Path) -> tuple[str, str]:
    """Return ``(kind, explanation)`` where kind is 'physics' or 'launch'."""
    recorded = failed_rules(sample_dir)
    if not recorded:
        # No rule body ever ran, yet snakemake returned non-zero: the DAG could
        # not be built or no job could be submitted.  This is the case that must
        # never be mistaken for a design result -- a bad container path, an
        # unbound filesystem or a rejected sbatch all land exactly here.
        return (
            "launch",
            "snakemake failed without any rule reaching its shell body "
            "(no .workflow/failed_rule marker): the DAG could not be built or no "
            "job could be submitted",
        )
    names = [rule for rule, _ in recorded]
    environmental = [rule for rule in names if rule in _ENVIRONMENT_ONLY_RULES]
    if environmental:
        return (
            "launch",
            "failed in rule(s) that cannot fail for a physical reason: "
            + ", ".join(environmental),
        )
    return "physics", "failed in rule(s): " + ", ".join(names)


def sbatch_spec(profile_dir: Path) -> dict | None:
    """The profile's `sbatch.yaml`, or None if it runs the DAG inline.

    The presence of the file IS the switch: profiles/local has none and runs
    snakemake in this process; profiles/slurm has one and submits.
    """
    path = Path(profile_dir) / "sbatch.yaml"
    if not path.is_file():
        return None
    return yaml.safe_load(path.read_text()) or {}


def _resolve_image(spec: dict) -> Path:
    image = Path(spec.get("image", "apptainer/padm.sif"))
    return image if image.is_absolute() else (REPO_ROOT / image).resolve()


def _snakemake_command(sample_dir: Path, np_: int, profile_dir: Path,
                       python_bin: str | None, extra_config) -> list[str]:
    cfg = [f"results_dir={sample_dir}", f"np={int(np_)}"]
    if python_bin:
        cfg.append(f"python_bin={python_bin}")
    cfg.extend(extra_config)
    return [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--directory", str(sample_dir),
        "--workflow-profile", str(profile_dir),
        "--config", *cfg,
    ]


def _write_nss(target: Path) -> tuple[Path, Path] | None:
    """Resolve the accounts a container needs, HERE, where NSS works.

    Site accounts are directory-based: `grep $USER /etc/passwd` finds nothing
    while `getent` resolves it.  Inside the image neither sssd nor its libnss
    module exists, so anything that calls getpwuid() dies --

        KeyError: 'getpwuid(): uid not found: 643395244'

    which is exactly what snakemake's own info_header does on startup.  Writing
    the two or three entries that matter is far more robust than plumbing sssd
    sockets and RHEL libnss modules into an Ubuntu image, and it does not rely on
    Apptainer's own passwd synthesis being enabled at the site.
    """
    def _getent(db: str, key: str) -> str:
        out = subprocess.run(["getent", db, key], capture_output=True, text=True)
        return out.stdout if out.returncode == 0 else ""

    user = os.environ.get("USER") or str(os.getuid())
    group = subprocess.run(["id", "-gn"], capture_output=True, text=True).stdout.strip()
    passwd = _getent("passwd", "root") + _getent("passwd", user)
    groups = _getent("group", "root") + (_getent("group", group) if group else "")
    if not passwd.strip():
        return None                       # nothing to add; leave the image alone
    target.mkdir(parents=True, exist_ok=True)
    (target / "passwd").write_text(passwd)
    (target / "group").write_text(groups)
    return target / "passwd", target / "group"


def _job_state(job_id: str) -> tuple[str, int] | None:
    """(state, exit code) once the job is terminal, else None.

    Liveness comes from squeue AND sacct together, never squeue alone: an empty
    squeue during a controller hiccup reads exactly like "the job finished".
    """
    live = subprocess.run(["squeue", "-h", "-j", job_id, "-o", "%T"],
                          capture_output=True, text=True)
    if live.returncode == 0 and live.stdout.strip():
        return None
    acct = subprocess.run(
        ["sacct", "-j", job_id, "-n", "-P", "-X", "-o", "State,ExitCode"],
        capture_output=True, text=True)
    line = (acct.stdout or "").strip().splitlines()
    if not line:
        return None          # sacct has not caught up; keep waiting
    state, _, code = line[0].partition("|")
    state = state.strip().split()[0] if state.strip() else ""
    if state in ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", ""):
        return None
    try:
        rc = int(code.split(":")[0])
    except (ValueError, IndexError):
        rc = 0 if state == "COMPLETED" else 1
    if state != "COMPLETED" and rc == 0:
        rc = 1               # CANCELLED / TIMEOUT / NODE_FAIL are not successes
    return state, rc


def _submit_and_wait(sample_dir: Path, np_: int, profile_dir: Path, spec: dict,
                     python_bin, extra_config, poll_seconds: int = 20):
    """Submit ONE job that runs this design's whole DAG in the container."""
    image = _resolve_image(spec)
    workflow = sample_dir / ".workflow"
    workflow.mkdir(parents=True, exist_ok=True)
    script = workflow / "design.sbatch"
    bind_args = [f"--bind {shlex.quote(str(b))}" for b in (spec.get("binds") or [])]
    nss = _write_nss(workflow / "nss")
    if nss:
        bind_args += [f"--bind {shlex.quote(str(nss[0]))}:/etc/passwd:ro",
                      f"--bind {shlex.quote(str(nss[1]))}:/etc/group:ro"]
    binds = " ".join(bind_args)
    inner = " ".join(shlex.quote(part) for part in
                     _snakemake_command(sample_dir, np_, profile_dir, python_bin, extra_config))

    directives = [
        f"#SBATCH -J {spec.get('job_name', 'padm-design')}",
        f"#SBATCH -A {spec['account']}",
        f"#SBATCH -N {int(spec.get('nodes', 1))}",
        # ONE task with np CPUs, not np tasks.  This is what lets the container's
        # own mpirun see np slots in the job's cpuset with no srun nesting -- and
        # it is why the launcher stays `mpirun` on both backends.
        "#SBATCH -n 1",
        f"#SBATCH --cpus-per-task={int(np_)}",
        f"#SBATCH --mem-per-cpu={int(spec.get('mem_mb_per_cpu', 4000))}",
        f"#SBATCH -t {int(spec.get('runtime_minutes', 480))}",
        f"#SBATCH -o {workflow / 'design.%j.out'}",
        f"#SBATCH -e {workflow / 'design.%j.err'}",
    ]
    directives += [f"#SBATCH {d}" for d in (spec.get("extra_directives") or [])]

    script.write_text(
        "#!/bin/bash\n"
        + "\n".join(directives)
        + "\n"
        # No `set -u`: the image sources OpenFOAM's etc/bashrc, which reads
        # unset variables.
        + "\nset -o pipefail\n"
        # cd into the STUDY directory.  A workflow profile's `configfile:` is
        # resolved against the process cwd, not against the profile directory, and
        # a batch job starts in whatever directory sbatch was called from.  Without
        # this the job dies with
        #   FileNotFoundError: ... 'config/config.yaml'
        # before any rule runs -- which the failure classifier correctly reports as
        # a LAUNCH failure, but which is trivially avoidable.
        + f"cd {shlex.quote(str(CASE_ROOT))} || exit 1\n"
        # --cleanenv: the host environment is inherited by default and $HOME is
        # bind-mounted, so without it a host conda on PATH shadows the image's
        # /opt/venv and an exported `module` function floods every job with
        # "/opt/lmod/.../lmod: No such file or directory".
        + f"\napptainer exec --cleanenv {binds} {shlex.quote(str(image))} {inner}\n"
        + "rc=$?\n"
        + 'echo "[design] snakemake exit=$rc"\n'
        + "exit $rc\n"
    )
    script.chmod(0o755)

    submitted = subprocess.run(["sbatch", "--parsable", str(script)],
                               capture_output=True, text=True)
    if submitted.returncode != 0:
        raise LaunchFailure(
            "sbatch rejected the design job: "
            + (submitted.stderr or submitted.stdout).strip()
        )
    job_id = submitted.stdout.strip().split(";")[0]

    # Record EVERY id we submit.  Cancelling must be done from this list -- on a
    # shared account `scancel -u $USER` destroys other people's work.
    with open(sample_dir.parent / ".padm_jobs", "a") as fh:
        fh.write(f"{job_id}\t{sample_dir.name}\n")
    print(f"[design] {sample_dir.name}: submitted job {job_id} "
          f"(np={np_}, {spec.get('runtime_minutes', 480)} min)", flush=True)

    while True:
        terminal = _job_state(job_id)
        if terminal is not None:
            state, rc = terminal
            print(f"[design] {sample_dir.name}: job {job_id} {state} (exit {rc})",
                  flush=True)
            return subprocess.CompletedProcess([f"sbatch:{job_id}"], rc)
        time.sleep(poll_seconds)


def run_design(
    sample_dir: Path,
    np_: int,
    profile_dir: Path,
    python_bin: str | None = None,
    extra_config: tuple[str, ...] = (),
) -> subprocess.CompletedProcess:
    """Run the per-design DAG once. Does not raise on a physics failure.

    Inline for profiles/local; one sbatch job for a profile carrying sbatch.yaml.
    Either way the DAG, the rules and the groups are identical -- which is the
    whole point: a discrepancy between a laptop run and a cluster run is a real
    result, not a difference in the driver.
    """
    sample_dir = Path(sample_dir)
    # A stale marker from a previous attempt would misclassify this one.
    marker = sample_dir / ".workflow" / "failed_rule"
    if marker.exists():
        marker.unlink()

    spec = sbatch_spec(profile_dir)
    if spec is not None:
        return _submit_and_wait(sample_dir, np_, profile_dir, spec,
                                python_bin, extra_config)
    return subprocess.run(
        _snakemake_command(sample_dir, np_, profile_dir, python_bin, extra_config))
