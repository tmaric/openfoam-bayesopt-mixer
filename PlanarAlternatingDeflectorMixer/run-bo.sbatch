#!/bin/bash
#SBATCH -J padm-bo
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=4000
#SBATCH -t 7-00:00:00
#SBATCH -o padm-bo.%j.out
#SBATCH -e padm-bo.%j.err
#
# ORCHESTRATOR for the PADM Bayesian-optimisation campaign on Lichtenberg.
#
#   cd /work/scratch/tm83tomy/openfoam-bayesopt-mixer/PlanarAlternatingDeflectorMixer
#   sbatch --parsable --export=ALL,N=1 run-bo.sbatch >> .padm_jobs
#
# WHY A BATCH JOB AND NOT `python research_sequence.py &` ON A LOGIN NODE:
# login nodes reap long-lived user processes.  A multi-hour driver started there
# is killed mid-campaign.  This job costs ONE core and 4 GB for its whole life --
# it only runs the BO loop and snakemake, which submit the real CFD work as
# separate jobs.
#
# RESUMING IS FREE.  Snakemake skips designs whose outputs exist and the BO
# driver resumes from results/<campaign>/, so re-submitting after any
# interruption continues where it stopped.
#
# NOTE: no `set -u`.  OpenFOAM's etc/bashrc reads unset variables and the image
# sources it on entry.
set -o pipefail

STUDY_DIR="${STUDY_DIR:-$SLURM_SUBMIT_DIR}"
REPO_ROOT="$(cd "${STUDY_DIR}/.." && pwd)"
SIF="${SIF:-${REPO_ROOT}/apptainer/padm.sif}"
N="${N:-1}"
PROFILE="${PROFILE:-profiles/slurm}"
COMMAND="${COMMAND:-next}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/work/scratch/$USER}"

cd "${STUDY_DIR}" || exit 1
[ -f "${SIF}" ] || { echo "ERROR: no image at ${SIF} -- build it first (apptainer/build.sh)"; exit 1; }

# ---------------------------------------------------------------------------
# NSS for the SLURM client inside the image.
#
# The orchestrator runs the BO driver INSIDE the container, and the driver shells
# out to sbatch/squeue/sacct to submit one job per design and wait for it, so the
# SLURM client has to work in there.  Two MEASURED traps:
#
#   * Binding the HOST /etc/passwd BREAKS it.  Site accounts are directory-based
#     (`grep tm83tomy /etc/passwd` finds nothing, `getent` does), so the bind
#     replaces Apptainer's synthesized entry for the current user with a file
#     that does not contain them -- "squeue: error: Invalid user".
#   * Synthesizing only the current user is not enough either: slurm.conf names
#     SlurmUser=hrzlslurm, and an unresolvable SlurmUser is fatal --
#     "Invalid user for SlurmUser hrzlslurm" / "Unable to process configuration
#     file".
#
# So resolve exactly the accounts SLURM needs on the HOST, where NSS works, and
# bind that minimal file.  Cheaper and far more robust than plumbing the sssd
# sockets and libnss modules into an Ubuntu image.
# ---------------------------------------------------------------------------
NSS_DIR="$(mktemp -d "${SCRATCH_ROOT}/.padm-nss.XXXXXX")"
trap 'rm -rf "${NSS_DIR}"' EXIT
SLURM_USER="$(awk -F= '/^SlurmUser=/{print $2}' /opt/slurm/current/etc/slurm.conf | tr -d ' \r')"
{ getent passwd root
  getent passwd "$USER"
  [ -n "$SLURM_USER" ] && getent passwd "$SLURM_USER"
} > "${NSS_DIR}/passwd" 2>/dev/null
{ getent group root
  getent group "$(id -gn)"
  [ -n "$SLURM_USER" ] && getent group "$SLURM_USER"
} > "${NSS_DIR}/group" 2>/dev/null

BINDS=(
    # Apptainer does not bind /work; without this the case directory is absent.
    --bind "${SCRATCH_ROOT}"
    # Real SLURM client binaries and libs (/shared/bin/sbatch is a site wrapper).
    --bind /opt/slurm
    # slurm.conf: /opt/slurm/current/etc is a SYMLINK to /shared/etc/slurm, which
    # lies outside the bind above -- without this every slurm command dies with
    # "cannot stat file .../slurm.conf".
    --bind /shared/etc/slurm
    # AuthType=auth/munge: without the socket nothing authenticates.
    --bind /run/munge
    --bind "${NSS_DIR}/passwd":/etc/passwd:ro
    --bind "${NSS_DIR}/group":/etc/group:ro
)

# The SLURM client must be on PATH inside the image.  PREPEND_PATH rather than a
# wholesale PATH override: the image's %environment sources OpenFOAM and
# prepends /opt/venv afterwards, and replacing PATH outright would fight both.
#
# The apptainer binary itself is NOT needed in here: the driver writes
# `apptainer exec` into each design's job script, and the COMPUTE NODE runs that
# with the host binary.  (Binding the host one does not work anyway -- it needs
# libsubid.so.3, which an Ubuntu image does not have.)
export APPTAINERENV_SLURM_CONF="${SLURM_CONF:-/opt/slurm/current/etc/slurm.conf}"
export APPTAINERENV_PREPEND_PATH="/opt/slurm/current/bin"

echo "[orchestrator] job ${SLURM_JOB_ID} on $(hostname)"
echo "[orchestrator] study     ${STUDY_DIR}"
echo "[orchestrator] image     ${SIF}"
echo "[orchestrator] profile   ${PROFILE}, ${N} new evaluation(s), command '${COMMAND}'"
echo "[orchestrator] SlurmUser ${SLURM_USER}"

apptainer exec "${BINDS[@]}" "${SIF}" \
    python3 research_sequence.py "${COMMAND}" \
        --max-new-evaluations "${N}" \
        --profile "${PROFILE}"
rc=$?

echo "[orchestrator] research_sequence.py exit=${rc}"
exit "${rc}"
