#!/usr/bin/env bash
#
# Build apptainer/padm.sif.
#
#   ./apptainer/build.sh              # build here
#   ./apptainer/build.sh --remote     # build on Lichtenberg, fetch the .sif back
#
# WHICH PATH DO YOU HAVE?  An unprivileged Apptainer build needs one of:
#   (a) newuidmap/newgidmap  -- Ubuntu/Debian package `uidmap`, plus a subuid
#       range for your user (this machine already has /etc/subuid), or
#   (b) a setuid Apptainer installation (`starter-suid`), which is what the
#       system packages on RHEL-family clusters install.
# A user-local Apptainer under $HOME has neither by default, which is why
# --remote exists: Lichtenberg's /usr/libexec/apptainer/bin/starter-suid is
# setuid and builds there work unprivileged.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEF="${REPO_ROOT}/apptainer/padm.def"
SIF="${REPO_ROOT}/apptainer/padm.sif"

REMOTE_HOST="${PADM_REMOTE_HOST:-tm83tomy@lcluster1.hrz.tu-darmstadt.de}"
REMOTE_DIR="${PADM_REMOTE_DIR:-/work/scratch/tm83tomy/padm-image}"

die() { echo "ERROR: $*" >&2; exit 1; }

build_remote() {
    echo "==> building on ${REMOTE_HOST}:${REMOTE_DIR}"
    ssh -n "${REMOTE_HOST}" "mkdir -p '${REMOTE_DIR}/tmp' '${REMOTE_DIR}/cache' '${REMOTE_DIR}/apptainer'"
    # The def file reads apptainer/requirements-container.txt through %files, so
    # the remote build needs the same relative layout.
    scp -q "${DEF}" "${REMOTE_HOST}:${REMOTE_DIR}/apptainer/padm.def"
    scp -q "${REPO_ROOT}/apptainer/requirements-container.txt" \
           "${REMOTE_HOST}:${REMOTE_DIR}/apptainer/requirements-container.txt"
    ssh -n "${REMOTE_HOST}" "
        set -eu
        cd '${REMOTE_DIR}'
        export APPTAINER_TMPDIR='${REMOTE_DIR}/tmp'
        export APPTAINER_CACHEDIR='${REMOTE_DIR}/cache'
        apptainer build --force apptainer/padm.sif apptainer/padm.def
        apptainer test apptainer/padm.sif
    "
    echo "==> fetching the image"
    rsync -avh --progress "${REMOTE_HOST}:${REMOTE_DIR}/apptainer/padm.sif" "${SIF}"
    echo "==> ${SIF}"
    echo "    the cluster copy stays at ${REMOTE_DIR}/apptainer/padm.sif"
}

build_local() {
    command -v apptainer >/dev/null 2>&1 \
        || die "apptainer is not on PATH (try: export PATH=\"\$HOME/apptainer/bin:\$PATH\")"

    if ! command -v newuidmap >/dev/null 2>&1 \
       && [ ! -u "$(dirname "$(command -v apptainer)")/../libexec/apptainer/bin/starter-suid" ] 2>/dev/null \
       && [ ! -u /usr/libexec/apptainer/bin/starter-suid ]; then
        cat >&2 <<'MSG'
ERROR: this Apptainer cannot build unprivileged here.

  newuidmap/newgidmap are absent and there is no setuid starter, so both
  `apptainer build` and `apptainer build --fakeroot` will fail with
  "newuidmap was not found in PATH".

  Fix it with EITHER of:

    sudo apt install uidmap          # then re-run this script

    ./apptainer/build.sh --remote    # build on the cluster, rsync the .sif back
MSG
        exit 1
    fi

    echo "==> building ${SIF}"
    APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${TMPDIR:-/tmp}}" \
        apptainer build --force "${SIF}" "${DEF}"
    apptainer test "${SIF}"
    echo "==> ${SIF}"
}

case "${1:-}" in
    --remote) build_remote ;;
    ""|--local) build_local ;;
    *) die "usage: $0 [--local|--remote]" ;;
esac

cat <<MSG

Next: build the study's OpenFOAM function-object library IN this environment.
It is dlopen'ed by the container's OpenFOAM, so a host-built copy of it will not
necessarily load -- and both builds use the same WM_OPTIONS directory, so this
overwrites the host one.

    apptainer exec --bind "${REPO_ROOT}" "${SIF}" ${REPO_ROOT}/Allwmake
MSG
