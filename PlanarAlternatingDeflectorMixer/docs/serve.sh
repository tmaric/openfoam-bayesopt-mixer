#!/bin/bash
set -euo pipefail

DOCS_DIR="$(cd "${0%/*}" && pwd)"
PORT="${1:-8000}"

cd "$DOCS_DIR"
python3 -m http.server "$PORT"
