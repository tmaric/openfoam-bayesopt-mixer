#!/bin/bash
# Serve this folder so the decks can be opened from another device.
#
# NOT required to view them: both decks are standalone HTML with local CSS/JS
# and no CDN, so double-clicking either file works, offline, in any browser.
set -euo pipefail

DOCS_DIR="$(cd "${0%/*}" && pwd)"
PORT="${1:-8000}"

cd "$DOCS_DIR"
echo "Serving on http://localhost:${PORT}/"
echo "  theory   : http://localhost:${PORT}/bayesian-optimization-cfd-theory.html"
echo "  tutorial : http://localhost:${PORT}/bayesian-optimization-cfd-tutorial.html"
python3 -m http.server "$PORT"
