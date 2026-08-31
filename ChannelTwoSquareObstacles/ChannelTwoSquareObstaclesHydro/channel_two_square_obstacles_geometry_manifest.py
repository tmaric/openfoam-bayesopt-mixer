#!/usr/bin/env python3
"""Write the fully expanded geometry manifest without exporting STL."""

from __future__ import annotations

from channel_two_square_obstacles_common import compute_resolved_geometry
from channel_two_square_obstacles_common import load_config
from channel_two_square_obstacles_common import write_geometry_manifest


def main() -> None:
    raw = load_config()
    mode = str(raw.get("cad_mode", "constrained")).strip().lower()
    resolved = compute_resolved_geometry(raw, mode)
    manifest_path = write_geometry_manifest(resolved)
    print(f"Written geometry manifest: {manifest_path}")


if __name__ == "__main__":
    main()
