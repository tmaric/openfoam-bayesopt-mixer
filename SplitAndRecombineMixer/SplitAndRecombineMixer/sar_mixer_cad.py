#!/usr/bin/env python3
"""
SAR Lamination Ladder Mixer - CADQuery geometry script.

Generates the 2D mixer fluid-domain boundary (extruded to a thin slab) and
exports it as a single ASCII STL file for use with cfMesh cartesian2DMesh:

  constant/triSurface/sar_mixer.stl

The STL contains four named solid regions that cfMesh turns into patches:
  inlet       – face at x = 0
  outlet      – face at x = TOTAL_L
  walls       – channel top/bottom walls + all internal obstacle surfaces
  frontAndBack – slab faces at z = 0 and z = DEPTH  (type empty in 2-D)

Parameter names follow the sketch in:
  docs/obsidian/06 Resources/sar_lamination_ladder_mixer_sketch.py
"""

import os
import cadquery as cq
import yaml

# ---------------------------------------------------------------------------
# Geometry parameters – loaded from sar_mixer_cad.yaml
# All values in the YAML are in normalised units (H_norm = 1).
# SCALE converts to SI metres.
# ---------------------------------------------------------------------------
_yaml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sar_mixer_cad.yaml")
with open(_yaml_path) as _f:
    _p = yaml.safe_load(_f)

SCALE  = _p["scale"]      # 1 normalised unit → SI metres

H      = _p["H"]      * SCALE
L0     = _p["L0"]     * SCALE
N      = _p["N"]               # integer – no scaling
L_cell = _p["L_cell"] * SCALE
L_s    = _p["L_s"]    * SCALE
L_m    = _p["L_m"]    * SCALE

w_s    = _p["w_s"]   * H
t_s    = _p["t_s"]   * H
t_m    = _p["t_m"]   * H

delta  = _p["delta"] * H

# Extrusion depth (thin slab for 2-D OpenFOAM simulation).
# Must satisfy span_z / span_x < ~0.001 for cfMesh cartesian2DMesh to
# classify the surface as 2D.  With TOTAL_L = 24e-3 m the limit is
# DEPTH < 2.4e-5 m.  0.01*H = 1e-5 m → ratio = 4.17e-4 (safe).
DEPTH  = 0.01 * H

L_c = L_cell - L_s - L_m
assert L_c > 0, "L_cell must be greater than L_s + L_m"

h_d = 0.5 * H - w_s   # deflector intrusion height from each wall
assert h_d > 0, "Need w_s < H/2"

TOTAL_L = 2 * L0 + N * L_cell

# ---------------------------------------------------------------------------
# Helper: cosine-envelope deflector profile as a list of (x, y) points
# ---------------------------------------------------------------------------
N_PTS = 120   # points along the cosine curve

def cosine_bump_points(x_start, x_end, amp, from_top=False, bias=0.0):
    """
    Return a polygon (list of (x,y)) for a cosine-shaped deflector.

    The polygon follows the cosine surface from x_start to x_end, then
    closes back along the wall (y=0 for bottom, y=H for top).

    Closing wall points are added only when the cosine curve endpoints do
    NOT already lie on the wall.  The cosine envelope is exactly 0 at both
    ends, so for the bottom deflector (bias=0) both endpoints are at y=0
    and NO closing points are needed – make_extruded_polygon's close() call
    draws the wall-return segment automatically.  For the top deflector with
    bias=delta>0 the endpoints sit at y=H-delta, so explicit wall points at
    y=H are required.
    """
    import math
    Lc = x_end - x_start
    wall_y = H if from_top else 0.0
    pts = []
    for i in range(N_PTS + 1):
        xi = i * Lc / N_PTS
        env = 0.5 * (1.0 - math.cos(2.0 * math.pi * xi / Lc))
        bval = amp * env
        x = x_start + xi
        if from_top:
            y = H - min(bval + bias, H)
        else:
            y = bval
        pts.append((x, y))
    # Only add explicit closing wall points when the curve does not already
    # end on the wall (avoids zero-length edges that OCC rejects).
    if abs(pts[-1][1] - wall_y) > 1e-10:
        pts.append((x_end, wall_y))
    if abs(pts[0][1] - wall_y) > 1e-10:
        pts.append((x_start, wall_y))
    return pts


# ---------------------------------------------------------------------------
# Build solid obstacle shapes
# ---------------------------------------------------------------------------
def make_extruded_polygon(points_2d, depth):
    """Extrude a closed 2-D polygon into a solid slab."""
    wire = cq.Workplane("XY").polyline(points_2d).close()
    return wire.extrude(depth)


def rect_points(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


# ---------------------------------------------------------------------------
# Collect obstacle solids per unit cell
# ---------------------------------------------------------------------------
obstacle_solids = []

for k in range(N):
    xk = L0 + k * L_cell

    # --- Split splitter (thin rectangle centred in channel) ---
    xs0, xs1 = xk, xk + L_s
    y0s = (H - t_s) / 2.0
    y1s = (H + t_s) / 2.0
    pts = rect_points(xs0, xs1, y0s, y1s)
    obstacle_solids.append(make_extruded_polygon(pts, DEPTH))

    # --- Bottom cosine deflector ---
    xc0, xc1 = xk + L_s, xk + L_s + L_c
    pts_bot = cosine_bump_points(xc0, xc1, h_d, from_top=False, bias=0.0)
    obstacle_solids.append(make_extruded_polygon(pts_bot, DEPTH))

    # --- Top cosine deflector (with delta bias) ---
    pts_top = cosine_bump_points(xc0, xc1, h_d, from_top=True, bias=delta)
    obstacle_solids.append(make_extruded_polygon(pts_top, DEPTH))

    # --- Merge splitter (thin rectangle centred in channel) ---
    xm0, xm1 = xk + L_s + L_c, xk + L_cell
    y0m = (H - t_m) / 2.0
    y1m = (H + t_m) / 2.0
    pts = rect_points(xm0, xm1, y0m, y1m)
    obstacle_solids.append(make_extruded_polygon(pts, DEPTH))


# ---------------------------------------------------------------------------
# Channel outer box and fluid domain
# ---------------------------------------------------------------------------
channel_box = (
    cq.Workplane("XY")
    .box(TOTAL_L, H, DEPTH)
    .translate((TOTAL_L / 2.0, H / 2.0, DEPTH / 2.0))
)

print("Building fluid domain (channel minus obstacles) ...")
fluid = channel_box
for obs in obstacle_solids:
    fluid = fluid.cut(obs)

# ---------------------------------------------------------------------------
# Face classification by outward normal direction
# ---------------------------------------------------------------------------
_TOL = 1e-6

def _is_inlet(f):   n = f.normalAt(); return n.x < -(1.0 - _TOL)
def _is_outlet(f):  n = f.normalAt(); return n.x >  (1.0 - _TOL)
def _is_wall(f):
    n = f.normalAt()
    return not (abs(n.x) > (1.0 - _TOL)) and not (abs(n.z) > (1.0 - _TOL))

# ---------------------------------------------------------------------------
# Named-solid STL export for cfMesh (solid name → OpenFOAM patch name)
# ---------------------------------------------------------------------------
# exporters.export() dropped ascii=True in newer CadQuery; write ASCII STL
# directly from the shape tessellation to avoid the dependency.

out_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "constant", "triSurface"
)
os.makedirs(out_dir, exist_ok=True)


def _stl_block(faces, solid_name, tolerance=1e-4, angular_tol=0.1):
    """Tessellate a face list and return a named ASCII STL solid block."""
    if not faces:
        print(f"  WARNING: no faces found for patch '{solid_name}'")
        return ""
    compound = cq.Compound.makeCompound(faces)
    verts, tris = compound.tessellate(tolerance, angular_tol)
    lines = [f"solid {solid_name}"]
    for tri in tris:
        v0, v1, v2 = verts[tri[0]], verts[tri[1]], verts[tri[2]]
        # Cross product for outward face normal (no CadQuery arithmetic needed)
        e1x, e1y, e1z = v1.x - v0.x, v1.y - v0.y, v1.z - v0.z
        e2x, e2y, e2z = v2.x - v0.x, v2.y - v0.y, v2.z - v0.z
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        nl = (nx*nx + ny*ny + nz*nz) ** 0.5
        if nl > 1e-30:
            nx, ny, nz = nx / nl, ny / nl, nz / nl
        lines += [
            f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}",
            f"    outer loop",
            f"      vertex {v0.x:.6e} {v0.y:.6e} {v0.z:.6e}",
            f"      vertex {v1.x:.6e} {v1.y:.6e} {v1.z:.6e}",
            f"      vertex {v2.x:.6e} {v2.y:.6e} {v2.z:.6e}",
            f"    endloop",
            f"  endfacet",
        ]
    lines.append(f"endsolid {solid_name}")
    return "\n".join(lines) + "\n"


all_faces = fluid.val().Faces()

# NOTE: frontAndBack (z-normal) faces are excluded from the STL.
# cartesian2DMesh requires that ALL face normals lie in a single plane
# (checked via covariance-matrix eigenvalue == 0).  Z-normal faces break
# that check.  cfMesh reads z_min / z_max from vertex coords and creates
# the frontAndBack patch automatically.
patch_faces = {
    "inlet":  [f for f in all_faces if _is_inlet(f)],
    "outlet": [f for f in all_faces if _is_outlet(f)],
    "walls":  [f for f in all_faces if _is_wall(f)],
}

print("Exporting patch STL solids ...")
stl_parts = []
for name, faces in patch_faces.items():
    print(f"  Patch '{name}': {len(faces)} face(s)")
    stl_parts.append(_stl_block(faces, name))

stl_path = os.path.join(out_dir, "sar_mixer.stl")
with open(stl_path, "w") as fh:
    fh.write("".join(stl_parts))
print(f"  Written: {stl_path}")

print()
print("Geometry summary:")
print(f"  Total length  L    = {TOTAL_L:.4f}")
print(f"  Channel height H   = {H:.4f}")
print(f"  Extrusion depth    = {DEPTH:.4f}")
print(f"  Unit cells    N    = {N}")
print(f"  Unit-cell len      = {L_cell:.4f}")
print(f"  Splitter thick t_s = {t_s:.4f}")
print(f"  Splitter thick t_m = {t_m:.4f}")
print(f"  Deflector height h_d = {h_d:.4f}")
print(f"  Shuffle offset delta = {delta:.4f}")
print(f"  Interaction length L_c = {L_c:.4f}")
print("Done.")
