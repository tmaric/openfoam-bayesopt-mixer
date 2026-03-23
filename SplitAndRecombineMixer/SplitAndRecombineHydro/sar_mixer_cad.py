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
# Small overshoot applied to all obstacle extrusions so that no obstacle face
# is coplanar with a channel-box face.  OpenCASCADE Boolean cuts fail (null
# TopoDS_Shape) when the tool and workpiece share an exactly coincident face.
_EPS = DEPTH * 0.01

L_c = L_cell - L_s - L_m
h_d = 0.5 * H - w_s   # deflector intrusion height from each wall

TOTAL_L = 2 * L0 + N * L_cell

# ---------------------------------------------------------------------------
# Geometry validation  (runs before any OCC operation)
# ---------------------------------------------------------------------------
# Minimum feature size: 1 % of H.  For the default SCALE=1e-3 and H=1.0 this
# is 1e-5 m (10 µm), which equals the wall-refinement cell size in meshDict
# and is consistent with TM_MARGIN=0.01 used in bayes_optimize_sequential.py.
_MESH_MIN = 0.01 * H

# x-gap inserted between each merge-splitter outlet (x = xk) and the
# following split-splitter inlet (x = xk + DELTA_X).  Without this gap,
# both features share the exact face x = xk, creating a compound 270°
# re-entrant corner where 3 STL triangles meet at one edge.  cfMesh
# resolves that topology by duplicating a mesh vertex → zero-area faces.
# Separating the two outlets by DELTA_X > 4 fine cells gives cfMesh two
# independent single re-entrant corners, which it handles cleanly.
DELTA_X = 4 * _MESH_MIN   # 4e-5 m  <<  L_s_min ≈ 8e-4 m

def _check_geometry() -> None:
    """Validate all CAD parameters before any OCC operation is attempted.

    Raises ValueError listing every violated constraint so the caller gets a
    diagnostic message instead of a cryptic 'Null TopoDS_Shape object' error
    from deep inside an OpenCASCADE Boolean operation.

    Each check is labelled Gn and maps to a concrete failure mode:

    G1  L_c >= _MESH_MIN
        Interaction region must have positive x-extent.
        Violated → 'assert L_c > 0' fires, or degenerate cosine polygon.

    G2  h_d >= _MESH_MIN
        Deflector intrusion height h_d = H/2 - w_s must be positive.
        Near-zero h_d → cosine bump polygon has near-zero area → OCC error.

    G3a t_s >= _MESH_MIN
    G3b t_m >= _MESH_MIN
        Both splitters must be thick enough for at least one mesh cell.

    G4  t_s - t_m >= _MESH_MIN
        Split splitter must be strictly wider than the merge splitter so
        that the step at each cell boundary has positive y-extent.
        The split and merge splitters are separated by DELTA_X = 4·_MESH_MIN
        in x, so their STL boundaries never coincide.  Without this gap,
        the compound 270° re-entrant corner (merge-end + split-start at the
        same x) causes cfMesh to create duplicate mesh vertices → zero-area
        faces.  The underlying rule: adjacent obstacle boundaries that are
        within a few mesh cells of each other in y must be separated by
        ≥ DELTA_X in x.  t_s > t_m ensures the step has real y-extent.

    G5  w_s - 0.5*t_s >= _MESH_MIN
        Minimum fluid gap in the split section between the deflector peak
        and the nearest splitter surface.
        Bottom deflector peak y = h_d = H/2 - w_s;
        split-splitter bottom face y = (H - t_s)/2.
        Gap = w_s - t_s/2.  Near-zero → channel pinches → meshing fails.

    G6  2*w_s - delta >= 3*_MESH_MIN
        Minimum fluid gap between top and bottom deflectors at the peak of
        the cosine in the interaction region.
        The cosine bval has a floor of _MESH_MIN (see cosine_bump_points), so
        the effective bottom deflector peak is h_d + _MESH_MIN, and the top
        deflector minimum y is H - (h_d + _MESH_MIN) - delta.
        Gap = H - 2*(h_d + _MESH_MIN) - delta = 2*w_s - delta - 2*_MESH_MIN.
        For Gap >= _MESH_MIN:  2*w_s - delta >= 3*_MESH_MIN.
        Zero or negative → the two deflector solids overlap → self-intersection.
    """
    failures = []

    def _chk(ok: bool, label: str, detail: str) -> None:
        if not ok:
            failures.append(f"{label}: {detail}")

    _chk(L_c >= _MESH_MIN,
         "G1",
         f"interaction length L_c = {L_c:.3e} m  (need >= {_MESH_MIN:.2e} m = 0.01·H)  "
         f"→ reduce L_s+L_m in normalised units to < {(L_cell - _MESH_MIN) / SCALE:.4f}")

    _chk(h_d >= _MESH_MIN,
         "G2",
         f"deflector height h_d = {h_d:.3e} m  (need >= {_MESH_MIN:.2e} m = 0.01·H)  "
         f"→ w_s must be < {(0.5 * H - _MESH_MIN) / SCALE:.4f} normalised  "
         f"(currently w_s = {w_s / SCALE:.4f})")

    _chk(t_s >= _MESH_MIN,
         "G3a",
         f"split splitter t_s = {t_s:.3e} m  (need >= {_MESH_MIN:.2e} m = 0.01·H)")

    _chk(t_m >= _MESH_MIN,
         "G3b",
         f"merge splitter t_m = {t_m:.3e} m  (need >= {_MESH_MIN:.2e} m = 0.01·H)")

    _chk(t_s - t_m >= _MESH_MIN,
         "G4",
         f"t_s - t_m = {(t_s - t_m):.3e} m  (need >= {_MESH_MIN:.2e} m = 0.01·H)  "
         f"→ split splitter must be strictly wider than merge splitter so that "
         f"the split-splitter solid has real material to cut at the cell boundary")

    _chk(w_s - 0.5 * t_s >= _MESH_MIN,
         "G5",
         f"split-section gap w_s - t_s/2 = {(w_s - 0.5*t_s):.3e} m  "
         f"(need >= {_MESH_MIN:.2e} m = 0.01·H)  "
         f"→ fluid channel between deflector peak and splitter surface is too narrow")

    _chk(2 * w_s - delta >= 3 * _MESH_MIN,
         "G6",
         f"deflector gap 2·w_s - delta - 2·_MESH_MIN = {(2*w_s - delta - 2*_MESH_MIN):.3e} m  "
         f"(need >= {_MESH_MIN:.2e} m = 0.01·H;  effective peak = h_d + _MESH_MIN)  "
         f"→ top and bottom deflectors overlap at peak of cosine in interaction region")

    if failures:
        raise ValueError(
            "Geometry validation failed — violated constraints:\n" +
            "\n".join(f"  {f}" for f in failures)
        )


_check_geometry()

# ---------------------------------------------------------------------------
# Helper: cosine-envelope deflector profile as a list of (x, y) points
# ---------------------------------------------------------------------------
N_PTS = 120   # points along the cosine curve

def cosine_bump_points(x_start, x_end, amp, from_top=False, bias=0.0):
    """
    Return a polygon (list of (x,y)) for a cosine-shaped deflector.

    The polygon follows the cosine surface from x_start to x_end, then
    closes back along the wall (y=0 for bottom, y=H for top).

    wall_y is set to -_EPS (bottom) or H+_EPS (top) so the closing edge lies
    just outside the channel boundary.  This prevents the obstacle face from
    being coplanar with the channel-box wall face, which would cause an OCC
    Boolean cut to return a null shape.

    bval has a minimum floor of _MESH_MIN so that the deflector solid always
    protrudes into the channel by at least _MESH_MIN even at the endpoints
    (where the cosine envelope is zero).  Without this floor the solid has
    zero cross-section at xc0/xc1, producing knife-edge geometry in OCC that
    tessellates to zero-area faces → cfMesh creates zero-area mesh faces →
    OpenFOAM deltaCoeffs() raises a floating-point exception at solver start.
    """
    import math
    Lc = x_end - x_start
    wall_y = (H + _EPS) if from_top else -_EPS
    pts = []
    for i in range(N_PTS + 1):
        xi = i * Lc / N_PTS
        env = 0.5 * (1.0 - math.cos(2.0 * math.pi * xi / Lc))
        bval = amp * env + _MESH_MIN   # floor: always protrude >= _MESH_MIN
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
    """Extrude a closed 2-D polygon into a solid slab.

    The workplane is shifted by -2*_EPS in z and the extrusion depth is
    extended by 4*_EPS so the resulting solid spans z = -2·_EPS … depth+2·_EPS.
    The channel box spans z = -_EPS … depth+_EPS (one _EPS inset on each side).
    Obstacle z-faces therefore lie strictly outside the channel-box z-faces,
    avoiding the coplanar-face failure in OpenCASCADE Boolean cuts.
    After the Boolean cut the fluid-domain STL has exactly two distinct
    z-values (-_EPS and depth+_EPS), satisfying cartesian2DMesh's 2D-surface
    uniformity check and eliminating the "z coordinates not uniform" warning.
    """
    wire = (cq.Workplane("XY")
            .workplane(offset=-2 * _EPS)
            .polyline(points_2d)
            .close())
    return wire.extrude(depth + 4 * _EPS)


def rect_points(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


# ---------------------------------------------------------------------------
# Collect obstacle solids per unit cell
# ---------------------------------------------------------------------------
obstacle_solids = []

for k in range(N):
    xk = L0 + k * L_cell

    # --- Split splitter (thin rectangle centred in channel) ---
    xs0, xs1 = xk + DELTA_X, xk + L_s
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
# Channel box spans z = -_EPS … DEPTH+_EPS so its z-faces sit between the
# obstacle z-faces (z = -2·_EPS … DEPTH+2·_EPS).  After Boolean cuts the
# fluid STL has exactly two z-values: -_EPS (front) and DEPTH+_EPS (back).
channel_box = (
    cq.Workplane("XY")
    .box(TOTAL_L, H, DEPTH + 2 * _EPS)
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
