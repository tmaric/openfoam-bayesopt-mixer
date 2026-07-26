#!/usr/bin/env python3
"""
Planar Alternating-Deflector Micromixer - CADQuery geometry script.

Generates the planar mixer fluid-domain boundary (extruded to a thin slab) and
exports it as a single ASCII STL file for use with cfMesh cartesian2DMesh:

  constant/triSurface/alternating_deflector_mixer.stl

The STL contains four named solid regions that cfMesh turns into patches:
  inlet       – face at x = 0
  outlet      – face at x = TOTAL_L
  walls       – channel top/bottom walls + all internal obstacle surfaces
  frontAndBack – slab faces at z = 0 and z = DEPTH  (type empty in 2-D)

The device contains a centre baffle and alternating top/bottom cosine wall
deflectors. It stretches and diffuses the inlet interface but does not route
separate branches out of plane, so it is not described as a true SAR mixer.
"""

import json
import math
import os
import cadquery as cq
import yaml

# ---------------------------------------------------------------------------
# Geometry parameters – loaded from alternating_deflector_cad.yaml
# All values in the YAML are in normalised units (H_norm = 1).
# SCALE converts to SI metres.
# ---------------------------------------------------------------------------
_yaml_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "alternating_deflector_cad.yaml",
)
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
k_slope = _p.get("k", 0.0) * H
TOPOLOGY = _p.get("topology", "alternating_deflector")
if TOPOLOGY not in {"alternating_deflector", "straight"}:
    raise ValueError(
        "topology must be either 'alternating_deflector' or 'straight' "
        f"(received {TOPOLOGY!r})"
    )

# Extrusion depth (thin slab for 2-D OpenFOAM simulation).
# Must satisfy span_z / span_x < ~0.001 for cfMesh cartesian2DMesh to
# classify the surface as 2D.  With TOTAL_L = 24e-3 m the limit is
# DEPTH < 2.4e-5 m.  0.01*H = 1e-5 m → ratio = 4.17e-4 (safe).
DEPTH  = 0.01 * H
# Small boundary offset used in the 2-D profiles and shared z-span so the
# fluid STL remains strictly planar while the y-wall closures avoid exact
# coplanarity with the channel box.
_EPS = DEPTH * 0.01

L_c = L_cell - L_s - L_m
h_d = 0.5 * H - w_s   # deflector intrusion height from each wall

TOTAL_L = 2 * L0 + N * L_cell


def interaction_midpoint_x(cell_idx: int) -> float:
    """Return the global x-location of the cell's deflector midpoint."""
    return L0 + cell_idx * L_cell + L_s + 0.5 * L_c


def interaction_midpoint_xhat(cell_idx: int) -> float:
    """Normalised midpoint used for the piecewise-constant delta profile."""
    return interaction_midpoint_x(cell_idx) / TOTAL_L


CELL_XHATS = [interaction_midpoint_xhat(idx) for idx in range(N)]
CELL_DELTAS = [delta + k_slope * xhat for xhat in CELL_XHATS]

# ---------------------------------------------------------------------------
# Geometry validation  (runs before any OCC operation)
# ---------------------------------------------------------------------------
# Minimum feature size: 1 % of H.  For the default SCALE=1e-3 and H=1.0 this
# is 1e-5 m (10 µm), which equals the wall-refinement cell size in meshDict
# and is consistent with TM_MARGIN=0.01 used in bayes_optimize_sequential.py.
_MESH_MIN = 0.01 * H
_DEFLECTOR_FLOOR = 2.0 * _MESH_MIN

def _check_geometry() -> None:
    """Validate all CAD parameters before any OCC operation is attempted.

    Raises ValueError listing every violated constraint so the caller gets a
    diagnostic message instead of a cryptic 'Null TopoDS_Shape object' error
    from deep inside an OpenCASCADE Boolean operation.

    Each check is labelled Gn and maps to a concrete failure mode:

    G1  L_c >= _MESH_MIN
        Interaction region must have positive x-extent.
        Violated → 'assert L_c > 0' fires, or degenerate cosine polygon.

    G2  h_d >= 0
        The analytic cosine amplitude may be zero because
        cosine_bump_points() adds a mesh-resolved floor. Negative amplitudes
        are invalid.

    G3a t_s >= _MESH_MIN
    G3b t_m >= _MESH_MIN
        Both splitters must be thick enough for at least one mesh cell.

    G4  t_s - t_m >= _MESH_MIN
        Split splitter must be strictly wider than the merge splitter so
        that the connected merge→split polygon has a real thickness step
        at each cell boundary.

    G5  w_s - 0.5*t_s >= _MESH_MIN
        Minimum fluid gap in the split section between the deflector peak
        and the nearest splitter surface.
        Bottom deflector peak y = h_d = H/2 - w_s;
        split-splitter bottom face y = (H - t_s)/2.
        Gap = w_s - t_s/2.  Near-zero → channel pinches → meshing fails.

    G6a min(delta_i) >= 0
        The alternating wall-bias amplitude is interpreted as an inward shift
        magnitude. Negative realised values would reverse that meaning.

    G6b 2*w_s - max(delta_i) >= 5*_MESH_MIN
        Minimum fluid gap between top and bottom deflectors at the peak of
        the cosine in the interaction region.
        The cosine bval has a floor of 2*_MESH_MIN (see cosine_bump_points), so
        the effective peak intrusion is h_d + 2*_MESH_MIN + delta_i on the
        shifted wall. The worst case therefore depends on max(delta_i).
        Gap = 2*w_s - max(delta_i) - 4*_MESH_MIN.
        For Gap >= _MESH_MIN:  2*w_s - max(delta_i) >= 5*_MESH_MIN.
        Zero or negative → the two deflector solids overlap → self-intersection.
    """
    failures = []
    min_delta = min(CELL_DELTAS) if CELL_DELTAS else delta
    max_delta = max(CELL_DELTAS) if CELL_DELTAS else delta

    def _chk(ok: bool, label: str, detail: str) -> None:
        if not ok:
            failures.append(f"{label}: {detail}")

    _chk(N >= 1,
         "G0",
         "number of unit cells N must be >= 1")

    _chk(L_c >= _MESH_MIN,
         "G1",
         f"interaction length L_c = {L_c:.3e} m  (need >= {_MESH_MIN:.2e} m = 0.01·H)  "
         f"→ reduce L_s+L_m in normalised units to < {(L_cell - _MESH_MIN) / SCALE:.4f}")

    _chk(h_d >= 0.0,
         "G2",
         f"deflector height h_d = {h_d:.3e} m  (need >= 0)  "
         f"→ w_s must be <= {0.5 * H / SCALE:.4f} normalised  "
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

    _chk(min_delta >= 0.0,
         "G6a",
         f"minimum realised deflector bias min(delta_i) = {min_delta:.3e} m  "
         f"(need >= 0)  → linear slope k drives at least one cell to a negative wall-bias amplitude")

    _chk(2 * w_s - max_delta >= 5 * _MESH_MIN,
         "G6b",
         f"deflector gap 2·w_s - max(delta_i) - 4·_MESH_MIN = {(2*w_s - max_delta - 4*_MESH_MIN):.3e} m  "
         f"(need >= {_MESH_MIN:.2e} m = 0.01·H;  effective peak = h_d + 2·_MESH_MIN)  "
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
# Keep the STL resolution commensurate with the 12.5 um finest CFD cells.
# The former 120-segment profile produced nanometre-scale y increments near
# the zero-slope cosine endpoints, which cfMesh converted into short edges and
# highly skew/concave cells. Forty-eight segments still resolve the shortest
# admissible 0.4 mm interaction region at roughly the CFD cell spacing.
N_PTS = 48

def cosine_bump_points(x_start, x_end, amp, from_top=False, bias=0.0):
    """
    Return a polygon (list of (x,y)) for a cosine-shaped deflector.

    The polygon follows the cosine surface from x_start to x_end, then
    closes back along the wall (y=0 for bottom, y=H for top).

    wall_y is set to -_EPS (bottom) or H+_EPS (top) so the closing edge lies
    just outside the channel boundary.  This prevents the obstacle face from
    being coplanar with the channel-box wall face, which would cause an OCC
    Boolean cut to return a null shape.

    bval has a minimum floor of 2*_MESH_MIN so that the deflector solid always
    protrudes into the channel by more than one fine cell at the endpoints
    (where the cosine envelope is zero).  Without this floor the solid has
    zero cross-section at xc0/xc1, producing knife-edge geometry in OCC that
    tessellates to zero-area faces → cfMesh creates zero-area mesh faces →
    OpenFOAM deltaCoeffs() raises a floating-point exception at solver start.
    """
    Lc = x_end - x_start
    wall_y = (H + _EPS) if from_top else -_EPS
    pts = []
    for i in range(N_PTS + 1):
        xi = i * Lc / N_PTS
        env = 0.5 * (1.0 - math.cos(2.0 * math.pi * xi / Lc))
        bval = amp * env + _DEFLECTOR_FLOOR
        x = x_start + xi
        if from_top:
            y = H - min(bval + bias, H)
        else:
            y = min(bval + bias, H)
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

    The obstacle slab spans the exact same z-range as the channel box:
    z = -_EPS … depth+_EPS. This keeps the exported surface strictly 2-D,
    with exactly two z-levels, which cartesian2DMesh can classify without
    generating sliver triangles at the front/back planes.
    """
    wire = (cq.Workplane("XY")
            .workplane(offset=-_EPS)
            .polyline(points_2d)
            .close())
    return wire.extrude(depth + 2 * _EPS)


def rect_points(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def centered_band_points(x0, x1, thickness):
    y0 = (H - thickness) / 2.0
    y1 = (H + thickness) / 2.0
    return rect_points(x0, x1, y0, y1)


def centered_step_points(x0, x_step, x1, left_thickness, right_thickness):
    """Connected polygon for a merge-to-next-split center plate transition.

    The thickness transition is spread over a small x-ramp instead of a
    zero-width vertical step. cfMesh repeatedly left inverted/zero-area
    faces at the exact step corners, while a one-mesh-cell taper retains
    topological connectivity and meshes robustly.
    """
    y0l = (H - left_thickness) / 2.0
    y1l = (H + left_thickness) / 2.0
    y0r = (H - right_thickness) / 2.0
    y1r = (H + right_thickness) / 2.0
    # Resolve the thickness transition over several fine cells.  A one-cell
    # taper repeatedly intersected the Cartesian grid in nanometre-scale edge
    # fragments and produced two concave cells at every unit-cell boundary.
    ramp_dx = min(
        4.0 * _MESH_MIN,
        0.25 * max(x_step - x0, 0.0),
        0.25 * max(x1 - x_step, 0.0),
    )
    xl = x_step - 0.5 * ramp_dx
    xr = x_step + 0.5 * ramp_dx
    return [
        (x0, y0l),
        (xl, y0l),
        (xr, y0r),
        (x1, y0r),
        (x1, y1r),
        (xr, y1r),
        (xl, y1l),
        (x0, y1l),
    ]


# ---------------------------------------------------------------------------
# Collect obstacle solids per unit cell
# ---------------------------------------------------------------------------
obstacle_solids = []

# The center plate is built from connected polygons so there is never a
# topological hole between the merge plate of cell i and the split plate of
# cell i+1. The direct thickness step lives inside one polygon, not as two
# adjacent rectangles with an x-gap.
if TOPOLOGY == "alternating_deflector":
    first_split_x0 = L0
    first_split_x1 = L0 + L_s
    obstacle_solids.append(
        make_extruded_polygon(centered_band_points(first_split_x0, first_split_x1, t_s), DEPTH)
    )

    for idx in range(N - 1):
        boundary_x = L0 + (idx + 1) * L_cell
        obstacle_solids.append(
            make_extruded_polygon(
                centered_step_points(
                    boundary_x - L_m,
                    boundary_x,
                    boundary_x + L_s,
                    t_m,
                    t_s,
                ),
                DEPTH,
            )
        )

    last_merge_x0 = L0 + (N - 1) * L_cell + L_s + L_c
    last_merge_x1 = L0 + N * L_cell
    obstacle_solids.append(
        make_extruded_polygon(centered_band_points(last_merge_x0, last_merge_x1, t_m), DEPTH)
    )

    for idx in range(N):
        xk = L0 + idx * L_cell
        xc0, xc1 = xk + L_s, xk + L_s + L_c
        delta_i = CELL_DELTAS[idx]
        bot_bias = delta_i if idx % 2 == 1 else 0.0
        top_bias = delta_i if idx % 2 == 0 else 0.0

        pts_bot = cosine_bump_points(xc0, xc1, h_d, from_top=False, bias=bot_bias)
        obstacle_solids.append(make_extruded_polygon(pts_bot, DEPTH))

        pts_top = cosine_bump_points(xc0, xc1, h_d, from_top=True, bias=top_bias)
        obstacle_solids.append(make_extruded_polygon(pts_top, DEPTH))


# ---------------------------------------------------------------------------
# Channel outer box and fluid domain
# ---------------------------------------------------------------------------
# Channel box and obstacle slabs share the same z-span: z = -_EPS … DEPTH+_EPS.
# The exported surface therefore has exactly two z-values, which lets cfMesh
# detect it as a true 2-D x-y surface.
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
# Face classification by physical location
# ---------------------------------------------------------------------------
_NORMAL_TOL = 1e-6
_POSITION_TOL = max(1e-12, TOTAL_L * 1e-9)


def _face_is_on_x_plane(face, x_value):
    """Return True only when the complete face lies on the requested x-plane.

    Normal-based classification is unsafe here: every upstream-facing
    obstacle face also has a negative x normal, and every downstream-facing
    obstacle face has a positive x normal.  Requiring both x bounds to lie on
    the exterior plane makes it impossible for an interior deflector surface
    to become an inlet or outlet.
    """
    bounds = face.BoundingBox()
    return (
        abs(bounds.xmin - x_value) <= _POSITION_TOL
        and abs(bounds.xmax - x_value) <= _POSITION_TOL
    )


def _is_inlet(face):
    return _face_is_on_x_plane(face, 0.0)


def _is_outlet(face):
    return _face_is_on_x_plane(face, TOTAL_L)


def _is_front_or_back(face):
    normal = face.normalAt()
    return abs(normal.z) > (1.0 - _NORMAL_TOL)


def _is_wall(f):
    return not _is_inlet(f) and not _is_outlet(f) and not _is_front_or_back(f)

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


def _combined_bounds(faces):
    boxes = [face.BoundingBox() for face in faces]
    return {
        "xmin": min(box.xmin for box in boxes),
        "xmax": max(box.xmax for box in boxes),
        "ymin": min(box.ymin for box in boxes),
        "ymax": max(box.ymax for box in boxes),
        "zmin": min(box.zmin for box in boxes),
        "zmax": max(box.zmax for box in boxes),
    }


def _validate_patch_partition():
    """Fail CAD generation if physical boundary ownership is inconsistent."""
    inlet_faces = patch_faces["inlet"]
    outlet_faces = patch_faces["outlet"]
    wall_faces = patch_faces["walls"]
    front_back_faces = [face for face in all_faces if _is_front_or_back(face)]

    if not inlet_faces or not outlet_faces or not wall_faces or not front_back_faces:
        raise RuntimeError(
            "CAD patch classification produced an empty physical boundary: "
            f"inlet={len(inlet_faces)}, outlet={len(outlet_faces)}, "
            f"walls={len(wall_faces)}, front/back={len(front_back_faces)}"
        )

    groups = (inlet_faces, outlet_faces, wall_faces, front_back_faces)
    identities = [id(face) for group in groups for face in group]
    if len(identities) != len(set(identities)) or len(identities) != len(all_faces):
        raise RuntimeError(
            "CAD patch classification must assign every face to exactly one boundary"
        )

    inlet_bounds = _combined_bounds(inlet_faces)
    outlet_bounds = _combined_bounds(outlet_faces)
    if (
        abs(inlet_bounds["xmin"]) > _POSITION_TOL
        or abs(inlet_bounds["xmax"]) > _POSITION_TOL
    ):
        raise RuntimeError(f"inlet is not confined to x=0: {inlet_bounds}")
    if (
        abs(outlet_bounds["xmin"] - TOTAL_L) > _POSITION_TOL
        or abs(outlet_bounds["xmax"] - TOTAL_L) > _POSITION_TOL
    ):
        raise RuntimeError(f"outlet is not confined to x=L: {outlet_bounds}")

    expected_area = H * (DEPTH + 2.0 * _EPS)
    inlet_area = sum(face.Area() for face in inlet_faces)
    outlet_area = sum(face.Area() for face in outlet_faces)
    for name, area in (("inlet", inlet_area), ("outlet", outlet_area)):
        if not math.isclose(area, expected_area, rel_tol=1e-7, abs_tol=1e-16):
            raise RuntimeError(
                f"{name} area {area:.12e} m^2 differs from expected "
                f"H*depth={expected_area:.12e} m^2"
            )

    return {
        "schema_version": 1,
        "topology": TOPOLOGY,
        "total_length_m": TOTAL_L,
        "channel_height_m": H,
        "extrusion_depth_m": DEPTH + 2.0 * _EPS,
        "position_tolerance_m": _POSITION_TOL,
        "patches": {
            "inlet": {
                "cad_face_count": len(inlet_faces),
                "area_m2": inlet_area,
                "bounds_m": inlet_bounds,
            },
            "outlet": {
                "cad_face_count": len(outlet_faces),
                "area_m2": outlet_area,
                "bounds_m": outlet_bounds,
            },
            "walls": {
                "cad_face_count": len(wall_faces),
                "area_m2": sum(face.Area() for face in wall_faces),
                "bounds_m": _combined_bounds(wall_faces),
            },
        },
    }


geometry_manifest = _validate_patch_partition()

print("Exporting patch STL solids ...")
stl_parts = []
for name, faces in patch_faces.items():
    print(f"  Patch '{name}': {len(faces)} face(s)")
    stl_parts.append(_stl_block(faces, name))

stl_path = os.path.join(out_dir, "alternating_deflector_mixer.stl")
with open(stl_path, "w") as fh:
    fh.write("".join(stl_parts))
print(f"  Written: {stl_path}")

manifest_path = os.path.join(out_dir, "geometry_manifest.json")
with open(manifest_path, "w") as fh:
    json.dump(geometry_manifest, fh, indent=2, sort_keys=True)
    fh.write("\n")
print(f"  Written: {manifest_path}")

print()
print("Geometry summary:")
print(f"  Topology           = {TOPOLOGY}")
print(f"  Total length  L    = {TOTAL_L:.4f}")
print(f"  Channel height H   = {H:.4f}")
print(f"  Extrusion depth    = {DEPTH:.4f}")
print(f"  Unit cells    N    = {N}")
print(f"  Unit-cell len      = {L_cell:.4f}")
print(f"  Splitter thick t_s = {t_s:.4f}")
print(f"  Splitter thick t_m = {t_m:.4f}")
print(f"  Deflector height h_d = {h_d:.4f}")
print(f"  Base deflector bias delta = {delta:.6e}")
print(f"  Linear delta slope k = {k_slope:.6e}")
print(f"  Realised delta_i range = [{min(CELL_DELTAS):.6e}, {max(CELL_DELTAS):.6e}]")
print("  Wall bias pattern = top, bottom, top, ...")
print(f"  Interaction length L_c = {L_c:.4f}")
print("Done.")
