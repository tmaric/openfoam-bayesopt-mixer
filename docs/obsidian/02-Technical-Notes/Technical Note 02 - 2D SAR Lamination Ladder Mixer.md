# Technical Note 02: 2D Split-and-Recombine (SAR) Lamination Ladder Mixer

## Motivation and Design Principle
Split-and-recombine (SAR) micromixers improve mixing by repeatedly splitting inlet streams into substreams and recombining them to laminate interfaces into thinner layers, which accelerates diffusion.

This note defines a planar (2D) SAR lamination ladder mixer family that:
- is straightforward to parameterize in Gmsh Python API,
- gives effective lamination in steady laminar flow,
- yields a non-convex design space for multi-objective Bayesian optimization (BO).

## Physical Model
Steady incompressible Navier-Stokes in fluid domain `Omega_f(theta) subset R^2`:
- `rho (u · grad)u = -grad p + mu laplacian(u)`
- `div(u) = 0`

Passive scalar advection-diffusion:
- `u · grad c = D laplacian(c)`

Dimensionless groups:
- `Re = rho U0 H / mu`
- `Pe = U0 H / D`

## Boundary Conditions
Let:
- `Gamma_in = {0} x (0,H)`
- `Gamma_out = {L} x (0,H)`
- `Gamma_w` = all solid walls (outer + internal)

Inlet:
- `u = (U0, 0)`
- `c(0,y)=1` for `y > H/2`, `c(0,y)=0` for `y <= H/2`

Outlet:
- `p = 0`
- `grad(c) · n = 0`

Walls:
- `u = 0`
- `grad(c) · n = 0`

## Objectives (Multi-objective BO)
Mixing defect at outlet:
- `c_bar = (1 / |Gamma_out|) * integral_{Gamma_out} c dS`
- `J_mix(theta) = (1 / |Gamma_out|) * integral_{Gamma_out} (c - c_bar)^2 dS`

For symmetric inlet streams, `c_bar ~ 0.5`.

Pressure drop:
- `J_delta_p(theta) = <p>_{Gamma_in} - <p>_{Gamma_out}`

Optimization problem:
- `min_{theta in D} (J_mix(theta), J_delta_p(theta))`

## Geometry: SAR Lamination Ladder (AI-agent-readable CAD Parameterization)

### Reference Dimensions
Microfluidic-scale reference setup:
- `H = 200 um`
- `L_cell = 800 um`
- `L = 2L0 + N L_cell`
- `L0 = 400 um`
- `N` unit cells (example: `N=6`)

### High-level Concept
Each unit cell performs:
- split,
- shuffle,
- recombine.

Repeated cells increase lamination count and reduce striation thickness.

### Design Region and Cell Placement
Full channel:
- `Omega = (0,L) x (0,H)`

Cell `k` (`k = 1..N`):
- `I_k = [x_k, x_k + L_cell]`
- `x_k = L0 + (k-1)L_cell`

### Parameter Vector
- `theta = (w_s, t_s, L_s, L_m, delta, r)`

Where:
- `w_s`: nominal subchannel width after splitting
- `t_s`: splitter thickness
- `L_s`: split length
- `L_m`: merge length
- `delta`: shuffle offset controlling recombination misalignment
- `r`: fillet radius at sharp corners

Remaining interaction length:
- `L_c = L_cell - L_s - L_m`, with `L_c > 0`

### Constraints (CAD Validity)
- `0 < t_s <= 0.15H`
- `0.25H <= w_s <= 0.45H`
- `0.2L_cell <= L_s <= 0.45L_cell`
- `0.2L_cell <= L_m <= 0.45L_cell`
- `L_c = L_cell - L_s - L_m > 0`
- `|delta| <= 0.15H`
- `0 <= r <= 0.05H`

### Obstacle Construction per Cell
Define splitter centerline at `y = H/2`.

1. Splitter wall (`S_k`) on split section:
- `x in [x_k, x_k + L_s]`
- `y in [(H - t_s)/2, (H + t_s)/2]`

2. Shuffle deflectors on `[x_k + L_s, x_k + L_s + L_c]`:
- `h_d = H/2 - w_s`, require `h_d > 0`
- envelope `g(xi) = 0.5 * (1 - cos(2pi xi / L_c))`, `xi in [0, L_c]`
- `xi = x - (x_k + L_s)`

Bottom deflector `D_k_bot`:
- `x in [x_k + L_s, x_k + L_s + L_c]`
- `y in [0, h_d g(xi)]`

Top deflector `D_k_top`:
- `x in [x_k + L_s, x_k + L_s + L_c]`
- `y in [H - h_d g(xi), H]`
- apply vertical offset `delta` in CAD control points before surface creation and clipping.

3. Merge splitter (`M_k`) on merge section:
- `x in [x_k + L_s + L_c, x_k + L_cell]`
- `y in [(H - t_m)/2, (H + t_m)/2]`
- recommended `t_m = 0.05H` (or `t_m = t_s`)

### Final Fluid Domain
Total solid set:
- `B(theta) = union_{k=1..N} (S_k U D_k_bot U D_k_top(delta) U M_k)`

Fluid domain:
- `Omega_f(theta) = Omega \ B(theta)`

### Filleting for Mesh Quality
Apply fillet radius `r` at obstacle corners and wall junctions (Gmsh OCC fillet) to reduce sliver cells and model manufacturability constraints.

## Recommended Operating Point and Robustness Check
Practical benchmark range:
- `Re in [1, 50]`
- `Pe in [10^2, 10^4]`

Suggested BO strategy:
- optimize at nominal `(Re0, Pe0)`
- validate Pareto set on a small `(Re, Pe)` grid.

## Remarks
This 2D SAR lamination ladder provides:
- robust low-dimensional CAD parameterization,
- clear trade-off between mixing and pressure drop,
- geometry aligned with SAR lamination principles.

## References
- Kim et al. (2005), serpentine laminating micromixer with split/recombine and advection.
- Hossain et al. (2015), 3D serpentine split-and-recombine mixing analysis.
- Taheri et al. (2019), SAR micromixer mixing performance.
- Juraeva et al. (2020), cross-channel SAR micromixer.
- Raza et al. (2019), asymmetrical SAR with baffles.
- Nishu et al. (2023), SAR micromixer with dislocated connecting channels.
