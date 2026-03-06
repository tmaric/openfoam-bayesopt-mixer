# Benchmark-Problem-Specification

Source converted from user-provided mixed overview + technical note into Markdown.

## Problem Overview
Steady, incompressible, laminar 2D channel flow with internal baffles. Two inlet streams with different scalar concentrations mix via advection and diffusion. Baffle geometry is parameterized and optimized with multi-objective Bayesian optimization.

### Objectives
- Minimize scalar mixing defect at outlet.
- Minimize pressure drop across channel.

### Design Vector
- `theta in R^d` maps directly to CAD geometry via Gmsh Python API.

## Computational Domain
- Channel domain: `Omega = {(x,y) | 0 < x < L, 0 < y < H}`.
- Boundary decomposition:
  - `Gamma_in`
  - `Gamma_out`
  - `Gamma_walls`
  - `Gamma_baffles`

## Governing Equations
### Flow
Steady incompressible Navier-Stokes in `Omega`:
- `rho * (u · grad u) = -grad p + mu * laplacian(u)`
- `div(u) = 0`

Variables:
- `u = (u_x, u_y)`: velocity
- `p`: pressure
- `rho`: density
- `mu`: dynamic viscosity

### Scalar Transport
Passive scalar `c(x,y)` in `Omega`:
- `u · grad c = D * laplacian(c)`

`D` is diffusion coefficient.

## Boundary Conditions
### Inlet (`Gamma_in = {0} x (0,H)`)
- `u = (U0, 0)`
- Scalar split:
  - `c = 1` for `y > H/2`
  - `c = 0` for `y <= H/2`

### Outlet (`Gamma_out = {L} x (0,H)`)
- `p = 0`
- `grad(c) · n = 0`

### Walls and Baffles (`Gamma_walls U Gamma_baffles`)
- `u = 0` (no-slip)
- `grad(c) · n = 0`

## Objective Functions
### Mixing Defect
Outlet mean concentration:
- `c_bar = (1 / |Gamma_out|) * integral_{Gamma_out}(c dS)`

For symmetric inlet streams, target is `c_bar = 0.5`.

Mixing defect:
- `J_mix(theta) = (1 / |Gamma_out|) * integral_{Gamma_out} (c - c_bar)^2 dS`

Perfect mixing: `J_mix = 0`.

### Pressure Drop
- `J_delta_p(theta) = <p>_{Gamma_in} - <p>_{Gamma_out}`

### Multi-objective Problem
- `min over theta in D: (J_mix(theta), J_delta_p(theta))`

## Geometry Parameterization
### Channel
- Rectangular base domain `Omega`.

### Baffle Family
- `N` alternating baffles along streamwise direction.
- Positions:
  - `x_k = x_0 + (k - 1) * delta_x`, `k = 1..N`
  - `delta_x = (L - 2*x_0)/(N - 1)`

Parameter vector:
- `theta = (h_1, ..., h_N)`
- `h_min <= h_k <= h_max`

Odd `k` (bottom intrusion):
- `B_k = {(x,y): |x - x_k| <= L_b/2, 0 <= y <= h_k * f_k(x)}`

Even `k` (top intrusion):
- `B_k = {(x,y): |x - x_k| <= L_b/2, H - h_k * f_k(x) <= y <= H}`

Smooth cosine tip shape:
- `f_k(x) = 0.5 * (1 + cos(pi * (x - x_k)/(L_b/2)))` for `|x - x_k| <= L_b/2`
- `f_k(x) = 0` outside support.

Fluid domain:
- `Omega_f(theta) = Omega \ union_{k=1..N} B_k`

## Design Constraints
- `0 < h_k < H - g_min`
- `L_b > 0`
- `x_{k+1} - x_k > L_b`

`g_min` ensures a minimum channel gap for manufacturability and mesh quality.

## Parameter Space
- `D = [h_min, h_max]^N`

Typical setup:
- `N = 6`
- `h_min = 0.1H`
- `h_max = 0.6H`

## Optimization Workflow
1. BO proposes `theta`.
2. Generate `Omega_f(theta)` via Gmsh Python API.
3. Generate mesh.
4. Solve steady flow.
5. Solve scalar transport.
6. Compute objectives.
7. Update BO model.

## Benchmark Qualities
- Fully analytical geometry parameterization.
- Stable steady laminar PDE system.
- Non-convex multi-objective landscape.
- Direct CAD-to-CFD integration.
- Reproducible workflow suitable for BO benchmarking.
