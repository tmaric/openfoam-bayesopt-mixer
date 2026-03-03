# Technical Note 01: Passive Mixer Benchmark (2D Laminar Two-Stream)

## Scope
This note captures the benchmark definition for the passive mixer case and includes a geometry/parameterization sketch used during implementation.

## Resources
- Python sketch script: [../resources/passive-mixer/passive_mixer_parameterization_sketch.py](../resources/passive-mixer/passive_mixer_parameterization_sketch.py)
- Rendered sketch PNG: [../resources/passive-mixer/passive_mixer_parameterization_sketch.png](../resources/passive-mixer/passive_mixer_parameterization_sketch.png)

![Passive mixer parameterization sketch](../resources/passive-mixer/passive_mixer_parameterization_sketch.png)

## Problem Overview
We consider a steady, incompressible, laminar two-dimensional flow in a channel with internal baffles. Two fluid streams with different scalar concentrations enter the channel and mix by advection and diffusion. The geometry of internal baffles is parameterized and optimized with multi-objective Bayesian optimization.

### Objectives
- Minimize scalar mixing defect at the outlet.
- Minimize pressure drop across the channel.

The geometry is analytically defined by a finite-dimensional parameter vector `theta in R^d`, mapped directly to CAD via the Gmsh Python API.

## Computational Domain
- Channel domain: `Omega = {(x,y) in R^2 | 0 < x < L, 0 < y < H}`.
- Boundary decomposition: `partial Omega = Gamma_in U Gamma_out U Gamma_walls U Gamma_baffles`.

## Governing Equations
### Flow
Steady incompressible Navier-Stokes in `Omega`:
- `rho (u · grad)u = -grad p + mu laplacian(u)`
- `div(u) = 0`

with velocity `u = (u_x, u_y)`, pressure `p`, density `rho`, and dynamic viscosity `mu`.

### Scalar Transport
Passive scalar `c(x,y)` in `Omega`:
- `u · grad c = D laplacian(c)`

with diffusion coefficient `D`.

## Boundary Conditions
### Inlet (`Gamma_in = {0} x (0,H)`)
- Velocity: `u = (U0, 0)`
- Scalar:
  - `c = 1` for `y > H/2`
  - `c = 0` for `y <= H/2`

### Outlet (`Gamma_out = {L} x (0,H)`)
- `p = 0`
- `grad(c) · n = 0`

### Walls and Baffles (`Gamma_walls U Gamma_baffles`)
- `u = 0`
- `grad(c) · n = 0`

## Objective Functions
### Mixing Defect
Outlet mean concentration:
- `c_bar = (1/|Gamma_out|) * integral_{Gamma_out} c dS`

For symmetric inlet streams, `c_bar = 0.5`.

Mixing defect:
- `J_mix(theta) = (1/|Gamma_out|) * integral_{Gamma_out} (c - c_bar)^2 dS`

Perfect mixing gives `J_mix = 0`.

### Pressure Drop
- `J_delta_p(theta) = <p>_{Gamma_in} - <p>_{Gamma_out}`

### Multi-objective Formulation
- `min_{theta in D} (J_mix(theta), J_delta_p(theta))`

## Geometry Parameterization
### Channel
Rectangular base domain `Omega`.

### Baffle Family
Let `N` alternating baffles be placed along the channel.

- Streamwise positions:
  - `x_k = x_0 + (k-1)Delta x`, `k = 1..N`
  - `Delta x = (L - 2x_0)/(N-1)`

- Parameter vector:
  - `theta = (h_1, h_2, ..., h_N)`
  - `h_min <= h_k <= h_max`

Odd `k` (bottom intrusion):
- `B_k = {(x,y): |x-x_k| <= L_b/2, 0 <= y <= h_k f_k(x)}`

Even `k` (top intrusion):
- `B_k = {(x,y): |x-x_k| <= L_b/2, H-h_k f_k(x) <= y <= H}`

Shape function (smooth cosine tip):
- `f_k(x) = 0.5 * [1 + cos(pi (x-x_k)/(L_b/2))]` for `|x-x_k| <= L_b/2`
- `f_k(x)=0` otherwise

Final fluid domain:
- `Omega_f(theta) = Omega \ union_{k=1}^N B_k`

## Design Constraints
- `0 < h_k < H - g_min`
- `L_b > 0`
- `x_{k+1} - x_k > L_b`

`g_min` enforces minimum channel gap for manufacturability and mesh quality.

## Parameter Space
- `D = [h_min, h_max]^N`

Typical example:
- `N = 6`
- `h_min = 0.1H`
- `h_max = 0.6H`

## Optimization Workflow
1. BO proposes `theta`.
2. Geometry `Omega_f(theta)` generated via Gmsh Python API.
3. Mesh generated.
4. Steady flow solved.
5. Scalar transport solved.
6. Objectives computed.
7. BO model updated.

## Summary
This benchmark provides:
- Analytical geometry parameterization.
- Stable steady laminar PDE system.
- Non-convex multi-objective landscape.
- Direct CAD-to-CFD integration.
- Reproducible BO benchmarking workflow.
