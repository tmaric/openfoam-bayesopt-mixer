# Six-unit review protocol audit

Status on 2026-07-28: the review metric and axial-station ambiguity have been
audited. Neither explains the fine-grid Re=1 discrepancy. The remaining issue
is the strong, monotone dependence of mixing on resolution, which points to
the reconstructed inter-layer connection topology or its mesh representation.
The numerical table in this document predates the production switch from
cfMesh to OpenFOAM-v2606 snappyHexMesh and remains historical comparison
evidence; the snappy mesh-only six-design gate is documented separately in
`mesh_qualification.md`.

## Source definition

The source is Raza, Hossain, and Kim, *Micromachines* 11 (2020) 455,
DOI `10.3390/mi11050455`. Its comparison protocol states:

- six M10 mixing units with the original unit dimensions;
- an exit 5050 um downstream of the start of the mixing unit;
- water properties at 25 C, including `mu=8.8e-4 Pa s` and
  `rho=997 kg/m3`;
- scalar diffusivity `1e-10 m2/s`;
- Reynolds number based on inlet mean velocity and inlet hydraulic diameter;
- uniform inlet velocities, zero outlet static pressure, and no-slip walls;
- M10 targets `M=0.915`, `DeltaP=16.3 Pa` at Re=1 and `M=0.901`,
  `DeltaP=390 Pa` at Re=20.

The source geometry figure gives `H=1070 um`, `P=640 um`, `b=150 um`,
`d=150 um`, and `w=300 um`. The implemented dimensions, topology, two
`0.15 x 0.30 mm` inlets, `0.30 x 0.30 mm` outlet, and six-unit count agree
with that schematic. The exact distribution of the small residual lead length
is not dimensioned in the source.

## Inter-layer aperture diagnostic

The CAD manifest now measures the projected open interface without changing
the geometry. Because all side walls are normal to the layer interface, the
area is obtained exactly by translating the upper solid through a thin probe
depth, intersecting it with the lower solid, and dividing overlap volume by
that depth.

For the six-unit reference it reports:

| Connection region | Open area |
|---|---:|
| Seven vertical segments | 1.12350 mm2 |
| Six X crossings outside the vertical segments | 0.557822 mm2 |
| Total mixing-core interface | 1.681322 mm2 |

The total is 37.36 times one inlet cross-sectional area. These numbers are the
full projected overlaps implied by the source schematic and are now covered by
a regression test. The publication does not provide an independent mask or
aperture dimension, so whether its numerical model opened the complete overlap
or a restricted subset cannot be proven from the review figure alone. This is
the most plausible remaining geometric source of the fine-grid mixing excess.

## Mixing-index equivalence

The review samples concentration at equally weighted points on a transverse
exit plane and defines

```text
sigma = sqrt(mean((c - 0.5)^2))
M = 1 - sigma/sigma_max
```

For equal inlet fractions, `sigma_max=0.5`. The function object's
area-weighted intensity is

```text
I_A = mean_A((T - 0.5)^2) / 0.25
```

and therefore its reported literature quantity `1-sqrt(I_A)` is exactly
`1-sigma/0.5` in the continuum limit. The historical accepted cfMesh outlet
is a uniform Cartesian face grid, so equal-point and area weighting are
equivalent to the accuracy relevant here. Flux weighting remains a separate
BO objective and is not used for the literature gate.

## Axial-station sensitivity

The six-unit CAD core is 4.89 mm. Symmetric 0.08 mm inlet and outlet leads give
the standardized 5.05 mm total axial extent and are the shortest clean leads:
an attempted 0.02 mm inlet intersected the first diagonal channel at the inlet
plane and correctly failed the two-face inlet topology check.

Two valid coarse cases added 80 um on either side to cover the plausible origin
of the source's axial coordinate. All use the production second-order scalar
scheme and four MPI ranks.

| Inlet lead (mm) | Outlet lead (mm) | Pressure (Pa) | Area MI | Change from symmetric |
|---:|---:|---:|---:|---:|
| 0.08 | 0.08 | 14.224 | 0.913232 | reference |
| 0.08 | 0.16 | 14.305 | 0.912678 | -0.000554 |
| 0.16 | 0.08 | 14.458 | 0.909899 | -0.003333 |

These changes are an order of magnitude smaller than the fine-grid target
error `+0.04133`. No fine lead-sensitivity run is justified.

## Resolution trend at Re=1

The corrected scheme produces a monotone trend across three resolutions:

| Nominal cell size | Pressure (Pa) | Area MI | MI error from 0.915 | Flux MI |
|---:|---:|---:|---:|---:|
| 30 um | 14.224 | 0.913232 | -0.001768 | 0.923576 |
| 20 um | 14.903 | 0.942197 | +0.027197 | 0.945942 |
| 14 um | 15.258 | 0.956333 | +0.041333 | 0.957138 |

Pressure moves toward the `16.3 Pa` target with refinement, but mixing moves
away from `0.915`. Agreement on the 30 um mesh is therefore accidental and
must not be used to clear the BO gate. This is not ordinary excessive
coarse-grid numerical diffusion: resolving the crossing topology increases
mixing. The next geometry audit must focus on the areas and shapes of the
openings between layers at the X nodes and vertical segments.

## Straight-channel pressure normalization

Every completed run now records a design-matched ideal straight reference with
the same total axial length, outlet cross-section `w x D`, mean velocity, and
fluid viscosity. For longer side `a` and shorter side `b`, the fully developed
rectangular-duct series is

```text
DeltaP_0 = 12 mu L U / (b^2 K)
K = 1 - 192 b/(pi^5 a) * sum(n odd, tanh(n pi a/(2b))/n^5)
```

Fifty odd terms are used. The implementation reproduces the square-duct Darcy
Poiseuille number `f Re = 56.91`. At the four-unit reference geometry and
Re=10, `DeltaP_0=43.5885 Pa`. For the six-unit review geometry, the ideal
values are `6.20061 Pa` at Re=1 and `124.012 Pa` at Re=20. This normalization
is evaluated per design, avoiding an extra CFD solve for every BO observation.

## Gate decision

The metric and lead-location questions are closed, and the straight pressure
objective is operational. BO remains gated because the target-fidelity scalar
result is not mesh-converged toward the review result. Before initialization:

1. verify from the original model, fabrication mask, or author data whether the
   complete projected interface areas above were open in the CEJ simulations;
2. decide explicitly whether the study claims an M10 reproduction or an
   M10-inspired reconstruction;
3. if it is a reconstruction, revise the literature gate and claims before a
   small paired coarse/fine rank-correlation pilot;
4. launch the 24+6 initialization only after that pilot supports a useful
   multifidelity correlation for both objectives.
