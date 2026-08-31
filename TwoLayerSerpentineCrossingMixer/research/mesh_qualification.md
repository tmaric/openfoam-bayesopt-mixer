# snappyHexMesh qualification (OpenFOAM-v2606)

The production mesh pipeline is now

```text
blockMesh -> surfaceFeatureExtract -> snappyHexMesh -> createPatch -> checkMesh
```

It uses the OpenFOAM-v2606 `snappyHexMesh` implementation on a uniform
hexahedral background, with no boundary layers and no surface-refinement level
change. Meshing uses at most four MPI ranks and cases are still evaluated one
at a time. The reconstructed mesh is hexahedral-dominant: oblique fitted walls
necessarily produce a thin prism/polyhedral cut-cell population, but actual
tetrahedra are forbidden.

The two STL inlet regions are merged into the solver-facing `inlet` patch by
`createPatch`; `outlet` and `walls` remain separate. The runner rejects a
non-empty leaked background or STL region patch.

## Nine-unit reference

| Nominal size | Cells | Strict hex fraction | Prisms | Tet wedges | Tetrahedra | Polyhedra | Decision |
|---:|---:|---:|---:|---:|---:|---:|---|
| 24 um | 84,528 | 93.562% | 4,991 | 0 | 0 | 451 | coarse accepted |
| 13 um | 579,096 | 96.352% | 19,765 | 1 | 0 | 1,361 | fine accepted |

The four-rank 24 um mesh has 38 concave cells, 19 warped faces, maximum
non-orthogonality `54.98 deg`, and maximum skewness `1.998`. The four-rank
13 um mesh has 117 concave cells, 16 warped faces, maximum non-orthogonality
`55.00 deg`, and maximum skewness `2.000`. Every explicit face-quality count
is zero at both resolutions.

## Six-design BO preflight

The fixed six paired Sobol anchors were meshed at both fidelities before any
new transport evaluation. The 12 meshes pass the same policy:

| Fidelity | Meshes | Cell range | Strict hex range | Tetrahedra | Max concave-cell fraction | Max warped-face fraction | Min face flatness |
|---|---:|---:|---:|---:|---:|---:|---:|
| coarse, 24 um | 6 | 32,915--49,581 | 91.972--98.987% | 0 | 0.07795% | 0.02110% | 0.5011 |
| fine, 13 um | 6 | 220,460--340,732 | 95.503--99.335% | 0 | 0.00729% | 0.00498% | 0.4359 |

All 12 have zero explicit errors for non-orthogonality, pyramid volume,
face-decomposition tet quality, concavity above the independent limit,
skewness, interpolation weight, neighbour-volume ratio, face twist, and cell
determinant. The ignored runtime evidence is under
`results/mesh_qualification/snappy_v2606_pilot_preflight_24um_13um/`.

## Complete solver smoke test

One four-rank nine-unit Re=20 case reused and independently rechecked the
accepted 24 um reference mesh, then completed `simpleFoam`, reconstruction,
second-order `scalarTransportFoam`, both repository-local function objects,
and every post-run gate under OpenFOAM-v2606:

| Pressure | Pressure ratio | Flux MI | Area MI | Mass-balance error | Scalar bounds |
|---:|---:|---:|---:|---:|---:|
| 493.282 Pa | 2.67831 | 0.990867 | 0.992376 | 2.96e-9 | [-2.38e-7, 1.0000007] |

The final-50 flux-segregation span is `2e-10` and the concentration mean is
stationary. This is an end-to-end compatibility test, not a fine-grid
literature comparison or mesh-convergence result.

## Fixed acceptance policy

`FlowCase/system/snappyMeshQualityDict` constrains mesh construction to:

- non-orthogonality at most 55 degrees;
- internal and boundary skewness at most 2;
- snapping concavity at most 30 degrees;
- determinant at least 0.01;
- interpolation weight at least 0.05;
- neighbouring-cell volume ratio at least 0.02.

`FlowCase/system/meshQualityDict` independently checks the reconstructed mesh.
Its concavity limit is 35 degrees because `snappyHexMesh` and `checkMesh`
evaluate fitted cut-face concavity differently; OpenFOAM's default is 80
degrees. `research_config.yaml` additionally requires:

- exactly zero tetrahedra at both fidelities;
- at least 90% strict hexahedra coarse and 95% fine;
- bounded tet-wedge, concave-cell, and warped-face populations;
- minimum face flatness of 0.30;
- no more than the single aggregate default concave/warped cut-cell warning.

The warped-face caps were fixed after the six-design qualification at
`2.5e-4` coarse and `6e-5` fine, leaving modest margins above observed maxima
of `2.1105e-4` and `4.9769e-5`. This does not relax any explicit
operator-oriented error count or permit tetrahedra.

## Superseded cfMesh evidence

The earlier source/review calculations used cfMesh `cartesianMesh`. Their
30/20/14 um convergence results remain valid historical numerical evidence
for that generator, and the rejected `pMesh` experiment remains useful audit
evidence. They are not the production BO meshes after this switch. In
particular, the former accepted 14 um reference contained 24 tetrahedral cut
transitions, so it does not satisfy the new zero-tetrahedra policy.

Mesh and bad-cell inspection uses foamlib/Python VTK tooling; no ParaView GUI
or headless ParaView rendering is required.
