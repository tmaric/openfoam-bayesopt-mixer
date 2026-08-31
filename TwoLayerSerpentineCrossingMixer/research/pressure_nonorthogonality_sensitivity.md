# Pressure non-orthogonality sensitivity

Status on 2026-07-28: **closed**. The production setting of one
non-orthogonal pressure correction is sufficient for the fine Re=1 review
case. Increasing the setting to two changes neither the reported pressure drop
nor the converged velocity and flux fields materially. Scalar transport was
therefore not rerun.

## Controlled comparison

The comparison reuses the exact validated 14 um, 330,440-cell review mesh and
changes only

```text
nNonOrthogonalCorrectors 1;
```

to

```text
nNonOrthogonalCorrectors 2;
```

Both flow solutions use `simpleFoam`, four MPI ranks, `consistent yes`,
`Gauss linear corrected` Laplacians, and corrected surface-normal gradients.
The mesh has maximum non-orthogonality `37.713794` degrees, average
non-orthogonality `4.0228356` degrees, and maximum skewness `1.4186347`.

| Quantity | One correction | Two corrections | Difference |
|---|---:|---:|---:|
| SIMPLE iterations | 112 | 111 | -1 |
| Pressure drop (m2/s2) | 0.015303792 | 0.015303792 | 0 |
| Pressure drop (Pa) | 15.257880624 | 15.257880624 | 0 |
| Patch mass-balance relative error | 5.03787e-9 | 1.95120e-9 | immaterial |
| Final global continuity error | 2.26503e-9 | 8.38292e-10 | immaterial |
| Internal `U` relative L2 difference | - | - | 8.72444e-7 |
| Internal `phi` relative L2 difference | - | - | 8.66416e-7 |
| Internal `U` relative maximum difference | - | - | 1.53907e-6 |
| Internal `phi` relative maximum difference | - | - | 1.42504e-6 |

The patch-integrated balance is the physical conservation check. The
iteration-history cumulative continuity value is not used to compare the two
runs because it accumulates signed solver corrections along different
convergence paths.

The two-correction pressure change is below the predeclared `0.5%`
materiality threshold; it is zero at the available log precision. The flow
field differences are approximately one part per million. Propagating this
change through a 2,400-iteration scalar solve would not be a responsible use
of the validation budget.

## Decision

Keep

```text
nNonOrthogonalCorrectors 1;
```

for reproduction and multifidelity BO evaluations. More correctors add work
to each SIMPLE iteration without addressing the remaining Re=1 mixing-index
discrepancy. The unresolved issue remains the reconstructed inter-layer
aperture/topology or the scientific framing of the geometry, not pressure
coupling convergence.

Runtime products are intentionally ignored under
`results/numerical_sensitivity/nonorthogonal_correctors/`. The comparison
summary is stored there as `Re1/fine/n2/flow_sensitivity.json`; the accepted
one-correction baseline remains under
`results/reproduction/review_second_order/Re1/fine/`.
