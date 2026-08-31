# Second-order passive-scalar validation

The M10 passive scalar is the dimensionless field `T`. Its production
discretization is now:

```text
gradT           cellLimited pointCellsLeastSquares 1;
div(phi,T)      Gauss linearUpwind gradT;
```

The earlier `bounded Gauss limitedLinear 1` result set is retained only as a
scheme-sensitivity reference. It must not be used as the production M10 scalar
result.

All cases below are the completed pre-switch validation dataset. They use the
original Hossain nine-unit protocol, diffusivity `1e-11 m2/s`, the accepted
14 um cfMesh Cartesian mesh, 2,400 scalar pseudo-iterations, and a final
50-iteration functional-stability gate. Re=1, 20, and 40 were completed
sequentially. The final matrix uses four MPI ranks; Re=1 was also run at two
ranks and its mixing index agreed within `1e-7`. These values preserve the
scheme study but are not snappyHexMesh convergence evidence; the production
mesh is now the qualified 13 um OpenFOAM-v2606 snappy mesh.

## Results

| Re | Pressure (Pa) | Area MI | Flux MI | Flux segregation | min(T) | max(T) | Paper criterion | Result |
|---:|---:|---:|---:|---:|---:|---:|---|---|
| 1 | 22.542 | 0.990087 | 0.990275 | 9.45710e-5 | -1.926e-5 | 1.000021 | approximately 0.99 for Re <= 10 | pass |
| 20 | 518.921 | 0.962381 | 0.968752 | 9.76435e-4 | -4.860e-6 | 1.000049 | at least 0.96 over Re=0.2--120 | pass |
| 40 | 1269.350 | 0.980422 | 0.979531 | 4.18964e-4 | -3.706e-4 | 1.000372 | at least 0.96 over Re=0.2--120 | pass |

All final-window intensity and mean-concentration spans pass the tightened
gate; the recorded spans are zero at Re=1 and 20 and `1e-9` at Re=40.

`linearUpwind` is not bounded. The largest source-benchmark excursion is at
Re=40: `T=-0.000371` to `1.000372`, or 0.0372% outside the physical interval.
The first 13 um design-space pilot anchor later reported one wall-adjacent
minimum of `-0.00109384`, maximum `1.0000149`, outlet mean `0.500000`, and a
stationary final-50 outlet objective. Before evaluating the remaining pilot
designs, the local excursion cap was declared as `0.002` (0.2%). Values beyond
that limit fail; the field is neither clipped nor replaced. Bounds continue to
be checked with `postProcess fieldMinMax(T)` after every run.

## Scheme effect

At Re=1, replacing the bounded scheme changes the fine area mixing index from
0.992734 to 0.990087. At Re=40 it changes from 0.987007 to 0.980422. The larger
high-Re reduction is consistent with less numerical mixing from the requested
second-order transport discretization.

The paper comparison above validates the reported mixing behavior. Exact
pressure-drop validation still requires digitizing the original paper's
pressure curve or obtaining its underlying data.

## Separate six-unit review protocol

The Raza-review comparison is a separate six-unit, 5.05 mm protocol with
diffusivity `1e-10 m2/s`; it must not be pooled with the nine-unit results
above. It was rerun independently with the same production scalar scheme, the
legacy accepted 14 um cfMesh, 2,400 pseudo-iterations, and four MPI ranks. The
table must be regenerated on the qualified snappy meshes before making a new
cross-generator comparison.

| Re | Pressure (Pa) | Area MI | Flux MI | min(T) | max(T) | Review target | Result |
|---:|---:|---:|---:|---:|---:|---|---|
| 1 | 15.258 | 0.956333 | 0.957138 | -5.552e-5 | 1.000062 | MI 0.915; 16.3 Pa | pressure pass; MI fail |
| 10 | 160.829 | 0.917259 | 0.943334 | -1.654e-4 | 1.000130 | no tabulated target | midpoint only |
| 20 | 354.132 | 0.891928 | 0.907025 | -4.336e-6 | 1.000005 | MI 0.901; 390 Pa | pass |

Pressure errors are -6.39% at Re=1 and -9.20% at Re=20, within the declared
10% tolerance. Mixing-index errors are +0.04133 and -0.00907 respectively;
only Re=20 is within the absolute 0.03 tolerance. All numerical acceptance
checks pass, so the remaining Re=1 issue is a geometry, outlet-station, or
metric-definition audit rather than an unconverged-run indication. See
`reproduction_status.md` for the resulting BO gate.

The subsequent audit closed the outlet-station and metric questions but found
a monotone resolution trend at Re=1: area MI is `0.913232`, `0.942197`, and
`0.956333` at 30, 20, and 14 um respectively. Pressure simultaneously moves
toward the review target. The target-fidelity mixing discrepancy must therefore
be treated as an inter-layer topology/mesh-representation issue. Full evidence
and the straight-channel normalization are in `review_protocol_audit.md`.
