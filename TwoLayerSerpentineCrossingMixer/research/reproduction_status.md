# Reproduction status and BO gate

The numerical reproduction results below were generated with the superseded
cfMesh pipeline. They remain an audit trail but must not be presented as
OpenFOAM-v2606 snappyHexMesh convergence results. The replacement mesh-only
reference and six-design qualification are complete; see
`mesh_qualification.md`.

Status on 2026-07-28: **BO remains gated**. The original-paper trend is
reproduced, and the separate six-unit review protocol has now been rerun with
the production second-order scalar scheme. Re=20 passes both review gates, but
the Re=1 area-weighted mixing-index error is `+0.04133`, outside the declared
absolute tolerance of `0.03`. Starting the 30-case multifidelity
initialization before resolving that discrepancy would spend the budget on an
incompletely validated CAD or measurement protocol.

The follow-up audit in `review_protocol_audit.md` confirms that the source's
mixing-index definition is equivalent to the implemented area-weighted metric
and that plausible 80 um shifts of the inlet/outlet station change coarse MI by
at most `0.00333`. Neither explains the discrepancy. A three-level Re=1 study
instead shows area MI increasing monotonically from `0.913232` at 30 um to
`0.942197` at 20 um and `0.956333` at 14 um. The coarse agreement with the
review target is accidental.

## Original Hossain protocol

The original protocol uses nine units and scalar diffusivity `1e-11 m2/s`.
The area and flux columns are kept separate; the BO objective is flux-weighted
segregation, while literature mixing-index comparisons use area weighting.

| Re | Fidelity | Pressure (Pa) | Area MI | Flux MI | Flux segregation |
|---:|---|---:|---:|---:|---:|
| 0.2 | coarse | 4.184 | 0.981271 | 0.986454 | 1.83483e-4 |
| 1 | coarse | 20.954 | 0.981000 | 0.986210 | 1.90175e-4 |
| 1 | fine | 22.542 | 0.992734 | 0.993341 | 4.43463e-5 |
| 10 | coarse | 220.818 | 0.988807 | 0.987936 | 1.45544e-4 |
| 10 | fine | 236.775 | 0.991881 | 0.995116 | 2.38506e-5 |
| 15 | coarse | 345.884 | 0.990231 | 0.992314 | 5.90705e-5 |
| 40 | coarse | 1171.625 | 0.998637 | 0.998732 | 1.60660e-6 |
| 40 | fine | 1269.350 | 0.987007 | 0.987250 | 1.62556e-4 |

The reported approximately 0.99 mixing for Re up to 10 is reproduced within
the declared absolute tolerance. Re=40 exposes severe coarse-grid numerical
mixing: pressure remains only 8.3% low, but coarse residual segregation is
about 101 times smaller than the fine value.

This evidence fixes the BO operating point at Re=10. Re=1 and Re=40 are
fine-only robustness validations, not unmodelled variations inside one GP.
Positive objectives are log transformed so near-unity mixing indices do not
hide orders-of-magnitude segregation errors.

## Separate Raza review protocol

This protocol uses six units, 5.05 mm axial length, and diffusivity
`1e-10 m2/s`. Errors below use the fine, area-weighted mixing index.

| Re | Target MI | Fine MI | MI error | Target pressure (Pa) | Fine pressure (Pa) | Pressure error |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.915 | 0.956333 | +0.041333 | 16.3 | 15.258 | -6.39% |
| 20 | 0.901 | 0.891928 | -0.009072 | 390.0 | 354.132 | -9.20% |

Pressure passes the 10% gate at both fine points. Re=20 mixing passes the
absolute 0.03 gate; Re=1 mixing does not. The Re=10 midpoint, for which the
review gives no tabulated target, produces area MI `0.917259`, flux MI
`0.943334`, and pressure drop `160.829 Pa`.

All three fine cases used the accepted 14 um hexahedral-dominant mesh, four MPI
ranks, 2,400 scalar pseudo-iterations, and the production `linearUpwind gradT`
scheme. Mass balance, final-window stability, and automated scalar-bound gates
pass. The largest scalar excursion is 0.0165% at Re=10, below the configured
0.1% cap. Compared with the superseded bounded-scheme results, the second-order
scheme moves both review mixing indices downward: Re=1 from `0.964434` to
`0.956333`, and Re=20 from `0.919337` to `0.891928`.

Correcting the function object's area weighting did not remove the low-Re
discrepancy because the outlet face areas are nearly uniform. The six-unit
protocol remains separate from the original nine-unit study and its results
must not be pooled with that data set.

## Pressure-correction sensitivity

The fine Re=1 review mesh has maximum non-orthogonality `37.713794` degrees
and average non-orthogonality `4.0228356` degrees. A controlled four-rank flow
rerun increased `nNonOrthogonalCorrectors` from one to two on the identical
mesh. Pressure drop remained `15.257880624 Pa`; the relative L2 changes in
internal `U` and `phi` were `8.72e-7` and `8.66e-7`, respectively. The
two-correction result passed patch mass balance at `1.95e-9`.

This closes pressure non-orthogonality correction as a cause of the low-Re
mixing discrepancy. One corrector remains the production setting, and the
scalar solve was not repeated because the pressure change was zero, well below
the predeclared `0.5%` materiality threshold. Full details are recorded in
`pressure_nonorthogonality_sensitivity.md`.

## Required next work

1. The reconstructed openings are now quantified as `1.12350 mm2` over seven
   vertical segments plus `0.557822 mm2` over six X crossings. Verify from the
   original model, fabrication mask, or author data whether the entire
   projected overlaps were open in the CEJ simulations.
2. Decide whether the study will claim an exact M10 reproduction or an
   M10-inspired reconstruction; the latter requires an explicit gate and claim
   revision rather than tuning numerics to the review value.
3. Run a small paired coarse/fine rank-correlation pilot only after that
   geometry decision. The current two fidelities differ materially at Re=1.
4. The design-matched analytical straight-channel pressure reference is now
   implemented and recorded in new objectives; validate its definition during
   the pilot rather than adding one CFD solve per BO sample.
5. Once all gates pass, generate the 24 coarse Sobol designs and six paired
   fine anchors at fixed Re=10, evaluate them strictly sequentially, and fit
   the explicit-fidelity GP. Final Pareto claims remain fine-only.

Run `python summarize_reproduction.py` to regenerate the complete result table
from ignored runtime products.
