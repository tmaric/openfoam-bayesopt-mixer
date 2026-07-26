# Research plan: corrected, gated sequential optimization

## Scientific status

This is a Planar Alternating-Deflector Micromixer (PADM), not a true
split-and-recombine mixer. The legacy and `verified_flux_sequential_v2`
campaigns are invalid: internal x-normal obstacle faces were assigned to the
inlet and outlet patches. Their pressure and scalar objectives must never be
combined with corrected data.

## Phase 0 — geometry and physics repair

Status: complete.

- Classify inlet/outlet by complete face location at `x=0` and `x=L`.
- Emit a CAD geometry manifest with boundary bounds and areas.
- Reconstruct mesh patch polygons independently and reject incorrect location,
  area, type, short edges, concave cells, or other failed `checkMesh` tests.
- Require `Q = U_mean A` within 3% and relative mass imbalance below `1e-5`.
- Use mesh-aware deflector endpoint floors, curve tessellation, and centre
  splitter transitions.
- Accept obstructed flow only after SIMPLE residual convergence. The straight
  baseline additionally permits a pressure/Ux convergence test because its
  physically zero transverse component has an ill-conditioned relative
  residual.
- Advance scalar transport in resumable chunks and require objective stability.

## Phase 1 — matched baselines

Status: complete at `Re=10`, `Sc=1000`, `H=1 mm`, and `L=24 mm`.

| Baseline | Pressure (Pa) | Pressure ratio | Mixing index |
|---|---:|---:|---:|
| Straight | 2.8737 | 1.000 | 0.1003 |
| Symmetric deflectors | 10.8116 | 3.762 | 0.0901 |
| Strong alternating | 30.5091 | 10.617 | 0.1450 |

The analytical straight-channel pressure drop is 2.8748 Pa; the CFD error is
0.039%. The strong alternating reference exceeds the predeclared 20 Pa budget
without approaching competitive mixing.

## Phase 2 — corrected topology feasibility screen

Status: complete — NO-GO. All twelve scrambled-Sobol designs passed every
numerical validation, with zero failures. The best mixing index was 0.1698 at
34.0596 Pa (pressure ratio 11.852). The best design below the 20 Pa budget was
`00003`, with mixing index 0.1126 at 16.4689 Pa. The gate required 0.60.

```bash
python research_sequence.py next --max-new-evaluations 1
```

The predeclared go/no-go gate requires:

- twelve successful corrected designs;
- best flux-weighted mixing index at least 0.60;
- failed-evaluation fraction no greater than 0.25.

The driver wrote `results/corrected_boundary_v3/screening_gate.json`; the
tracked result summary is `research/corrected_screening_summary.yaml`. The
no-go means the planar topology must be adapted before spending the full BO
budget.

## Phase 3 — topology adaptation after a no-go

Status: next action. The corrected screen has triggered this phase.

Do not merely enlarge the current amplitude bounds. Preserve the matched
operating point and first test one of these mechanisms:

1. an inclined/parallelogram-barrier or modified-Tesla planar unit with angle,
   barrier length, throat ratio, lateral offset, pitch, and unit count;
2. preferably, a genuine three-dimensional SAR/crossing-channel or staggered
   groove topology with groove/channel angle, depth, width, overlap,
   asymmetry, pitch, and unit count.

Run the same three baselines and twelve-point feasibility screen for each new
topology version. Do not use the current 2-D topology as a numerical
low-fidelity model for a 3-D mechanism unless cross-fidelity correlation is
demonstrated.

## Phase 4 — full sequential BO after a pass

The full stage must be explicitly requested:

```bash
python research_sequence.py optimization --max-new-evaluations 1
```

It expands the initialization to 32 successful Sobol designs, then performs
80 `q=1` qLogNEHVI evaluations using independent-output exact GPs, normalized
inputs, standardized outputs, and ARD Matérn-5/2 kernels. The pressure
objective is normalized by the validated straight channel. The acquisition
reference corresponds to 20 Pa and segregation intensity 1.0.

Failed CFD cases retain geometry and failure status but receive no artificial
objective penalty. If failures cluster after enough evidence exists, add an
explicit feasibility classifier before resuming acquisition.

## Phase 5 — numerical and robust verification

Select at least five Pareto designs spanning the trade-off, plus all baselines.
For each:

1. run at least three systematically refined meshes;
2. compare the checked-in bounded high-resolution scalar scheme with another
   bounded high-resolution scheme;
3. report mass conservation and scalar bounds;
4. quantify discretization uncertainty and re-evaluate Pareto dominance with
   uncertainty intervals;
5. perturb manufacturable dimensions and report expected and worst-case
   performance under fabrication tolerances;
6. if claiming operating robustness, repeat at predeclared Reynolds numbers
   such as 1, 5, 10, and 20.

## Publication gate

A corrected 2-D campaign is not automatically a Chemical Engineering Science
paper. A strong submission requires a physical mixing mechanism that exceeds
matched baselines, a verified mixing-pressure Pareto front, numerical
uncertainty, mechanistic flow/Lagrangian evidence, complete reproducibility,
and preferably 3-D or experimental validation. BO is the search method, not by
itself the physical novelty.
