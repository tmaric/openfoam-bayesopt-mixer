# Project Charter

## Project
`openfoam-bayesopt-mixer`

## Objective
Build a reproducible, working multi-objective Bayesian optimization workflow for a 2D laminar two-stream mixer with parameterized internal baffles.

## Final Milestone
A working end-to-end BO workflow that:
1. proposes design parameters,
2. generates geometry and mesh,
3. runs CFD + scalar transport,
4. computes both objectives,
5. updates a multi-objective BO model,
6. produces a Pareto set/front with reproducible outputs.

## Scope
- CFD benchmark in OpenFOAM.
- CAD/mesh generation via Gmsh Python API.
- Multi-objective optimization over baffle geometry parameters.
- Weekly execution tracking and milestone updates.

## Out of Scope (for now)
- Turbulent or 3D variants.
- Manufacturing optimization beyond minimum geometric constraints.
- Production-grade UI.

## Documentation Strategy
- Keep high-level intent in `01 Project`.
- Keep weekly activity in `02 Planning/Weekly`.
- Keep detailed technical implementation notes in `04 Technical Notes` and later consolidate into repository `README.md`.
