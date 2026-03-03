# Milestones

Final objective: working, reproducible multi-objective Bayesian optimization workflow by **2026-04-13**.

## Deadline Plan
1. M0 - Planning baseline (done)
- Status: Completed
- Date: 2026-03-02
- Outcome: benchmark spec + milestones + weekly tracking in place.

2. M1 - Geometry generator + constraints
- Status: Planned
- Target week: 2026-03-02 to 2026-03-06
- Deliverable: generate feasible geometry from `theta` in Gmsh Python API.

3. M2 - Mesh + OpenFOAM baseline case
- Status: Planned
- Target week: 2026-03-09 to 2026-03-13
- Deliverable: stable flow + scalar run for one design.

4. M3 - Objective extraction
- Status: Planned
- Target week: 2026-03-16 to 2026-03-20
- Deliverable: automated `J_mix` and `J_delta_p` computation from solver output.

5. M4 - Single-point automation
- Status: Planned
- Target week: 2026-03-23 to 2026-03-27
- Deliverable: one command evaluates a single `theta` end-to-end.

6. M5 - Multi-objective BO loop
- Status: Planned
- Target week: 2026-03-30 to 2026-04-03
- Deliverable: BO proposes points and updates model with objective values.

7. M6 - Stabilization + reproducibility + handoff
- Status: Planned
- Target window: 2026-04-06 to 2026-04-13
- Deliverable: stable end-to-end workflow, reproducibility notes, README consolidation.

## Weekly Alignment Rule
Every weekly note must include:
- `Primary milestone:`
- `Secondary milestone (optional):`
- `Milestone status change:`
