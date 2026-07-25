# Research plan: verified sequential optimization

## Decision already made

The device is studied as a **Planar Alternating-Deflector Micromixer (PADM)**,
not as a split-and-recombine mixer. The present 2-D topology does not split,
permute, and restack fluid layers.

The historical 28-point campaign is an exploratory archive. Its numerical
scheme, objective definition, design variables, and convergence policy differ
from the verified campaign, so old and new targets are not interchangeable.

## Phase 0 — workflow pilot

Status: complete for two distinct Sobol geometries.

- Use second-order bounded `limitedLinear 1` scalar convection.
- Require explicit SIMPLE flow convergence, with up to 2000 iterations.
- Run at least 200 and at most 600 scalar iterations, then require the final 50
  samples to have spans no greater than `1e-4` in flux-weighted mean and
  segregation intensity.
- Use PBiCGStab/DILU for the high-Peclet scalar matrix and equation relaxation
  0.7.
- Minimize kinematic pressure drop and flux-weighted intensity of segregation.
- Exclude failed or unconverged evaluations from the GP, retaining their
  geometry and failure status for later feasibility modelling.

Pilot observations are recorded in `../QUICKSTART.md`. The pilot establishes
workflow operation only; it does not estimate an optimum.

## Phase 1 — space-filling initialization

Run 32 scrambled-Sobol designs in the six-dimensional feasible box. Execute
strictly one design at a time (`q = 1`), with two OpenFOAM MPI ranks and one
Torch thread.

Advance in bounded invocations so machine use remains controllable:

```bash
python bayes_optimize_sequential.py --max-new-evaluations 1
```

After 8, 16, and 32 successful designs, review:

- mesh-generation and CFD failure rate;
- target distributions and duplicates;
- mass balance and scalar boundedness;
- whether the fixed reference point `(0.02 m²/s², 1.0)` is dominated by all
  valid observations;
- whether any parameter bound is repeatedly selected by the non-dominated set.

If failures cluster geometrically, fit a feasibility classifier or contract
the design bounds before BO. Do not replace failures with arbitrary objective
penalties.

## Phase 2 — sequential multi-objective BO

Fit independent-output exact GPs with normalized six-dimensional inputs,
standardized outputs, and an ARD Matérn-5/2 covariance. Select one candidate
per iteration using qLogNEHVI (`q = 1`), 32 multistarts, 1024 raw samples, and
the fixed physical reference point.

The configured budget is 80 BO observations after the 32-point initialization.
At every 10 successful BO observations, save and review:

- Pareto front and dominated hypervolume;
- GP length scales and standardized residuals;
- candidate distances from earlier designs;
- objective repeatability for at least one re-evaluated design;
- concentration of Pareto points at design-space bounds.

Stop early if the hypervolume gain is practically negligible for 15
consecutive successful evaluations and posterior uncertainty is already small
along the estimated front. Any stopping rule used for a paper must be fixed
before interpreting the final candidates.

## Phase 3 — numerical verification

Do not identify a publishable optimum directly from the BO mesh. Select at
least five Pareto designs spanning low pressure loss to high mixing, plus a
straight channel and a symmetric-deflector baseline.

For each selected design:

1. repeat on at least three systematically refined meshes;
2. compare `limitedLinear 1` with another bounded high-resolution scheme;
3. demonstrate mass conservation and outlet scalar boundedness;
4. repeat enough evaluations to estimate numerical/run-to-run uncertainty;
5. report flux-weighted mixing index, pressure loss, Reynolds number,
   Schmidt/Péclet number, channel length, and computational cost.

Use a refinement-based uncertainty estimate (for example GCI when the
asymptotic regime is demonstrated). Re-evaluate Pareto dominance using the
verified objectives and their uncertainty intervals.

## Phase 4 — scientific comparison

Compare all designs at matched hydraulic diameter, channel length or residence
time, Reynolds number, diffusivity, inlet condition, and mixing-index
definition. At minimum include:

- an unobstructed straight channel;
- a symmetric version of this deflector topology;
- one established planar passive-mixer benchmark reproducible with the same
  solver;
- literature data only where metric and operating-condition conversions are
  defensible.

A fair performance plot should show mixing index against pressure drop or
pumping power, not mixing alone. Report dimensional pressure drop in pascals
alongside the native kinematic value.

## Phase 5 — publication gate

The current legacy optimum is **not publication-ready as a global optimum**.
A Chemical Engineering Science-level claim would require, at minimum:

- completion or defensible early stopping of the verified sequential campaign;
- mesh/scheme independence and quantified uncertainty;
- meaningful passive-mixer baselines at matched conditions;
- a clear physical mechanism supported by velocity/scalar fields;
- reproducible code, parameter bounds, random seed, objective definitions, and
  complete successful/failed evaluation history;
- preferably 3-D confirmation or experimental validation, especially if the
  intended claim concerns SAR-like lamination.

Without those items, the defensible contribution is a reproducible 2-D
multi-objective numerical design study and a corrected methodology—not a
validated universal passive-mixer optimum.
