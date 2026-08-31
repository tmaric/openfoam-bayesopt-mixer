# M10-inspired two-layer serpentine-crossing micromixer

This study implements and optimizes a genuinely three-dimensional mixer
inspired by the M10 micromixer reported by Hossain et al. in *Chemical
Engineering Journal* 327 (2017), 268-277, DOI
`10.1016/j.cej.2017.06.106`. It does not claim exact geometric reproduction:
the published sources identify nominal dimensions and connection locations but
do not provide CAD/mask data, independent aperture dimensions, unique lead
transitions, or corner treatment.

The committed PADM study remains a negative planar-topology result. This new
folder is independent: no PADM result is used to train an M10 surrogate.

## Current milestone

The repository now contains:

- a parameterized, portable CadQuery M10-inspired fluid domain;
- separate original-paper and standardized-review validation protocols;
- a six-dimensional extension of the published three-ratio parameterization;
- an explicit coarse/fine multifidelity contract;
- a four-rank OpenFOAM-v2606 `snappyHexMesh` pipeline with zero tetrahedra;
- a passing 12-mesh coarse/fine preflight over the six paired BO anchors;
- a passing end-to-end coarse Re=20 flow/scalar smoke test under v2606;
- static tests that prevent the fidelity level from being treated as an
  after-the-fact mesh label.

The source-matched coarse curve and fine Re=1, 10, and 40 benchmarks are
complete on the superseded cfMesh pipeline. The separate six-unit review was
also rerun at fine resolution for Re=1, 10, and 20 with the production
second-order scalar scheme. Re=20 passes both literature gates; Re=1 pressure
passes but its area mixing-index error remains `+0.04133` against the allowed
`0.03`. Those values remain historical validation evidence. Because the study
is explicitly an M10-inspired reconstruction, literature values are external
benchmarks rather than an exact-geometry acceptance target. Production
snappyHexMesh transport and the paired-fidelity rank-correlation gate remain
to be completed. See `research/reproduction_status.md`.

The production passive-scalar discretization is second order:
`Gauss linearUpwind gradT`, with
`gradT cellLimited pointCellsLeastSquares 1`. Fine Re=1, 20, and 40 validation
is documented in `research/second_order_scalar_validation.md`.

The short inlet and outlet lead construction is an explicitly recorded
reconstruction assumption. The generated CAD has been audited for the intended
two-layer crossing topology; the exact low-Re outlet station and lead
construction still require source comparison before BO.

A controlled fine Re=1 sensitivity increased the pressure non-orthogonal
correctors from one to two on the identical mesh. Pressure drop was unchanged
and the relative L2 changes in `U` and `phi` were below `9e-7`, so one
corrector remains the production setting. See
`research/pressure_nonorthogonality_sensitivity.md`.

## Generate and inspect the reference CAD

Run from Ubuntu WSL with the repository's Python environment:

```bash
cd TwoLayerSerpentineCrossingMixer
python FlowCase/two_layer_serpentine_crossing_cad.py
python render_geometry.py
```

Outputs are written below `generated/reference/` using paths relative to this
study. The STEP file is intended for visual/topological inspection; the
multi-solid ASCII STL carries `inlet1`, `inlet2`, `outlet`, and `walls` patch
names into `snappyHexMesh`. The actual 3-D workflow uses a uniform `blockMesh`
background, explicit feature extraction, four-rank snapping, patch merging,
and transport-oriented `checkMesh` criteria. Both fidelities are
hexahedral-dominant and contain zero tetrahedra; body-fitted prism/polyhedral
cut cells are measured and capped. Details are in
`research/mesh_qualification.md`.

Validate the optimization and fidelity contract:

```bash
python multifidelity_design.py validate-config
python -m unittest discover -s tests -v
python multifidelity_design.py reference --fidelity coarse
python multifidelity_design.py reference --fidelity fine
```

Prepare a portable, fully materialized reproduction case without running CFD:

```bash
python run_case.py \
  --protocol original --fidelity coarse --reynolds 10 \
  --results-dir results/reproduction/original/Re10/coarse \
  --prepare-only
```

Remove `--prepare-only` only after sourcing OpenFOAM-v2606 and building the
repository-local function objects:

```bash
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
../Allwmake
```

Use `--mesh-only` for the first bounded mesh audit. A result directory can be
replaced only with the explicit `--force` flag, and the runner refuses paths
outside this study's `results/` directory. No repository path is embedded in
the generated case, so the checkout can be moved or cloned elsewhere.

## Fidelity is part of BO

The augmented surrogate input is

```text
[H/P, w/P, D/P, b/P, end_inset/w, crossing_phase/b, s]
```

with `s=0` for the 24 micrometre coarse discretization and `s=1` for the 13
micrometre target discretization. Each objective is modeled by a
`SingleTaskMultiFidelityGP`. Sequential, cost-aware
`qMultiFidelityKnowledgeGradient` selects a design and fidelity together after
randomized augmented-Chebyshev reduction of the two objectives. This is a
multifidelity BO model, not coarse screening followed by an unrelated fine
rerun.

The BO operating point is fixed at Re=10. Pressure ratio and residual
segregation are modeled on logarithmic scales; Re=1 and Re=40 are fine-only
robustness checks.

The pressure ratio now uses a per-design fully developed rectangular-duct
reference with the same `w x D` outlet cross-section, total axial length, mean
velocity, and viscosity. The exact series gives `43.5885 Pa` for the four-unit
reference design at Re=10 and is evaluated analytically, so it does not double
the CFD cost of every BO observation.

The initialization contains 24 coarse Sobol designs and six paired fine
anchors. The paired observations identify the coarse-to-fine discrepancy. The
campaign is strictly sequential (`q=1`), uses at most four OpenFOAM ranks and
one Torch thread, and is capped by a fine-equivalent cost budget. Eight final
Pareto candidates are reevaluated at fine resolution; only those fine values
may support scientific claims.

The mesh-only paired-anchor gate is complete and passes all 12 snappy meshes.
Run `python multifidelity_pilot.py preflight` to reapply the current fixed
policy to its retained runtime evidence before starting transport.

## Mandatory gates before BO

1. Confirm the reconstructed topology and lead transitions visually.
2. Preserve the original and six-unit literature comparisons as external
   benchmarks, without claiming exact geometric reproduction.
3. Pass the fixed snappyHexMesh reference and six-design coarse/fine mesh
   qualification (complete).
4. Demonstrate adequate coarse/fine rank correlation and quantify bias for
   both objectives.
5. Complete fine mesh convergence for the reference and selected designs.

The OpenFOAM reproduction workflow and straight-pressure normalization are
operational, and the work is explicitly framed as an M10-inspired
reconstruction. BO execution remains gated by production-snappy transport,
the paired coarse/fine objective-correlation pilot, and selected-design mesh
convergence. This avoids spending the campaign budget before the new mesh
generator is numerically qualified.
