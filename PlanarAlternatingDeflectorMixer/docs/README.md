# PADM Bayesian-optimization slides

This folder contains the Reveal.js presentation for the
Planar Alternating-Deflector Micromixer study.

Serve it from this directory so browsers load all relative assets correctly:

```bash
./serve.sh
```

Open <http://localhost:8000/>. Use the arrow keys to navigate, `Esc` for the
overview, `S` for speaker view, and `F` for full screen.

The deck is self-contained except for the pinned Reveal.js and KaTeX files,
which are loaded from jsDelivr. All study paths and image references are
relative; the folder can be moved with the repository.

Data provenance:

- corrected matched baselines under the ignored
  `../results/corrected_boundary_v3_baselines/` directory when generated;
- the corrected gated campaign under `../results/corrected_boundary_v3/`;
- the archived 28-point campaign in `../results/all_samples.csv`;
- invalid pilot observations under `../results/verified_flux_sequential_v2/`,
  retained only to document the former boundary-classification defect;
- scalar-field renders from samples `00022` and `00027`;
- current geometry, mesh, solver, and BO source under the parent study folder;
- historical SAR slide material under `../../docs/obsidian/02-Technical-Notes/`.

The checked-in campaign is strictly sequential (`q = 1`) and cannot enter full
BO until its corrected twelve-point feasibility gate passes. See
`RESEARCH_PLAN.md` for the topology, numerical-verification, and publication
gates.
