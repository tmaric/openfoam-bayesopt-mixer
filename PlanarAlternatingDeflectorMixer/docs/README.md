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

- 28 stored samples in `../results/all_samples.csv`;
- scalar-field renders from samples `00022` and `00027`;
- current geometry, mesh, solver, and BO source under the parent study folder;
- historical SAR slide material under `../../docs/obsidian/02-Technical-Notes/`.
