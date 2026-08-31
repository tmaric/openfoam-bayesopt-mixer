# Slide decks

Two standalone decks for the **Delft Workshop on Bayesian Optimization for CFD**.

| file | what |
|---|---|
| [`bayesian-optimization-cfd-theory.html`](bayesian-optimization-cfd-theory.html) | **Theory.** Gaussian distribution → multivariate Gaussian → Gaussian process → Bayes and GP conditioning → acquisition functions → multi-objective and hypervolume → PINN surrogates and concurrent workflows. |
| [`bayesian-optimization-cfd-tutorial.html`](bayesian-optimization-cfd-tutorial.html) | **Hands-on.** The planar alternating-deflector micromixer end to end: the problem, the method, building it with agents, the results, a six-exercise assignment with answers, and a reproduction appendix. |

## No installation, no plugins, no build step

Both decks are plain HTML plus a local `slides.css` and `slides.js`. There is
**no reveal.js, no plugin, no CDN and no package manager** — the slide engine is
about 90 lines at the top of `slides.js`. Consequences worth knowing:

- **double-click either file** and it opens; `file://` works, no server needed;
- they work **offline**, on any modern browser, on any machine;
- nothing to install for anyone you hand them to.

To serve the folder instead (useful when presenting from another device):

```bash
./serve.sh          # python3 -m http.server 8000
```

### Navigating

| key | |
|---|---|
| `→` `space` `PageDown` | next slide |
| `←` `PageUp` | previous slide |
| `Home` `End` | first / last |
| `N` | toggle speaker notes |

The URL carries the slide number (`#/23`), so any slide can be linked directly.

## Figures are computed, not drawn

`slides.js` renders every diagram at load time from real algebra or real data:

- **GP posterior and UCB** — exact squared-exponential posterior over five
  observations; the mean interpolates the data and the variance collapses there
  because the algebra is real, not because it was drawn that way.
- **The κ sweep** — the Forrester function, one shared posterior, three values
  of κ. κ = 0 picks 0.63 (next to a sample it already has), κ = 3 lands on the
  true optimum at 0.758.
- **GP prior samples** — five functions drawn by Cholesky-factorising the kernel
  matrix and multiplying by standard normals.
- **Unit-cell geometry** — drawn *from* the parameter vector shown above it,
  following the section layout in `../FlowCase/alternating_deflector_cad.py`.
- **Pareto front and hypervolume** — the dominated staircase against the
  declared reference point, and the area one candidate would add.
- **The campaign Pareto chart** — the twelve corrected designs, read from
  `../results/corrected_boundary_v3/all_samples.csv`.

Editing a figure means editing the data or the kernel, not moving pixels.

## The assignment

Part five of the tutorial deck is six exercises alternating **task → answer**,
with the exact commands. Tasks 1–4 run CFD (minutes each); tasks 5–6 run none,
reusing the twelve designs the campaign already paid for. The whole assignment
runs **inside the container** — building the image is the only step that cannot:

```bash
apptainer shell --bind "$PWD" apptainer/padm.sif
```

Reference solution for tasks 5 and 6: [`exercises/surrogate.py`](exercises/surrogate.py).

## Assets

`assets/pareto_animation.mp4` and `.gif` come from a workflow rule, not by hand:

```bash
cd .. && snakemake visualize --workflow-profile profiles/local \
    --config results_dir=results/corrected_boundary_v3
cp results/corrected_boundary_v3/visualizations/pareto_animation.* docs/assets/
```

## Source material

`2026-BayesianOptimization-CFD.pptx` and `.pdf` are the **input** lecture slides
the theory deck was derived from. They are deliberately **not tracked** (see
`.gitignore`) and are not modified by anything here.

## Data provenance

- **corrected** baselines: `../results/corrected_boundary_v3_baselines/`
- **corrected** gated campaign, the only data shown: `../results/corrected_boundary_v3/`
- the archived 28-point campaign and `../results/verified_flux_sequential_v2/`
  predate the boundary repair. Retained to document the defect; **never pooled**
  with the corrected results.

See `RESEARCH_PLAN.md` for the phase gates and `../README.md` for how to run the
study.
