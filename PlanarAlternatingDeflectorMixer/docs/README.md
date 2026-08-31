# PADM slide deck

Reveal.js presentation for the Planar Alternating-Deflector Micromixer study.

```bash
./serve.sh          # python3 -m http.server 8000
```

Open <http://localhost:8000/>. Arrow keys navigate, `Esc` is the overview, `S`
opens speaker notes, `F` is full screen.

## Structure

| Part | Slides | Content |
|---|---|---|
| **1 · The problem** | the device, why Re = 10 / Sc = 1000 makes mixing hard, the six design parameters, both objectives, the constraints |
| **2 · Bayesian optimization** | Bayes' theorem and the probability square, GPs as priors over functions, the closed-form mean and variance, UCB, optimising the acquisition, data efficiency |
| **3 · Building it** | development with agents, test-driven from the top down, Snakemake, local/SLURM profiles, BoTorch vs Ax vs scikit-learn, the cost asymmetry, batch BO, semi-automatism, multi-objective fronts |
| **4 · Results** | the campaign animation, the corrected Pareto front, the predeclared NO-GO |

## Figures are computed, not drawn

`slides.js` renders two SVGs at load time:

- **the GP panel** — a squared-exponential kernel and an *exact* posterior over
  five synthetic observations, plus its UCB acquisition. The mean interpolates
  the data and the variance collapses there because the algebra is real; the
  argmax it marks (x ≈ 0.771) is genuinely where UCB points.
- **the Pareto chart** — the twelve corrected screening designs, straight from
  `../results/corrected_boundary_v3/all_samples.csv`, with the three matched
  baselines and the predeclared 0.60 gate.

Editing either figure means editing the data or the kernel, not moving pixels.

## Assets

`assets/pareto_animation.mp4` and `.gif` are produced by a workflow rule, not by
hand:

```bash
cd .. && snakemake visualize --workflow-profile profiles/local \
    --config results_dir=results/corrected_boundary_v3
cp results/corrected_boundary_v3/visualizations/pareto_animation.* docs/assets/
```

MP4 export needs ffmpeg (present in the Apptainer image); the GIF is Pillow and
always works.

## Data provenance

- **corrected** baselines: `../results/corrected_boundary_v3_baselines/`;
- **corrected** gated campaign (the only data shown): `../results/corrected_boundary_v3/`;
- the archived 28-point campaign and `../results/verified_flux_sequential_v2/`
  predate the boundary repair. They are retained to document the defect and are
  **never pooled** with the corrected results.

The deck is self-contained apart from Reveal.js and KaTeX, pinned from jsDelivr.
All paths are relative, so the folder moves with the repository.

See `RESEARCH_PLAN.md` for the topology, verification and publication gates, and
`../README.md` for how to run the study.
