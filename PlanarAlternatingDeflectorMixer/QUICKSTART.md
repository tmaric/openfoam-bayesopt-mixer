# Quickstart — Planar Alternating-Deflector Micromixer

## What is being optimized

This is a two-dimensional passive micromixer with a centre baffle and five
repeated pairs of cosine wall deflectors. The stronger deflector alternates
between the top and bottom wall. The geometry bends and stretches the single
inlet scalar interface; it does not implement a three-dimensional
split-and-recombine network.

All CAD dimensions are normalized by `H = 1` and converted to SI units by
`scale = 1e-3`.

### Latent BO coordinates

| Coordinate | Bounds | CAD mapping |
|---|---:|---|
| `w_s` | 0.25–0.45 | direct half-gap parameter |
| `t_s` | 0.04–0.15 | direct centre-baffle thickness |
| `t_m_ratio` | 0–1 | maps to a mesh-safe `t_m` interval |
| `L_c` | 0.40–2.40 | direct interaction length |
| `L_s_ratio` | 0–1 | partitions the remaining length into `L_s`, `L_m` |
| `delta_ratio` | 0–1 | maps to an admissible base wall-bias interval |
| `k_ratio` | 0–1 | maps to the feasible linear bias-slope interval |

The physical geometry stored in each sample is

```text
(w_s, t_s, t_m, L_s, L_m, delta, k)
```

with `L_s + L_c + L_m = L_cell = 4`.

### Mesh-aware transforms

The 0–1 coordinates are transformed before CAD generation so proposals satisfy
the current cfMesh resolution assumptions:

- `t_m >= 0.025` and `t_s - t_m >= 0.0125`;
- `0.8 <= L_s, L_m <= 1.8`;
- realized `delta_i = delta + k*xhat_i` stays above `0.04` and below the
  available local gap;
- `max(delta_i) - min(delta_i) <= 0.06`;
- downstream ramp-up is at most `0.05`.

## Objectives

Both objectives are minimized.

### Kinematic pressure drop

```text
J_dp = <p>_inlet - <p>_outlet    [m²/s²]
```

`simpleFoam` uses kinematic pressure. Convert to pascals with
`DeltaP = rho * J_dp`. Legacy result files use the suffix `_Pa`; the value is
still kinematic and is read transparently by the current scripts.

### Intensity of segregation

The outlet function object computes an unweighted face-sample variance:

```text
mean(T) = sum(T_i) / N
I_s = variance(T) / (0.5 * 0.5)
```

The BO minimizes `mixing_intensity_of_segregation`. The same CSV also contains
flux-weighted metrics, but they are diagnostics rather than objectives.

## BO algorithm

```text
8 feasible Sobol designs
          ↓
SingleTaskGP on the negated two-objective response
Normalize inputs + standardize outputs
          ↓
qLogNoisyExpectedHypervolumeImprovement (qNEHVI fallback)
10 restarts, 256 raw samples, q = 1
          ↓
CAD → cfMesh → simpleFoam → scalarTransportFoam → objectives.csv
```

Failed CFD evaluations receive `(J_dp, I_s) = (1000, 1)`.

## Run from any clone location

Source OpenFOAM first; no installation path is hard-coded.

```bash
source /path/to/OpenFOAM/etc/bashrc
cd /path/to/openfoam-bayesopt-mixer
./Allwmake
cd PlanarAlternatingDeflectorMixer
python3 bayes_optimize_sequential.py
```

For a single configured case:

```bash
snakemake --cores 4
```

## Current 28-sample campaign

The stored campaign contains eight Sobol samples and twenty BO suggestions.
All 28 completed without penalty rows.

| Quantity | Best stored value | Sample |
|---|---:|---:|
| kinematic `J_dp` | `1.00487e-3 m²/s²` | `00022` |
| `I_s` | `0.499471` | `00027` |
| mixing quality `1-I_s` | `0.500529` | `00027` |

There are 15 non-dominated samples. Many BO suggestions lie on parameter
bounds (`w_s = 0.45`, `delta = 0.15`, and small baffle thicknesses), so the
front should be treated as a boundary-seeking result rather than a converged
global characterization.

## Interpretation limits

- The device does not demonstrate SAR layer multiplication in the scalar
  fields; it is an alternating-deflector mixer.
- Bounded upwind scalar convection adds numerical diffusion.
- The optimized mixing metric is face-count weighted, not area weighted.
- Pressure labels in pre-rename CSV files are dimensionally wrong.
- A mesh/order/refinement study is needed before reporting physical mixing
  performance.

See `docs/index.html` for the complete Reveal.js presentation.
