# Technical Note 03: Function Object Specification - `outletMixingQuality`

## Purpose
The function object `outletMixingQuality` evaluates concentration-based mixing
quality measures on a user-specified outlet boundary patch of a static mixer.
It is intended for miscible mixing studies where a transported scalar field
(e.g. species mass fraction) is available at the outlet.

## Scope and Quantities Reported
Given a scalar concentration field $a(\mathbf{x},t)$ defined on the outlet patch
$\Gamma_{\mathrm{out}}$, the function object reports, at each write time:

- Mean concentration $\bar{a}$
- (Unweighted) standard deviation $\sigma_a$ and coefficient of variation $\mathrm{CoV}$
- Flux-weighted mean $\bar{a}_F$, flux-weighted variance $\sigma^2_{F,a}$, and flux-weighted $\mathrm{CoV}_F$
- Intensity of segregation $I$
- Relative standard deviation $\mathrm{RSD}$
- Mixing indices $M_{\mathrm{RSD}}$ and $M_I$ (optional)
- Optional absolute deviation measures: $\delta_{\max}$ and $\delta$

## User Inputs (Dictionary Entries)
- `scalarField`: name of scalar concentration field $a$ (e.g. `Y_A`, `alphaA`)
- `patch`: outlet patch name (e.g. `outlet`)
- `weighting`: `none | phi | Un`
- `meanMode`: `fromField | fromInletRatio`
- `aMean`: prescribed mean concentration $\bar{a}$ (used if `meanMode=fromInletRatio`)
- `reportMixingIndex`: boolean
- `reportAbsoluteDeviations`: boolean

### Weighting options
Let $w_i$ denote the weight assigned to outlet face (or sample) $i$:

- `none`: $w_i = 1$
- `phi`: $w_i = \max(\phi_i,0)$ where $\phi$ is the volumetric flux through face $i$
- `Un`: $w_i = \max(U_{n,i} A_i,0)$ where $U_{n,i}=\mathbf{U}_i\cdot \mathbf{n}_i$

Only outflow contributions are included by using $\max(\cdot,0)$.

## Sampling and Discretization
The outlet patch $\Gamma_{\mathrm{out}}$ is discretized by its boundary faces.
For each face $i=1,\dots,N$ on the patch, define:

- Face area $A_i$
- Scalar value $a_i$ (face-interpolated from cell values or native patch values)
- Outflow weight $w_i$ according to the chosen weighting scheme

Define the total weight

$$
W = \sum_{i=1}^N w_i.
$$

## Definitions of Reported Measures
### Unweighted mean and variance
The (unweighted) mean concentration is

$$
\bar{a} = \frac{1}{N}\sum_{i=1}^N a_i.
$$

The (unweighted) variance and standard deviation are

$$
\sigma_a^2 = \frac{1}{N}\sum_{i=1}^N (a_i-\bar{a})^2,
\qquad
\sigma_a = \sqrt{\sigma_a^2}.
$$

### Coefficient of Variation (CoV)
The coefficient of variation is

$$
\mathrm{CoV} = \frac{\sigma_a}{\bar{a}}.
$$

If $\bar{a}$ is computed from the samples, $\bar{a}$ is used directly; if prescribed,
the prescribed $\bar{a}$ is used.

### Flux-weighted mean and variance
The flux-weighted mean concentration is

$$
\bar{a}_F = \frac{1}{W}\sum_{i=1}^N w_i a_i.
$$

The flux-weighted variance is

$$
\sigma^2_{F,a} = \frac{1}{W}\sum_{i=1}^N w_i (a_i-\bar{a}_F)^2,
\qquad
\sigma_{F,a} = \sqrt{\sigma^2_{F,a}}.
$$

The flux-weighted coefficient of variation is

$$
\mathrm{CoV}_F = \frac{\sigma_{F,a}}{\bar{a}_F}.
$$

### Intensity of Segregation
For a binary mixture with mean concentration $\bar{a}$, define the maximum (inlet) variance

$$
\sigma_0^2 = \bar{a}(1-\bar{a}).
$$

Then the intensity of segregation is

$$
I = \frac{\sigma_a^2}{\sigma_0^2}.
$$

If flux-weighting is enabled, the function object also reports

$$
I_F = \frac{\sigma^2_{F,a}}{\sigma_0^2}.
$$

### Relative Standard Deviation (RSD)
Define the inlet coefficient of variation

$$
\mathrm{CoV}_0 = \frac{\sigma_0}{\bar{a}}, \qquad \sigma_0=\sqrt{\bar{a}(1-\bar{a})}.
$$

Then the relative standard deviation is

$$
\mathrm{RSD} = \frac{\mathrm{CoV}}{\mathrm{CoV}_0}
= \frac{\sigma_a}{\sigma_0}
= \sqrt{I}.
$$

Analogously (optional),

$$
\mathrm{RSD}_F = \frac{\sigma_{F,a}}{\sigma_0}=\sqrt{I_F}.
$$

### Mixing Indices
Two common mixing indices are reported when enabled:

$$
M_{\mathrm{RSD}} = 1-\mathrm{RSD},
\qquad
M_I = 1-I.
$$

### Absolute Deviation Measures (Optional)
Let

$$
a_{\min} = \min_i a_i, \qquad a_{\max}=\max_i a_i.
$$

Define the maximum absolute relative deviation

$$
\delta_{\max} = \frac{\max(\,a_{\max}-\bar{a},\,\bar{a}-a_{\min}\,)}{\bar{a}}.
$$

Define the mean absolute relative deviation

$$
\delta = \frac{1}{\bar{a}}\frac{1}{N}\sum_{i=1}^N |a_i-\bar{a}|.
$$

(Optional flux-weighted variants replace $\frac{1}{N}\sum$ by $\frac{1}{W}\sum w_i(\cdot)$.)

## Outputs
At each write time, the function object writes:

- A single line to the log with the scalar metrics.
- A `.dat` file in `postProcessing/outletMixingQuality/` with columns:

$$
t,\ \bar{a},\ \sigma_a,\ \mathrm{CoV},\ \bar{a}_F,\ \sigma_{F,a},\ \mathrm{CoV}_F,\ I,\ I_F,\ \mathrm{RSD},\ \mathrm{RSD}_F,\ M_{\mathrm{RSD}},\ M_I,\ \delta_{\max},\ \delta
$$

(columns included depending on options).

## Numerical Notes
- The reported $\mathrm{CoV}$ depends on sampling resolution (patch faces or optional binning).
- Only outflow faces are included in weighted statistics by truncating negative flux weights.
- If $\bar{a}\to 0$, $\mathrm{CoV}$ becomes ill-conditioned; the implementation should guard divisions by small $\bar{a}$ with a user-configurable $\epsilon$.

## Reference
- Von Damnitz, Lukas, and Denis Anders. "A Review on the Mixing Quality of Static Mixers." *ChemEngineering* 9, no. 6 (2025): 128. https://doi.org/10.3390/chemengineering9060128
