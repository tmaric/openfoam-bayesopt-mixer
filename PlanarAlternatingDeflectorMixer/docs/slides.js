/* ===================================================================
   Minimal slide engine -- the whole of what reveal.js was doing here.
   No CDN, no plugins, no build step: these decks open from a web server
   or from a bare file:// double-click, online or offline.

   Keys:  right / space / page-down  next      left / page-up  previous
          home / end                 first / last
          n                          toggle speaker notes
   =================================================================== */

/* ===================================================================
   6. Figures for the theory deck: the Gaussian, iso-probability
   contours, and samples drawn from a GP prior.  All computed.
   =================================================================== */
function seededRandom(seed) {           /* deterministic, so the slide never changes */
  let a = seed >>> 0;
  return () => {
    a = (a * 1664525 + 1013904223) >>> 0;
    return a / 4294967296;
  };
}
function gaussPair(rnd) {               /* Box-Muller */
  const u = Math.max(rnd(), 1e-12), v = rnd();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function cholesky(A) {
  const n = A.length;
  const L = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j <= i; j += 1) {
      let sum = A[i][j];
      for (let k = 0; k < j; k += 1) sum -= L[i][k] * L[j][k];
      if (i === j) L[i][i] = Math.sqrt(Math.max(sum, 1e-12));
      else L[i][j] = sum / L[j][j];
    }
  }
  return L;
}

function renderGaussian1D(svg) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 760, H = 330, m = { left: 60, right: 40, top: 30, bottom: 56 };
  const mu = 0, sd = 1, lo = -3.6, hi = 3.6;
  const pdf = (x) => Math.exp(-((x - mu) ** 2) / (2 * sd * sd)) / (sd * Math.sqrt(2 * Math.PI));
  const X = (v) => m.left + (v - lo) / (hi - lo) * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - (v / pdf(mu)) * (H - m.top - m.bottom);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  const pts = [];
  for (let i = 0; i <= 300; i += 1) { const x = lo + (i / 300) * (hi - lo); pts.push(`${X(x)},${Y(pdf(x))}`); }
  /* +/- one sigma, shaded */
  const band = [`${X(-sd)},${Y(0)}`];
  for (let i = 0; i <= 120; i += 1) { const x = -sd + (i / 120) * 2 * sd; band.push(`${X(x)},${Y(pdf(x))}`); }
  band.push(`${X(sd)},${Y(0)}`);
  svg.appendChild(svgElement('polygon', { points: band.join(' '), fill: BLUE, opacity: 0.18 }));
  svg.appendChild(svgElement('polyline', { points: pts.join(' '), fill: 'none', stroke: BLUE, 'stroke-width': 4 }));
  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(0), y2: Y(0), stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', {
    x1: X(mu), x2: X(mu), y1: Y(0), y2: Y(pdf(mu)), stroke: INK, 'stroke-width': 2, 'stroke-dasharray': '7 5' }));
  svg.appendChild(svgSym(X(mu), Y(0) + 30, 'μ', '', { 'text-anchor': 'middle', fill: INK, 'font-size': 24, 'font-weight': 700 }));
  [-sd, sd].forEach((v) => {
    svg.appendChild(svgElement('line', { x1: X(v), x2: X(v), y1: Y(0), y2: Y(pdf(v)), stroke: MUTED, 'stroke-width': 2 }));
    svg.appendChild(svgSym(X(v), Y(0) + 30, v < 0 ? 'μ−σ' : 'μ+σ', '',
      { 'text-anchor': 'middle', fill: MUTED, 'font-size': 20 }));
  });
  svg.appendChild(svgText(X(0), Y(pdf(0)) - 14, '68 % of the mass lies within one σ',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 20, 'font-weight': 700 }));
}

function renderIsoContours(svg, mode) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 560, H = 420, m = 54;
  const cx = W / 2, cy = H / 2 - 6, unit = 46;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.appendChild(svgElement('line', { x1: m, x2: W - m, y1: cy, y2: cy, stroke: GRID_C, 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', { x1: cx, x2: cx, y1: m - 20, y2: H - m, stroke: GRID_C, 'stroke-width': 2 }));

  /* dependent case: a linear map tilts and stretches the circles into ellipses */
  const rot = mode === 'dependent' ? -28 : 0;
  const sx = mode === 'dependent' ? 1.55 : 1;
  const sy = mode === 'dependent' ? 0.62 : 1;
  [1, 2, 3].forEach((k, i) => {
    const e = svgElement('ellipse', {
      cx, cy, rx: unit * k * sx, ry: unit * k * sy,
      fill: 'none', stroke: BLUE, 'stroke-width': 3, opacity: 1 - i * 0.24,
    });
    e.setAttribute('transform', `rotate(${rot} ${cx} ${cy})`);
    svg.appendChild(e);
  });
  svg.appendChild(svgElement('circle', { cx, cy, r: 6, fill: INK }));
  svg.appendChild(svgSym(cx + 12, cy - 10, 'μ', '', { fill: INK, 'font-size': 21, 'font-weight': 700 }));
  svg.appendChild(svgSym(W - m + 6, cy + 6, 'y', '1', { fill: MUTED, 'font-size': 20 }));
  svg.appendChild(svgSym(cx + 8, m - 24, 'y', '2', { fill: MUTED, 'font-size': 20 }));
  svg.appendChild(svgText(W / 2, H - 12,
    mode === 'dependent' ? 'off-diagonal Σ: tilted ellipses' : 'diagonal Σ: circles',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 20, 'font-weight': 700 }));
}

function renderGPPrior(svg, lengthScale, narrow = false) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  /* aspect chosen to match the CSS box, otherwise the drawing is letterboxed:
     wide for a full-slide figure, near-square for one cell of a panel row */
  const W = narrow ? 480 : 1450, H = narrow ? 330 : 340;
  const m = narrow ? { left: 20, right: 16, top: 32, bottom: 16 }
                   : { left: 60, right: 50, top: 30, bottom: 34 };
  const n = 70;
  const xs = Array.from({ length: n }, (_, i) => i / (n - 1));
  const k = (a, b) => Math.exp(-((a - b) ** 2) / (2 * lengthScale * lengthScale));
  const K = xs.map((a, i) => xs.map((b, j) => k(a, b) + (i === j ? 1e-8 : 0)));
  const L = cholesky(K);
  const X = (v) => m.left + v * (W - m.left - m.right);
  const Y = (v) => H / 2 - v * ((H - m.top - m.bottom) / 6.2);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  svg.appendChild(svgElement('rect', {
    x: m.left, y: Y(2), width: W - m.left - m.right, height: Y(-2) - Y(2),
    fill: BLUE, opacity: 0.12 }));
  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(0), y2: Y(0), stroke: BLUE, 'stroke-width': 3 }));

  const colours = [INK, ORANGE, GREEN, '#7a2fa8', '#c84b4b'];
  for (let s = 0; s < 5; s += 1) {
    const rnd = seededRandom(1234 + s * 977);
    const z = Array.from({ length: n }, () => gaussPair(rnd));
    const f = L.map((row) => row.reduce((acc, v, j) => acc + v * z[j], 0));
    svg.appendChild(svgElement('polyline', {
      points: xs.map((x, i) => `${X(x)},${Y(f[i])}`).join(' '),
      fill: 'none', stroke: colours[s], 'stroke-width': 2.5, opacity: 0.85,
    }));
  }
  svg.appendChild(svgText(m.left + 4, narrow ? m.top - 12 : m.top + 4,
    narrow ? `ℓ = ${lengthScale}` : `five functions drawn from the prior   ℓ = ${lengthScale}`,
    { fill: narrow ? INK : MUTED, 'font-size': narrow ? 23 : 19, 'font-weight': narrow ? 700 : 400 }));
  if (!narrow) svg.appendChild(svgText(W - m.right - 4, Y(0) - 10, 'prior mean μ = 0',
    { 'text-anchor': 'end', fill: BLUE, 'font-size': 18 }));
}

function renderTheoryFigures() {
  document.querySelectorAll('svg[data-fig]').forEach((svg) => {
    if (svg.dataset.rendered) return;
    svg.dataset.rendered = 'true';
    const kind = svg.dataset.fig;
    if (kind === 'gauss1d') renderGaussian1D(svg);
    else if (kind === 'iso-independent') renderIsoContours(svg, 'independent');
    else if (kind === 'iso-dependent') renderIsoContours(svg, 'dependent');
    else if (kind === 'gp-prior') renderGPPrior(svg, parseFloat(svg.dataset.ell || '0.15'), 'narrow' in svg.dataset);
  });
}


/* ===================================================================
   7. The didactic core of the theory deck.  A reusable GP, then four
   figures that SHOW what the algebra claims:
     - conditioning tightening the posterior as data arrives
     - what the kernel choice actually decides
     - the BO loop running, iteration by iteration
     - BO against random search and a grid, on the same budget
   =================================================================== */
const KERNELS = {
  rbf: (l, sf = 1) => (a, b) => sf * sf * Math.exp(-((a - b) ** 2) / (2 * l * l)),
};

/* Exact GP posterior. Returns a predictor over a grid. */
function fitGP(X, Y, k, noise = 1e-6) {
  const K = X.map((a, i) => X.map((b, j) => k(a, b) + (i === j ? noise : 0)));
  const alpha = X.length ? solve(K, Y) : [];
  return (xs) => xs.map((x) => {
    if (!X.length) return { x, mean: 0, sd: Math.sqrt(k(x, x)) };
    const ks = X.map((xi) => k(x, xi));
    const mean = ks.reduce((acc, v, i) => acc + v * alpha[i], 0);
    const v = solve(K, ks);
    const sd = Math.sqrt(Math.max(k(x, x) - ks.reduce((acc, kk, i) => acc + kk * v[i], 0), 1e-12));
    return { x, mean, sd };
  });
}

const GRID = Array.from({ length: 221 }, (_, i) => i / 220);

/* Shared axes helper for the 1-D figures. */
function axes1d(svg, W, H, m, lo, hi) {
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const X = (v) => m.left + v * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - ((v - lo) / (hi - lo)) * (H - m.top - m.bottom);
  return { X, Y };
}

function drawPosterior(svg, pred, X, Y, opts = {}) {
  const band = pred.map((p) => `${X(p.x)},${Y(p.mean + 2 * p.sd)}`)
    .concat(pred.slice().reverse().map((p) => `${X(p.x)},${Y(p.mean - 2 * p.sd)}`)).join(' ');
  svg.appendChild(svgElement('polygon', { points: band, fill: BLUE, opacity: opts.bandOpacity || 0.16 }));
  svg.appendChild(svgElement('polyline', {
    points: pred.map((p) => `${X(p.x)},${Y(p.mean)}`).join(' '),
    fill: 'none', stroke: BLUE, 'stroke-width': opts.meanWidth || 3.5,
  }));
}

/* ---- A. conditioning: the posterior tightening as data arrives ---------- */
const COND_X = [0.10, 0.34, 0.52, 0.72, 0.90, 0.22];
function renderConditioning(svg, n) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 730, H = 280, m = { left: 26, right: 22, top: 34, bottom: 18 };
  const { X, Y } = axes1d(svg, W, H, m, -2.9, 2.5);
  const k = KERNELS.rbf(0.14, 1.1);
  const xs = COND_X.slice(0, n);
  const ys = xs.map(forrester);
  const pred = fitGP(xs, ys, k, 1e-6)(GRID);

  svg.appendChild(svgElement('polyline', {
    points: GRID.map((x) => `${X(x)},${Y(forrester(x))}`).join(' '),
    fill: 'none', stroke: MUTED, 'stroke-width': 2, 'stroke-dasharray': '7 5',
  }));
  drawPosterior(svg, pred, X, Y);
  xs.forEach((xi, i) => svg.appendChild(svgElement('circle', { cx: X(xi), cy: Y(ys[i]), r: 7, fill: INK })));

  const avgSd = pred.reduce((a, p) => a + p.sd, 0) / pred.length;
  svg.appendChild(svgText(m.left + 4, m.top - 12,
    n === 0 ? 'no data yet' : `${n} observation${n > 1 ? 's' : ''}`,
    { fill: INK, 'font-size': 21, 'font-weight': 700 }));
  svg.appendChild(svgText(W - m.right - 4, m.top - 12, `average σ = ${avgSd.toFixed(2)}`,
    { 'text-anchor': 'end', fill: MUTED, 'font-size': 20 }));
}

/* ---- B. what the kernel decides ---------------------------------------- */
function renderKernel(svg, ell) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 480, H = 330, m = { left: 20, right: 16, top: 32, bottom: 16 };
  /* an over-smooth prior overshoots hard where there is no data, down to about
     -4, so the shared range has to hold it -- clipping it would hide the lesson */
  const { X, Y } = axes1d(svg, W, H, m, -4.5, 2.8);
  const k = KERNELS.rbf(ell, 1.1);
  const xs = COND_X.slice(0, 5), ys = xs.map(forrester);
  const pred = fitGP(xs, ys, k, 1e-6)(GRID);
  svg.appendChild(svgElement('polyline', {
    points: GRID.map((x) => `${X(x)},${Y(forrester(x))}`).join(' '),
    fill: 'none', stroke: MUTED, 'stroke-width': 2, 'stroke-dasharray': '7 5' }));
  drawPosterior(svg, pred, X, Y);
  xs.forEach((xi, i) => svg.appendChild(svgElement('circle', { cx: X(xi), cy: Y(ys[i]), r: 7, fill: INK })));
  const err = Math.max(...pred.map((p) => Math.abs(p.mean - forrester(p.x))));
  svg.appendChild(svgText(m.left + 4, m.top - 12, `ℓ = ${ell.toFixed(2)}`,
    { fill: INK, 'font-size': 23, 'font-weight': 700 }));
  svg.appendChild(svgText(W - m.right - 4, m.top - 12, `worst error ${err.toFixed(2)}`,
    { 'text-anchor': 'end', fill: MUTED, 'font-size': 20 }));
}

/* ---- C. the loop running ------------------------------------------------ */
function boRun(steps, kappa = 2.0) {
  const k = KERNELS.rbf(0.14, 1.1);
  const xs = [0.10, 0.52, 0.90];
  const ys = xs.map(forrester);
  const picks = [];
  for (let t = 0; t < steps; t += 1) {
    const pred = fitGP(xs, ys, k, 1e-6)(GRID);
    const best = pred.reduce((a, p) => ((p.mean + kappa * p.sd) > (a.mean + kappa * a.sd) ? p : a));
    picks.push(best.x);
    xs.push(best.x); ys.push(forrester(best.x));
  }
  return { xs, ys, picks, k, kappa };
}
function renderBOStep(svg, step) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 730, H = 280, m = { left: 26, right: 22, top: 36, bottom: 18 };
  const { xs, ys, k, kappa } = boRun(step + 1);
  const shown = xs.slice(0, 3 + step), shownY = ys.slice(0, 3 + step);
  const pred = fitGP(shown, shownY, k, 1e-6)(GRID);
  const acq = pred.map((p) => p.mean + kappa * p.sd);
  const lo = -2.9, hi = Math.max(3.0, Math.max(...acq));
  const { X, Y } = axes1d(svg, W, H, m, lo, hi);

  svg.appendChild(svgElement('polyline', {
    points: GRID.map((x) => `${X(x)},${Y(forrester(x))}`).join(' '),
    fill: 'none', stroke: MUTED, 'stroke-width': 2, 'stroke-dasharray': '7 5' }));
  drawPosterior(svg, pred, X, Y);
  svg.appendChild(svgElement('polyline', {
    points: pred.map((p, i) => `${X(p.x)},${Y(acq[i])}`).join(' '),
    fill: 'none', stroke: ORANGE, 'stroke-width': 3 }));
  shown.forEach((xi, i) => svg.appendChild(svgElement('circle', {
    cx: X(xi), cy: Y(shownY[i]), r: 7, fill: INK })));

  const nextX = xs[3 + step];
  if (nextX !== undefined) {
    svg.appendChild(svgElement('line', {
      x1: X(nextX), x2: X(nextX), y1: m.top, y2: H - m.bottom,
      stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '8 5' }));
    svg.appendChild(svgElement('circle', {
      cx: X(nextX), cy: Y(forrester(nextX)), r: 9,
      fill: 'none', stroke: ORANGE, 'stroke-width': 3.5 }));
  }
  const bestSoFar = Math.max(...shownY);
  svg.appendChild(svgText(m.left + 4, m.top - 12,
    `${shown.length} simulations spent`, { fill: INK, 'font-size': 21, 'font-weight': 700 }));
  svg.appendChild(svgText(W - m.right - 4, m.top - 12,
    `best so far  f = ${bestSoFar.toFixed(2)}  of  1.00`,
    { 'text-anchor': 'end', fill: MUTED, 'font-size': 20 }));
}

/* ---- D. BO vs random vs grid, same budget ------------------------------- */
function efficiencyCurves(budget = 18) {
  const k = KERNELS.rbf(0.14, 1.1);
  const seed = [0.10, 0.52, 0.90];
  const runBO = (kappa) => {
    const xs = seed.slice(), ys = xs.map(forrester);
    const curve = [Math.max(...ys)];
    for (let t = 3; t < budget; t += 1) {
      const pred = fitGP(xs, ys, k, 1e-6)(GRID);
      const best = pred.reduce((a, p) => ((p.mean + kappa * p.sd) > (a.mean + kappa * a.sd) ? p : a));
      xs.push(best.x); ys.push(forrester(best.x));
      curve.push(Math.max(...ys));
    }
    return curve;
  };
  /* random search averaged over 40 seeds, so the curve is not one lucky draw */
  const rand = new Array(budget - 2).fill(0);
  const R = 40;
  for (let r = 0; r < R; r += 1) {
    const rnd = seededRandom(97 + r * 31);
    let best = -Infinity; const row = [];
    for (let t = 0; t < budget; t += 1) {
      best = Math.max(best, forrester(rnd()));
      if (t >= 2) row.push(best);
    }
    row.forEach((v, i) => { rand[i] += v / R; });
  }
  return { bo: runBO(2.0), greedy: runBO(0.0), rand, budget };
}

function renderEfficiency(svg) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const { bo, greedy, rand, budget } = efficiencyCurves();
  const W = 1400, H = 470, m = { left: 96, right: 392, top: 30, bottom: 74 };
  const lo = 0, hi = 1.12;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const X = (n) => m.left + ((n - 3) / (budget - 3)) * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - ((v - lo) / (hi - lo)) * (H - m.top - m.bottom);

  [0, 0.25, 0.5, 0.75, 1].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(t), y2: Y(t), stroke: GRID_C }));
    svg.appendChild(svgText(m.left - 14, Y(t) + 7, t.toFixed(2), { 'text-anchor': 'end', 'font-size': 19, fill: MUTED }));
  });
  svg.appendChild(svgElement('line', {
    x1: m.left, x2: W - m.right, y1: Y(1.0032), y2: Y(1.0032),
    stroke: GREEN, 'stroke-width': 2.5, 'stroke-dasharray': '9 6' }));
  svg.appendChild(svgText(W - m.right - 8, Y(1.0032) - 12, 'true optimum',
    { 'text-anchor': 'end', fill: GREEN, 'font-size': 19, 'font-weight': 700 }));

  const series = [
    [bo, ORANGE, 'Bayesian optimization, κ = 2'],
    [rand, BLUE, 'random search (mean of 40 runs)'],
    [greedy, RED, 'Bayesian optimization, κ = 0'],
  ];
  /* greedy is dashed: for the first segment it lies exactly under the kappa = 2
     run, and a solid line there would simply hide one of the two */
  series.forEach(([arr, colour], i) => {
    svg.appendChild(svgElement('polyline', Object.assign({
      points: arr.map((v, j) => `${X(j + 3)},${Y(v)}`).join(' '),
      fill: 'none', stroke: colour, 'stroke-width': 4.5,
      'stroke-linejoin': 'round', 'stroke-linecap': 'round',
    }, i === 2 ? { 'stroke-dasharray': '11 7' } : {})));
  });

  /* the two moments worth naming */
  const hit = bo.findIndex((v) => v > 1.0);
  if (hit >= 0) {
    svg.appendChild(svgElement('circle', { cx: X(hit + 3), cy: Y(bo[hit]), r: 9, fill: ORANGE }));
    svg.appendChild(svgText(X(hit + 3) + 16, Y(bo[hit]) + 30,
      `optimum found on evaluation ${hit + 3}`, { fill: ORANGE, 'font-size': 20, 'font-weight': 700 }));
  }
  svg.appendChild(svgText(X(budget) - 6, Y(greedy[greedy.length - 1]) - 16,
    `stalls at ${greedy[greedy.length - 1].toFixed(2)}`,
    { 'text-anchor': 'end', fill: RED, 'font-size': 20, 'font-weight': 700 }));

  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(lo), y2: Y(lo), stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', { x1: m.left, x2: m.left, y1: m.top, y2: Y(lo), stroke: INK, 'stroke-width': 2 }));
  [3, 6, 9, 12, 15, 18].forEach((n) => svg.appendChild(
    svgText(X(n), H - m.bottom + 28, String(n), { 'text-anchor': 'middle', 'font-size': 19, fill: MUTED })));
  svg.appendChild(svgText((m.left + W - m.right) / 2, H - 16, 'expensive evaluations spent',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 22 }));
  const yl = svgText(32, H / 2, 'best value found so far', { 'text-anchor': 'middle', fill: INK, 'font-size': 22 });
  yl.setAttribute('transform', `rotate(-90 32 ${H / 2})`);
  svg.appendChild(yl);

  series.forEach(([, colour, label], i) => {
    const ly = m.top + 30 + i * 38;
    svg.appendChild(svgElement('rect', Object.assign(
      { x: W - m.right + 16, y: ly - 14, width: 30, height: 7, fill: colour, rx: 3 },
      i === 2 ? { width: 12 } : {})));
    if (i === 2) svg.appendChild(svgElement('rect', { x: W - m.right + 34, y: ly - 14, width: 12, height: 7, fill: colour, rx: 3 }));
    svg.appendChild(svgText(W - m.right + 56, ly, label, { 'font-size': 20, fill: INK }));
  });
}

/* ---- E. maximum likelihood: what it fits, and what it does NOT fix -------
   Ten samples, not five.  With five the profile likelihood is nearly flat and
   the optimiser runs to the lower bound on the length scale -- true, and worth
   saying out loud, but it does not make the figure the slide needs. */
const MLE_X = Array.from({ length: 10 }, (_, i) => (i + 0.5) / 10);
const MLE_Y = MLE_X.map(forrester);
const MLE_NOISE = 0.02;

function logMarginalLikelihood(X, Y, l, sf, sn) {
  const n = X.length;
  const K = X.map((a, i) => X.map((b, j) =>
    sf * sf * Math.exp(-((a - b) ** 2) / (2 * l * l)) + (i === j ? sn * sn : 0)));
  const L = cholesky(K);
  const alpha = solve(K, Y);
  let quad = 0, logdet = 0;
  for (let i = 0; i < n; i += 1) { quad += Y[i] * alpha[i]; logdet += 2 * Math.log(L[i][i]); }
  return -0.5 * quad - 0.5 * logdet - 0.5 * n * Math.log(2 * Math.PI);
}

/* profile likelihood: for each length scale, the best signal variance */
function mleProfile() {
  const out = [];
  for (let i = 0; i <= 110; i += 1) {
    const l = 0.02 + i * 0.003;
    let best = -Infinity, sf = 0;
    for (let j = 0; j < 300; j += 1) {
      const s = 0.1 + j * 0.02;
      const v = logMarginalLikelihood(MLE_X, MLE_Y, l, s, MLE_NOISE);
      if (v > best) { best = v; sf = s; }
    }
    out.push({ l, logL: best, sf });
  }
  const star = out.reduce((a, b) => (b.logL > a.logL ? b : a));
  return { out, star };
}
let MLE_CACHE = null;
const mle = () => (MLE_CACHE || (MLE_CACHE = mleProfile()));

function renderMLEProfile(svg) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const { out, star } = mle();
  const W = 720, H = 330, m = { left: 76, right: 26, top: 34, bottom: 54 };
  const lo = -15, hi = -5.4;
  const shown = out.filter((p) => p.logL >= lo - 1);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const X = (v) => m.left + ((v - 0.02) / 0.36) * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - ((v - lo) / (hi - lo)) * (H - m.top - m.bottom);
  [-14, -12, -10, -8, -6].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(t), y2: Y(t), stroke: GRID_C }));
    svg.appendChild(svgText(m.left - 12, Y(t) + 6, String(t), { 'text-anchor': 'end', fill: MUTED, 'font-size': 18 }));
  });
  svg.appendChild(svgElement('polyline', {
    points: shown.map((p) => `${X(p.l)},${Y(Math.max(p.logL, lo))}`).join(' '),
    fill: 'none', stroke: BLUE, 'stroke-width': 4 }));
  svg.appendChild(svgElement('line', {
    x1: X(star.l), x2: X(star.l), y1: Y(star.logL), y2: H - m.bottom,
    stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '8 5' }));
  svg.appendChild(svgElement('circle', { cx: X(star.l), cy: Y(star.logL), r: 8, fill: ORANGE }));
  svg.appendChild(svgText(X(star.l) + 14, Y(star.logL) + 6,
    `ℓ* = ${star.l.toFixed(2)},  σ_f* = ${star.sf.toFixed(2)}`.replace('σ_f', 'σf'),
    { fill: ORANGE, 'font-size': 21, 'font-weight': 700 }));
  [0.05, 0.1, 0.2, 0.3].forEach((t) => svg.appendChild(
    svgText(X(t), H - m.bottom + 26, t.toFixed(2), { 'text-anchor': 'middle', fill: MUTED, 'font-size': 18 })));
  svg.appendChild(svgText((m.left + W - m.right) / 2, H - 12, 'length scale ℓ',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 21 }));
  const yl = svgText(24, H / 2, 'log p(y | θ)', { 'text-anchor': 'middle', fill: INK, 'font-size': 21 });
  yl.setAttribute('transform', `rotate(-90 24 ${H / 2})`);
  svg.appendChild(yl);
  svg.appendChild(svgText(m.left + 6, m.top - 12, 'the fit is a one-dimensional hill-climb',
    { fill: MUTED, 'font-size': 19 }));
}

function renderMLEPrior(svg) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const { star } = mle();
  const W = 720, H = 330, m = { left: 26, right: 22, top: 34, bottom: 20 };
  const { X, Y } = axes1d(svg, W, H, m, -4.2, 4.2);
  const n = 90;
  const xs = Array.from({ length: n }, (_, i) => i / (n - 1));
  const k = (a, b) => star.sf * star.sf * Math.exp(-((a - b) ** 2) / (2 * star.l * star.l));
  const L = cholesky(xs.map((a, i) => xs.map((b, j) => k(a, b) + (i === j ? 1e-8 : 0))));
  svg.appendChild(svgElement('rect', {
    x: m.left, y: Y(2 * star.sf), width: W - m.left - m.right,
    height: Y(-2 * star.sf) - Y(2 * star.sf), fill: BLUE, opacity: 0.12 }));
  svg.appendChild(svgElement('line', {
    x1: m.left, x2: W - m.right, y1: Y(0), y2: Y(0), stroke: BLUE, 'stroke-width': 3.5 }));
  const colours = [MUTED, ORANGE, GREEN, '#7a2fa8', RED];
  for (let sIdx = 0; sIdx < 5; sIdx += 1) {
    const rnd = seededRandom(4242 + sIdx * 613);
    const z = Array.from({ length: n }, () => gaussPair(rnd));
    const f = L.map((row) => row.reduce((acc, v, j) => acc + v * z[j], 0));
    svg.appendChild(svgElement('polyline', {
      points: xs.map((x, i) => `${X(x)},${Y(f[i])}`).join(' '),
      fill: 'none', stroke: colours[sIdx], 'stroke-width': 2.2, opacity: 0.7 }));
  }
  MLE_X.forEach((xi, i) => svg.appendChild(svgElement('circle', {
    cx: X(xi), cy: Y(MLE_Y[i]), r: 7, fill: INK })));
  svg.appendChild(svgText(m.left + 6, m.top - 12,
    `prior with ℓ* = ${star.l.toFixed(2)}, σf* = ${star.sf.toFixed(2)}`,
    { fill: INK, 'font-size': 21, 'font-weight': 700 }));
  svg.appendChild(svgText(W - m.right - 4, H - m.bottom - 8, 'mean still 0 — no data used yet',
    { 'text-anchor': 'end', fill: BLUE, 'font-size': 19 }));
}

/* ---- F. marginalise versus condition, on one joint Gaussian --------------
   The two operations people conflate.  Same joint, same picture: marginalising
   projects all the mass onto one axis, conditioning takes a single slice
   through it and renormalises.  rho = 0.8 so the difference is unmissable. */
const MC_RHO = 0.8;
const MC_C = 1.5;                       /* the value y2 is observed to take */

function renderMarginalConditional(svg, mode) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 720, H = 380, m = { left: 66, right: 30, top: 118, bottom: 56 };
  const lim = 3.2;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const X = (v) => m.left + ((v + lim) / (2 * lim)) * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - ((v + lim) / (2 * lim)) * (H - m.top - m.bottom);
  /* the density strip lives above the panel, drawn downwards from its baseline */
  const base = m.top - 14;
  const D = (d) => base - d * 84;

  [-3, -2, -1, 0, 1, 2, 3].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: X(t), x2: X(t), y1: m.top, y2: H - m.bottom, stroke: GRID_C }));
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(t), y2: Y(t), stroke: GRID_C }));
    svg.appendChild(svgText(X(t), H - m.bottom + 24, String(t), { 'text-anchor': 'middle', fill: MUTED, 'font-size': 17 }));
    svg.appendChild(svgText(m.left - 12, Y(t) + 6, String(t), { 'text-anchor': 'end', fill: MUTED, 'font-size': 17 }));
  });

  /* iso-probability ellipses of N(0, [[1,rho],[rho,1]]) via its Cholesky factor */
  const a = Math.sqrt(1 - MC_RHO * MC_RHO);
  [1, 2, 3].forEach((r) => {
    const pts = [];
    for (let i = 0; i <= 120; i += 1) {
      const t = (i / 120) * 2 * Math.PI;
      pts.push(`${X(r * Math.cos(t))},${Y(MC_RHO * r * Math.cos(t) + a * r * Math.sin(t))}`);
    }
    svg.appendChild(svgElement('polygon', {
      points: pts.join(' '), fill: BLUE, opacity: 0.10, stroke: BLUE,
      'stroke-width': 1.6, 'stroke-opacity': 0.5 }));
  });
  svg.appendChild(svgText(W - m.right - 6, H - m.bottom + 24, 'y₁',
    { 'text-anchor': 'end', fill: INK, 'font-size': 21, 'font-weight': 700 }));
  const y2lab = svgText(22, (m.top + H - m.bottom) / 2, 'y₂',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 21, 'font-weight': 700 });
  y2lab.setAttribute('transform', `rotate(-90 22 ${(m.top + H - m.bottom) / 2})`);
  svg.appendChild(y2lab);

  const gauss = (x, mu, sd) => Math.exp(-((x - mu) ** 2) / (2 * sd * sd)) / (sd * Math.sqrt(2 * Math.PI));
  const curve = (mu, sd, colour) => {
    const pts = [];
    for (let i = 0; i <= 200; i += 1) {
      const x = -lim + (i / 200) * 2 * lim;
      pts.push(`${X(x)},${D(gauss(x, mu, sd))}`);
    }
    svg.appendChild(svgElement('polyline', { points: `${X(-lim)},${D(0)} ` + pts.join(' ') + ` ${X(lim)},${D(0)}`,
      fill: colour, opacity: 0.16, stroke: 'none' }));
    svg.appendChild(svgElement('polyline', { points: pts.join(' '), fill: 'none', stroke: colour, 'stroke-width': 3.5 }));
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: D(0), y2: D(0), stroke: MUTED, 'stroke-width': 1.5 }));
  };

  if (mode === 'marginal') {
    /* every column of the joint, summed: N(0, 1) */
    for (let i = 0; i <= 26; i += 1) {
      const x = -lim + (i / 26) * 2 * lim;
      svg.appendChild(svgElement('line', {
        x1: X(x), x2: X(x), y1: m.top, y2: H - m.bottom,
        stroke: ORANGE, 'stroke-width': 1.4, opacity: 0.35 }));
    }
    curve(0, 1, ORANGE);
    svg.appendChild(svgText(m.left, 26, 'p(y₁) = ∫ p(y₁, y₂) dy₂',
      { fill: ORANGE, 'font-size': 22, 'font-weight': 700 }));
    svg.appendChild(svgText(m.left, 50, 'mean 0,  sd 1.00  —  as wide as it started',
      { fill: MUTED, 'font-size': 19 }));
  } else {
    const mu = MC_RHO * MC_C, sd = Math.sqrt(1 - MC_RHO * MC_RHO);
    svg.appendChild(svgElement('line', {
      x1: m.left, x2: W - m.right, y1: Y(MC_C), y2: Y(MC_C),
      stroke: ORANGE, 'stroke-width': 3.5 }));
    svg.appendChild(svgText(W - m.right - 6, Y(MC_C) - 10, `y₂ = ${MC_C}  observed`,
      { 'text-anchor': 'end', fill: ORANGE, 'font-size': 19, 'font-weight': 700 }));
    curve(mu, sd, ORANGE);
    svg.appendChild(svgElement('line', {
      x1: X(mu), x2: X(mu), y1: D(0), y2: D(gauss(mu, mu, sd)),
      stroke: ORANGE, 'stroke-width': 2, 'stroke-dasharray': '6 4' }));
    svg.appendChild(svgText(m.left, 26, 'p(y₁ | y₂ = 1.5)',
      { fill: ORANGE, 'font-size': 22, 'font-weight': 700 }));
    svg.appendChild(svgText(m.left, 50, `mean ${mu.toFixed(2)},  sd ${sd.toFixed(2)}  —  moved, and 40% as wide`,
      { fill: MUTED, 'font-size': 19 }));
  }
}

/* ---- G. why a weighted sum cannot reach a non-convex front ---------------
   Both objectives minimised.  A weighted sum picks whichever front point its
   iso-cost line touches first, so it can only ever select vertices of the
   LOWER CONVEX HULL.  Sweep the weights across a concave stretch and the
   chosen design jumps over it -- every design in between is unselectable, for
   every weighting.  Front, hull and jump are all computed here, not drawn. */
const SC_FRONT = (t) => 1 - Math.sqrt(Math.max(2 * t - t * t, 0))
                          + 0.14 * Math.exp(-(((t - 0.5) / 0.12) ** 2));

function renderScalarise(svg) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 1340, H = 470, m = { left: 96, right: 300, top: 34, bottom: 74 };
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const X = (v) => m.left + (v / 1.06) * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - (v / 1.06) * (H - m.top - m.bottom);

  const N = 400;
  const P = Array.from({ length: N + 1 }, (_, i) => [i / N, SC_FRONT(i / N)]);
  /* lower convex hull, monotone chain */
  const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const hull = [];
  P.forEach((pt) => {
    while (hull.length >= 2 && cross(hull[hull.length - 2], hull[hull.length - 1], pt) <= 0) hull.pop();
    hull.push(pt);
  });
  /* the one hull edge that bridges a gap is the concave stretch */
  let bridge = [hull[0], hull[1]];
  for (let i = 0; i < hull.length - 1; i += 1) {
    if (hull[i + 1][0] - hull[i][0] > bridge[1][0] - bridge[0][0]) bridge = [hull[i], hull[i + 1]];
  }
  const [A, B] = bridge;                          /* A is the left end, B the right */
  const inside = (pt) => pt[0] > A[0] + 1e-9 && pt[0] < B[0] - 1e-9;

  [0, 0.25, 0.5, 0.75, 1].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(t), y2: Y(t), stroke: GRID_C }));
    svg.appendChild(svgElement('line', { x1: X(t), x2: X(t), y1: m.top, y2: H - m.bottom, stroke: GRID_C }));
    svg.appendChild(svgText(m.left - 14, Y(t) + 6, t.toFixed(2), { 'text-anchor': 'end', fill: MUTED, 'font-size': 18 }));
    svg.appendChild(svgText(X(t), H - m.bottom + 26, t.toFixed(2), { 'text-anchor': 'middle', fill: MUTED, 'font-size': 18 }));
  });

  /* the weighted-sum iso-cost line through both ends of the bridge */
  const w1 = (B[1] - A[1]) / ((B[1] - A[1]) + (A[0] - B[0]));
  const cost = w1 * A[0] + (1 - w1) * A[1];
  const isoAt = (x) => (cost - w1 * x) / (1 - w1);
  const isoLo = Math.max(0, (cost - (1 - w1) * 1.06) / w1);
  const isoHi = Math.min(1.06, cost / w1);
  svg.appendChild(svgElement('line', {
    x1: X(isoLo), y1: Y(isoAt(isoLo)), x2: X(isoHi), y2: Y(isoAt(isoHi)),
    stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '10 6' }));

  /* the front: reachable in blue, unreachable in red */
  const seg = (pts, colour, width) => svg.appendChild(svgElement('polyline', {
    points: pts.map((q) => `${X(q[0])},${Y(q[1])}`).join(' '),
    fill: 'none', stroke: colour, 'stroke-width': width, 'stroke-linecap': 'round' }));
  seg(P.filter((q) => q[0] <= A[0]), BLUE, 5);
  seg(P.filter((q) => q[0] >= B[0]), BLUE, 5);
  seg(P.filter(inside), RED, 6);

  [[A, 'A'], [B, 'B']].forEach(([pt, lab]) => {
    svg.appendChild(svgElement('circle', { cx: X(pt[0]), cy: Y(pt[1]), r: 9, fill: ORANGE }));
    svg.appendChild(svgText(X(pt[0]) + (lab === 'A' ? -18 : 16), Y(pt[1]) - 14, lab,
      { 'text-anchor': lab === 'A' ? 'end' : 'start', fill: ORANGE, 'font-size': 24, 'font-weight': 700 }));
  });
  const mid = P.find((q) => q[0] > (A[0] + B[0]) / 2);
  svg.appendChild(svgText(X(mid[0]), Y(mid[1]) - 22, 'never selected, for any weights',
    { 'text-anchor': 'middle', fill: RED, 'font-size': 21, 'font-weight': 700 }));

  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: Y(0), y2: Y(0), stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', { x1: m.left, x2: m.left, y1: m.top, y2: Y(0), stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgText((m.left + W - m.right) / 2, H - 14, 'objective 1  (minimise →)',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 21 }));
  const yl = svgText(30, H / 2, 'objective 2  (minimise ↓)', { 'text-anchor': 'middle', fill: INK, 'font-size': 21 });
  yl.setAttribute('transform', `rotate(-90 30 ${H / 2})`);
  svg.appendChild(yl);

  const L = W - m.right + 18;
  const lines = [
    ['The weighted sum', INK, 23, 700],
    [`minimise  w₁f₁ + w₂f₂`, MUTED, 21, 400],
    ['', INK, 10, 400],
    [`Its iso-cost lines are straight, so it`, MUTED, 20, 400],
    ['can only ever touch the front at a', MUTED, 20, 400],
    ['vertex of the convex hull.', MUTED, 20, 400],
    ['', INK, 10, 400],
    [`Raise w₁ past ${w1.toFixed(2)} and the answer`, INK, 20, 640],
    [`jumps B → A, from (${B[0].toFixed(2)}, ${B[1].toFixed(2)})`, INK, 20, 640],
    [`to (${A[0].toFixed(2)}, ${A[1].toFixed(2)}) — with nothing`, INK, 20, 640],
    ['in between ever chosen.', INK, 20, 640],
  ];
  let y = m.top + 26;
  lines.forEach(([txt, fill, size, weight]) => {
    if (txt) svg.appendChild(svgText(L, y, txt, { fill, 'font-size': size, 'font-weight': weight }));
    y += size + 8;
  });
}

/* ---- H. what a launch failure costs when it is recorded as physics -------
   Same six designs, same kernel.  On the right, the two nearest the optimum
   did not fail physically -- the tool failed -- but they were written back as
   the penalty value every failure handler writes.  The GP believes the best
   region is the worst one, and the acquisition leaves and does not return. */
const FAIL_X = [0.05, 0.25, 0.45, 0.72, 0.80, 0.95];
const FAIL_IDX = [3, 4];
const FAIL_PENALTY = -2.5;

function renderFailureGP(svg, mode) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 720, H = 330, m = { left: 26, right: 22, top: 36, bottom: 20 };
  const { X, Y } = axes1d(svg, W, H, m, -3.6, 3.0);
  const poisoned = mode === 'poisoned';
  const ys = FAIL_X.map((x, i) => (poisoned && FAIL_IDX.includes(i) ? FAIL_PENALTY : forrester(x)));
  const pred = fitGP(FAIL_X, ys, KERNELS.rbf(0.14, 1.1), 1e-6)(GRID);
  const acq = pred.map((q) => q.mean + 2 * q.sd);
  const pick = pred[acq.indexOf(Math.max(...acq))];

  svg.appendChild(svgElement('polyline', {
    points: GRID.map((x) => `${X(x)},${Y(forrester(x))}`).join(' '),
    fill: 'none', stroke: MUTED, 'stroke-width': 2, 'stroke-dasharray': '7 5' }));
  drawPosterior(svg, pred, X, Y);
  svg.appendChild(svgElement('polyline', {
    points: pred.map((q, i) => `${X(q.x)},${Y(acq[i])}`).join(' '),
    fill: 'none', stroke: ORANGE, 'stroke-width': 2.5, opacity: 0.9 }));
  svg.appendChild(svgElement('line', {
    x1: X(pick.x), x2: X(pick.x), y1: m.top, y2: H - m.bottom,
    stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '8 5' }));

  /* the true optimum, so the eye has something to compare the pick against */
  svg.appendChild(svgElement('line', {
    x1: X(0.759), x2: X(0.759), y1: Y(1.0), y2: Y(-3.4),
    stroke: GREEN, 'stroke-width': 2, 'stroke-dasharray': '4 4', opacity: 0.8 }));
  svg.appendChild(svgText(X(0.759), Y(-3.4) + 2, 'true optimum',
    { 'text-anchor': 'middle', fill: GREEN, 'font-size': 17 }));

  FAIL_X.forEach((xi, i) => {
    const bad = poisoned && FAIL_IDX.includes(i);
    svg.appendChild(svgElement('circle', {
      cx: X(xi), cy: Y(ys[i]), r: bad ? 8 : 7, fill: bad ? RED : INK }));
    if (bad) {
      svg.appendChild(svgElement('line', { x1: X(xi) - 12, x2: X(xi) + 12, y1: Y(ys[i]) - 12, y2: Y(ys[i]) + 12, stroke: RED, 'stroke-width': 3 }));
      svg.appendChild(svgElement('line', { x1: X(xi) - 12, x2: X(xi) + 12, y1: Y(ys[i]) + 12, y2: Y(ys[i]) - 12, stroke: RED, 'stroke-width': 3 }));
    }
  });
  svg.appendChild(svgText(m.left + 4, m.top - 12,
    poisoned ? 'two recorded as failures' : 'all recorded correctly',
    { fill: poisoned ? RED : INK, 'font-size': 21, 'font-weight': 700 }));
  svg.appendChild(svgText(W - m.right - 4, m.top - 12, `next pick  x = ${pick.x.toFixed(2)}`,
    { 'text-anchor': 'end', fill: ORANGE, 'font-size': 20, 'font-weight': 700 }));
}

/* ---- I. two fidelities of one simulation ---------------------------------
   The standard multi-fidelity Forrester pair, in this deck's scaling.  The
   cheap model is BIASED, not noisy: its own maximum is in the wrong place.
   Co-kriging is the recursive form -- fit the cheap model densely, regress the
   expensive points onto it, and give the leftover discrepancy its own GP. */
const MF_LO = (x) => 0.5 * forrester(x) - (10 * x - 5) / 6 + 5 / 6;
const MF_XH = [0, 0.35, 0.65, 1];
const MF_XL = Array.from({ length: 11 }, (_, i) => i / 10);

function mfModels() {
  const YH = MF_XH.map(forrester);
  const hfOnly = fitGP(MF_XH, YH, KERNELS.rbf(0.20, 1.4), 1e-8);
  const lo = fitGP(MF_XL, MF_XL.map(MF_LO), KERNELS.rbf(0.12, 1.4), 1e-8);
  const ml = lo(MF_XH).map((q) => q.mean);
  const n = MF_XH.length;
  const sx = ml.reduce((a, b) => a + b, 0), sy = YH.reduce((a, b) => a + b, 0);
  const sxx = ml.reduce((a, b) => a + b * b, 0);
  const sxy = ml.reduce((a, b, i) => a + b * YH[i], 0);
  const rho = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const beta = (sy - rho * sx) / n;
  const dg = fitGP(MF_XH, MF_XH.map((x, i) => YH[i] - rho * ml[i] - beta),
                   KERNELS.rbf(0.40, 1.1), 1e-8);
  const co = (xs) => {
    const a = lo(xs), b = dg(xs);
    return xs.map((x, i) => ({ x, mean: rho * a[i].mean + beta + b[i].mean }));
  };
  return { hfOnly, lo, co, rho, beta };
}
let MF_CACHE = null;
const mf = () => (MF_CACHE || (MF_CACHE = mfModels()));

function renderMultiFidelity(svg, mode) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 720, H = 330, m = { left: 26, right: 22, top: 36, bottom: 20 };
  const rmse = (v) => Math.sqrt(GRID.reduce((a, x, i) => a + (v[i] - forrester(x)) ** 2, 0) / GRID.length);
  const line = (v, colour, width, dash) => svg.appendChild(svgElement('polyline', Object.assign({
    points: GRID.map((x, i) => `${X(x)},${Y(v[i])}`).join(' '),
    fill: 'none', stroke: colour, 'stroke-width': width }, dash ? { 'stroke-dasharray': dash } : {})));
  let X, Y;

  if (mode === 'functions') {
    ({ X, Y } = axes1d(svg, W, H, m, -3.0, 2.2));
    const hi = GRID.map(forrester), lo = GRID.map(MF_LO);
    line(hi, BLUE, 4);
    line(lo, ORANGE, 4);
    [[hi, BLUE, 'high'], [lo, ORANGE, 'low']].forEach(([v, c]) => {
      const i = v.indexOf(Math.max(...v));
      svg.appendChild(svgElement('circle', { cx: X(GRID[i]), cy: Y(v[i]), r: 8, fill: c }));
      const high = Y(v[i]) < m.top + 46;
      svg.appendChild(svgText(X(GRID[i]) + (i < 20 ? 16 : 0), Y(v[i]) + (high ? 28 : -16),
        `x = ${GRID[i].toFixed(2)}`,
        { 'text-anchor': i < 20 ? 'start' : 'middle', fill: c, 'font-size': 19, 'font-weight': 700 }));
    });
    svg.appendChild(svgText(m.left + 4, m.top - 12, 'high fidelity', { fill: BLUE, 'font-size': 21, 'font-weight': 700 }));
    svg.appendChild(svgText(m.left + 168, m.top - 12, 'low fidelity', { fill: ORANGE, 'font-size': 21, 'font-weight': 700 }));
    svg.appendChild(svgText(W - m.right - 4, m.top - 12, 'correlation 0.74', { 'text-anchor': 'end', fill: MUTED, 'font-size': 19 }));
  } else {
    const { hfOnly, co } = mf();
    const a = hfOnly(GRID).map((q) => q.mean), b = co(GRID).map((q) => q.mean);
    ({ X, Y } = axes1d(svg, W, H, m, -3.0, 2.2));
    line(GRID.map(forrester), MUTED, 2, '7 5');
    line(a, RED, 4);
    line(b, BLUE, 4);
    MF_XH.forEach((xi) => svg.appendChild(svgElement('circle', { cx: X(xi), cy: Y(forrester(xi)), r: 7, fill: INK })));
    svg.appendChild(svgText(m.left + 4, m.top - 12, `4 expensive runs only — RMSE ${rmse(a).toFixed(2)}`,
      { fill: RED, 'font-size': 20, 'font-weight': 700 }));
    svg.appendChild(svgText(W - m.right - 4, m.top - 12, `+ 11 cheap ones — RMSE ${rmse(b).toFixed(2)}`,
      { 'text-anchor': 'end', fill: BLUE, 'font-size': 20, 'font-weight': 700 }));
  }
}

function renderDidacticFigures() {
  document.querySelectorAll('svg[data-cond]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderConditioning(svg, parseInt(svg.dataset.cond, 10));
  });
  document.querySelectorAll('svg[data-kernel]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderKernel(svg, parseFloat(svg.dataset.kernel));
  });
  document.querySelectorAll('svg[data-bostep]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderBOStep(svg, parseInt(svg.dataset.bostep, 10));
  });
  document.querySelectorAll('svg[data-failgp]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderFailureGP(svg, svg.dataset.failgp);
  });
  document.querySelectorAll('svg[data-mf]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderMultiFidelity(svg, svg.dataset.mf);
  });
  document.querySelectorAll('svg[data-scalarise]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderScalarise(svg);
  });
  document.querySelectorAll('svg[data-mc]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderMarginalConditional(svg, svg.dataset.mc);
  });
  document.querySelectorAll('svg[data-mle]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    if (svg.dataset.mle === 'profile') renderMLEProfile(svg); else renderMLEPrior(svg);
  });
  document.querySelectorAll('svg[data-efficiency]').forEach((svg) => {
    if (svg.dataset.rendered) return; svg.dataset.rendered = 'true';
    renderEfficiency(svg);
  });
}

const Deck = (() => {
  let sections = [];
  let index = 0;
  const listeners = [];

  function scale() {
    const root = document.querySelector('.reveal');
    if (!root) return;
    const s = Math.min(window.innerWidth / 1600, window.innerHeight / 900);
    root.style.setProperty('--deck-scale', s);
  }

  function show(i, push = true) {
    index = Math.max(0, Math.min(i, sections.length - 1));
    sections.forEach((el, k) => el.classList.toggle('present', k === index));
    const counter = document.querySelector('.deck-counter');
    if (counter) counter.textContent = `${index + 1} / ${sections.length}`;
    const prev = document.querySelector('.deck-nav .prev');
    const next = document.querySelector('.deck-nav .next');
    if (prev) prev.disabled = index === 0;
    if (next) next.disabled = index === sections.length - 1;
    if (push) {
      const h = `#/${index}`;
      if (location.hash !== h) history.replaceState(null, '', h);
    }
    listeners.forEach((fn) => fn());
  }

  function fromHash() {
    const m = /^#\/(\d+)/.exec(location.hash);
    return m ? parseInt(m[1], 10) : 0;
  }

  function init() {
    sections = Array.from(document.querySelectorAll('.reveal .slides > section'));
    if (!sections.length) return;

    const reveal = document.querySelector('.reveal');
    const counter = document.createElement('div');
    counter.className = 'deck-counter';
    reveal.appendChild(counter);
    const nav = document.createElement('div');
    nav.className = 'deck-nav';
    nav.innerHTML = '<button class="prev" aria-label="previous slide">&#8249;</button>'
                  + '<button class="next" aria-label="next slide">&#8250;</button>';
    reveal.appendChild(nav);
    nav.querySelector('.prev').addEventListener('click', () => show(index - 1));
    nav.querySelector('.next').addEventListener('click', () => show(index + 1));

    document.addEventListener('keydown', (e) => {
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      switch (e.key) {
        case 'ArrowRight': case 'PageDown': case ' ': show(index + 1); e.preventDefault(); break;
        case 'ArrowLeft':  case 'PageUp':          show(index - 1); e.preventDefault(); break;
        case 'Home': show(0); e.preventDefault(); break;
        case 'End':  show(sections.length - 1); e.preventDefault(); break;
        case 'n': case 'N': reveal.classList.toggle('show-notes'); break;
        default: break;
      }
    });

    window.addEventListener('resize', scale);
    window.addEventListener('hashchange', () => show(fromHash(), false));
    scale();
    show(fromHash(), false);
  }

  return {
    init,
    on: (fn) => listeners.push(fn),
    current: () => sections[index],
    go: show,
  };
})();

/* ===================================================================
   Slide figures for the PADM deck.
   Rendered as SVG at run time so every number on a chart comes from the
   data or from real kernel algebra -- nothing here is a drawn cartoon.
   =================================================================== */

const INK = '#102a43';
const MUTED = '#5d7285';
const BLUE = '#2f78a8';
const ORANGE = '#f28e2b';
const GREEN = '#3a8f57';
const RED = '#c84b4b';
const GRID_C = '#d9e1e7';

function svgElement(name, attributes = {}) {
  const element = document.createElementNS('http://www.w3.org/2000/svg', name);
  Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
  return element;
}

function svgText(x, y, content, attrs = {}) {
  const t = svgElement('text', Object.assign({ x, y, fill: MUTED, 'font-size': 18 }, attrs));
  t.textContent = content;
  return t;
}

/* SVG has no <sub>, so a symbol like t_s has to be built from a tspan that is
   nudged down and set smaller.  Everything on the figures goes through here so
   the labels read as real symbols rather than "t s". */
function svgSym(x, y, base, sub, attrs = {}) {
  const size = Number(attrs['font-size'] || 19);
  const t = svgElement('text', Object.assign({ x, y, fill: MUTED, 'font-size': size }, attrs));
  t.appendChild(document.createTextNode(base));
  if (sub) {
    const s2 = svgElement('tspan', { dy: size * 0.28, 'font-size': size * 0.72 });
    s2.textContent = sub;
    t.appendChild(s2);
  }
  return t;
}

/* ===================================================================
   1. Corrected campaign: the twelve screening designs.
   Straight from results/corrected_boundary_v3/all_samples.csv.
   pa    = metric_pressure_drop_Pa
   mi    = metric_flux_weighted_mixing_index  (higher is better)
   =================================================================== */
const samples = [
  { id: '00000', pa: 13.1849, mi: 0.10482 },
  { id: '00001', pa: 21.2523, mi: 0.11516 },
  { id: '00002', pa: 13.6269, mi: 0.07851 },
  { id: '00003', pa: 16.4689, mi: 0.11258 },
  { id: '00004', pa: 11.1158, mi: 0.09550 },
  { id: '00005', pa: 15.3966, mi: 0.09814 },
  { id: '00006', pa: 12.8789, mi: 0.09212 },
  { id: '00007', pa: 34.0596, mi: 0.16975 },
  { id: '00008', pa: 13.3620, mi: 0.08903 },
  { id: '00009', pa: 14.6275, mi: 0.10363 },
  { id: '00010', pa: 12.3224, mi: 0.10327 },
  { id: '00011', pa: 19.0632, mi: 0.10534 },
];

const BASELINES = [
  { name: 'straight', pa: 2.874, mi: 0.1003 },
  { name: 'symmetric', pa: 10.812, mi: 0.0901 },
  { name: 'strong alt.', pa: 30.509, mi: 0.1450 },
];

/* Minimise pressure, maximise mixing index. */
function nonDominated(data) {
  return data.filter((p) => !data.some((o) =>
    o !== p && o.pa <= p.pa && o.mi >= p.mi && (o.pa < p.pa || o.mi > p.mi)
  ));
}

function renderParetoChart() {
  const svg = document.getElementById('pareto-chart');
  if (!svg || svg.dataset.rendered) return;
  svg.dataset.rendered = 'true';

  const W = 1320, H = 555;
  const m = { left: 118, right: 210, top: 30, bottom: 84 };
  const xMin = 0, xMax = 37, yMin = 0.0, yMax = 0.68;
  const x = (v) => m.left + (v - xMin) / (xMax - xMin) * (W - m.left - m.right);
  const y = (v) => H - m.bottom - (v - yMin) / (yMax - yMin) * (H - m.top - m.bottom);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  [0, 5, 10, 15, 20, 25, 30, 35].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: x(t), x2: x(t), y1: m.top, y2: H - m.bottom, stroke: GRID_C }));
    svg.appendChild(svgText(x(t), H - 46, String(t), { 'text-anchor': 'middle', 'font-size': 19 }));
  });
  [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: y(t), y2: y(t), stroke: GRID_C }));
    svg.appendChild(svgText(m.left - 18, y(t) + 7, t.toFixed(1), { 'text-anchor': 'end', 'font-size': 19 }));
  });

  /* The 0.60 gate and the 20 Pa budget: the two lines the screen was judged against. */
  svg.appendChild(svgElement('line', {
    x1: m.left, x2: W - m.right, y1: y(0.60), y2: y(0.60),
    stroke: GREEN, 'stroke-width': 3, 'stroke-dasharray': '10 7',
  }));
  svg.appendChild(svgText(m.left + 12, y(0.60) - 14, 'predeclared go / no-go gate: mixing index 0.60',
    { fill: GREEN, 'font-size': 21, 'font-weight': 700 }));
  svg.appendChild(svgElement('line', {
    x1: x(20), x2: x(20), y1: m.top, y2: H - m.bottom,
    stroke: '#c84b4b', 'stroke-width': 2, 'stroke-dasharray': '6 6',
  }));
  svg.appendChild(svgText(x(20) + 8, m.top + 84, '20 Pa budget', { fill: '#c84b4b', 'font-size': 19 }));

  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: H - m.bottom, y2: H - m.bottom, stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', { x1: m.left, x2: m.left, y1: m.top, y2: H - m.bottom, stroke: INK, 'stroke-width': 2 }));

  const front = nonDominated(samples).sort((a, b) => a.pa - b.pa);
  svg.appendChild(svgElement('polyline', {
    points: front.map((p) => `${x(p.pa)},${y(p.mi)}`).join(' '),
    fill: 'none', stroke: INK, 'stroke-width': 3,
  }));

  BASELINES.forEach((b) => {
    svg.appendChild(svgElement('rect', {
      x: x(b.pa) - 7, y: y(b.mi) - 7, width: 14, height: 14,
      fill: 'none', stroke: MUTED, 'stroke-width': 3,
    }));
    svg.appendChild(svgText(x(b.pa), y(b.mi) - 16, b.name, { 'text-anchor': 'middle', 'font-size': 18 }));
  });

  samples.forEach((p) => {
    const isFront = front.includes(p);
    const c = svgElement('circle', {
      cx: x(p.pa), cy: y(p.mi), r: isFront ? 9 : 7,
      fill: ORANGE, opacity: isFront ? 1 : 0.55,
      stroke: isFront ? INK : 'none', 'stroke-width': isFront ? 3 : 0,
    });
    const title = svgElement('title');
    title.textContent = `sample ${p.id}: ${p.pa.toFixed(2)} Pa, mixing index ${p.mi.toFixed(4)}`;
    c.appendChild(title);
    svg.appendChild(c);
  });

  const best = samples.reduce((a, b) => (b.mi > a.mi ? b : a));
  svg.appendChild(svgText(x(best.pa) - 14, y(best.mi) - 14, `best: ${best.id} (${best.mi.toFixed(4)})`,
    { 'text-anchor': 'end', fill: INK, 'font-size': 20, 'font-weight': 700 }));

  svg.appendChild(svgText((m.left + W - m.right) / 2, H - 10, 'Pressure drop Δp  [Pa]',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 22 }));
  const yl = svgSym(26, H / 2, 'Mixing index   1 \u2212 \u221AI', 's,flux',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 22 });
  yl.setAttribute('transform', `rotate(-90 26 ${H / 2})`);
  svg.appendChild(yl);
}

/* ===================================================================
   2. Gaussian-process posterior and a UCB acquisition, computed for
   real: squared-exponential kernel, exact posterior, argmax of UCB.
   =================================================================== */
const GP = {
  lengthScale: 0.115,
  signal: 1.0,
  noise: 1e-4,
  kappa: 2.0,
  /* Five evaluations of an expensive black box. */
  X: [0.06, 0.24, 0.42, 0.66, 0.92],
  Y: [-0.55, 0.30, -0.15, 0.72, 0.05],
};

function kernel(a, b) {
  const d = a - b;
  return GP.signal * GP.signal * Math.exp(-(d * d) / (2 * GP.lengthScale * GP.lengthScale));
}

/* Solve A z = b by Gauss-Jordan with partial pivoting. n is 5 here. */
function solve(A, b) {
  const n = b.length;
  const M = A.map((row, i) => row.concat([b[i]]));
  for (let c = 0; c < n; c += 1) {
    let piv = c;
    for (let r = c + 1; r < n; r += 1) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    const d = M[c][c];
    for (let k = c; k <= n; k += 1) M[c][k] /= d;
    for (let r = 0; r < n; r += 1) {
      if (r === c) continue;
      const f = M[r][c];
      for (let k = c; k <= n; k += 1) M[r][k] -= f * M[c][k];
    }
  }
  return M.map((row) => row[n]);
}

function gpPosterior() {
  const n = GP.X.length;
  const K = GP.X.map((xi, i) => GP.X.map((xj, j) => kernel(xi, xj) + (i === j ? GP.noise : 0)));
  const alpha = solve(K, GP.Y);
  const grid = [];
  for (let i = 0; i <= 240; i += 1) {
    const xs = i / 240;
    const ks = GP.X.map((xi) => kernel(xs, xi));
    const mean = ks.reduce((s, k, i2) => s + k * alpha[i2], 0);
    const v = solve(K, ks);
    let var_ = kernel(xs, xs) - ks.reduce((s, k, i2) => s + k * v[i2], 0);
    var_ = Math.max(var_, 1e-9);
    const sd = Math.sqrt(var_);
    grid.push({ x: xs, mean, sd, ucb: mean + GP.kappa * sd });
  }
  return grid;
}

function renderGpChart() {
  const svg = document.getElementById('gp-chart');
  const acq = document.getElementById('ucb-chart');
  if (!svg || svg.dataset.rendered) return;
  svg.dataset.rendered = 'true';

  const grid = gpPosterior();
  const next = grid.reduce((a, b) => (b.ucb > a.ucb ? b : a));

  /* ---- posterior panel ---- */
  const W = 1320, H = 400;
  const m = { left: 70, right: 330, top: 26, bottom: 46 };
  const yMin = -2.6, yMax = 3.2;
  const x = (v) => m.left + v * (W - m.left - m.right);
  const y = (v) => H - m.bottom - (v - yMin) / (yMax - yMin) * (H - m.top - m.bottom);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  const band = grid.map((p) => `${x(p.x)},${y(p.mean + 2 * p.sd)}`)
    .concat(grid.slice().reverse().map((p) => `${x(p.x)},${y(p.mean - 2 * p.sd)}`)).join(' ');
  svg.appendChild(svgElement('polygon', { points: band, fill: BLUE, opacity: 0.18 }));
  svg.appendChild(svgElement('polyline', {
    points: grid.map((p) => `${x(p.x)},${y(p.mean)}`).join(' '),
    fill: 'none', stroke: BLUE, 'stroke-width': 4,
  }));
  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: y(0), y2: y(0), stroke: GRID_C }));

  /* The acquisition, traced on the posterior itself.  With kappa = 2 it is
     exactly the upper edge of the +/-2 sigma band -- which is the clearest way
     to see what "mean plus two error bars of optimism" means. */
  svg.appendChild(svgElement('polyline', {
    points: grid.map((p) => `${x(p.x)},${y(p.ucb)}`).join(' '),
    fill: 'none', stroke: ORANGE, 'stroke-width': 4,
  }));

  GP.X.forEach((xi, i) => {
    svg.appendChild(svgElement('circle', { cx: x(xi), cy: y(GP.Y[i]), r: 9, fill: INK }));
  });

  svg.appendChild(svgElement('line', {
    x1: x(next.x), x2: x(next.x), y1: m.top, y2: H - m.bottom,
    stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '9 6',
  }));

  const legend = [
    [INK, 'the 5 runs we paid for'],
    [BLUE, 'mean μ(x) — best guess'],
    ['rgba(47,120,168,0.35)', '±2σ(x) — our ignorance'],
    [ORANGE, 'acquisition  \u03BC + \u03BA\u03C3'],
  ];
  legend.forEach(([colour, label], i) => {
    const ly = m.top + 16 + i * 34;
    svg.appendChild(svgElement('rect', { x: W - m.right + 14, y: ly - 12, width: 26, height: 14, fill: colour }));
    svg.appendChild(svgText(W - m.right + 50, ly, label, { 'font-size': 18, fill: INK }));
  });
  svg.appendChild(svgText(m.left, 22, 'GP posterior over an expensive objective', { 'font-size': 19 }));

  /* ---- acquisition panel ---- */
  if (!acq || acq.dataset.rendered) return;
  acq.dataset.rendered = 'true';
  const AH = 190;
  const am = { left: 70, right: 330, top: 26, bottom: 34 };
  const uMin = Math.min(...grid.map((p) => p.ucb));
  const uMax = Math.max(...grid.map((p) => p.ucb));
  const ay = (v) => AH - am.bottom - (v - uMin) / (uMax - uMin) * (AH - am.top - am.bottom);
  acq.setAttribute('viewBox', `0 0 ${W} ${AH}`);

  acq.appendChild(svgElement('polygon', {
    points: `${x(0)},${AH - am.bottom} ` + grid.map((p) => `${x(p.x)},${ay(p.ucb)}`).join(' ') + ` ${x(1)},${AH - am.bottom}`,
    fill: ORANGE, opacity: 0.20,
  }));
  acq.appendChild(svgElement('polyline', {
    points: grid.map((p) => `${x(p.x)},${ay(p.ucb)}`).join(' '),
    fill: 'none', stroke: ORANGE, 'stroke-width': 4,
  }));
  acq.appendChild(svgElement('line', {
    x1: x(next.x), x2: x(next.x), y1: am.top, y2: AH - am.bottom,
    stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '9 6',
  }));
  acq.appendChild(svgElement('circle', { cx: x(next.x), cy: ay(next.ucb), r: 10, fill: ORANGE, stroke: INK, 'stroke-width': 3 }));
  acq.appendChild(svgText(m.left, 20, 'acquisition  α(x) = μ(x) + κ·σ(x),  κ = 2', { 'font-size': 19 }));
  acq.appendChild(svgText(W - am.right + 14, AH / 2 + 6,
    `argmax at x = ${next.x.toFixed(3)}`, { 'font-size': 20, fill: INK, 'font-weight': 700 }));
  acq.appendChild(svgText(W - am.right + 14, AH / 2 + 32,
    'milliseconds, not core-hours', { 'font-size': 18 }));
}


/* ===================================================================
   3. One unit cell, drawn FROM the parameters.
   Layout follows FlowCase/alternating_deflector_cad.py exactly:
     [x0, x0+L_s]              centred baffle, thickness t_s
     [x0+L_s, x0+L_s+L_c]      cosine deflectors from both walls
     [x0+L_s+L_c, x0+L_cell]   centred baffle, thickness t_m
   Amplitudes a_weak (bottom) and a_strong (top) are wall intrusions in H.
   =================================================================== */
function cosinePath(x0, x1, amp, fromTop, X, Y, H) {
  /* Matches cosine_bump_points(): a raised-cosine bump of height amp. */
  const pts = [];
  const n = 44;
  for (let i = 0; i <= n; i += 1) {
    const t = i / n;
    const x = x0 + t * (x1 - x0);
    const h = amp * 0.5 * (1 - Math.cos(2 * Math.PI * t));
    pts.push(`${X(x)},${Y(fromTop ? H - h : h)}`);
  }
  return pts;
}

function renderCell(svg, p, opts = {}) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const H = 1.0;
  const Ls = p.L_s, Lc = p.L_c, Lm = p.L_m;
  const L = Ls + Lc + Lm;
  const W = 900, HH = 260;
  const m = { left: 58, right: 58, top: opts.annotate ? 56 : 26, bottom: opts.annotate ? 60 : 26 };
  const X = (v) => m.left + (v / L) * (W - m.left - m.right);
  const Y = (v) => HH - m.bottom - (v / H) * (HH - m.top - m.bottom);
  svg.setAttribute('viewBox', `0 0 ${W} ${HH}`);

  /* fluid domain */
  svg.appendChild(svgElement('rect', {
    x: X(0), y: Y(H), width: X(L) - X(0), height: Y(0) - Y(H),
    fill: '#eaf4fb', stroke: INK, 'stroke-width': 3,
  }));

  /* section shading */
  if (opts.sections) {
    [[0, Ls, 'split'], [Ls, Ls + Lc, 'interaction'], [Ls + Lc, L, 'merge']].forEach(([a, b, name], i) => {
      svg.appendChild(svgElement('rect', {
        x: X(a), y: Y(H), width: X(b) - X(a), height: Y(0) - Y(H),
        fill: i === 1 ? '#f7e6cf' : '#e3eef6', opacity: 0.75,
      }));
      svg.appendChild(svgText((X(a) + X(b)) / 2, Y(H) - 10, name,
        { 'text-anchor': 'middle', 'font-size': 17, fill: MUTED }));
    });
    svg.appendChild(svgElement('rect', {
      x: X(0), y: Y(H), width: X(L) - X(0), height: Y(0) - Y(H),
      fill: 'none', stroke: INK, 'stroke-width': 3,
    }));
  }

  /* centre baffles */
  const baffle = (a, b, t) => svg.appendChild(svgElement('rect', {
    x: X(a), y: Y((H + t) / 2), width: X(b) - X(a), height: Y((H - t) / 2) - Y((H + t) / 2),
    fill: '#f6b66e', stroke: '#ad5c10', 'stroke-width': 2,
  }));
  baffle(0, Ls, p.t_s);
  baffle(Ls + Lc, L, p.t_m);

  /* cosine deflectors: strong on the top wall, weak on the bottom */
  const bot = cosinePath(Ls, Ls + Lc, p.a_weak, false, X, Y, H);
  svg.appendChild(svgElement('polygon', {
    points: `${X(Ls)},${Y(0)} ` + bot.join(' ') + ` ${X(Ls + Lc)},${Y(0)}`,
    fill: '#8bc1df', stroke: BLUE, 'stroke-width': 2,
  }));
  const top = cosinePath(Ls, Ls + Lc, p.a_strong, true, X, Y, H);
  svg.appendChild(svgElement('polygon', {
    points: `${X(Ls)},${Y(H)} ` + top.join(' ') + ` ${X(Ls + Lc)},${Y(H)}`,
    fill: '#8bc1df', stroke: BLUE, 'stroke-width': 2,
  }));

  if (!opts.annotate) return;

  /* ---- dimension annotations ---- */
  const dimX = (a, b, yPix, label, colour = INK) => {
    svg.appendChild(svgElement('line', { x1: X(a), x2: X(b), y1: yPix, y2: yPix, stroke: colour, 'stroke-width': 2 }));
    [a, b].forEach((v) => svg.appendChild(svgElement('line', {
      x1: X(v), x2: X(v), y1: yPix - 6, y2: yPix + 6, stroke: colour, 'stroke-width': 2 })));
    svg.appendChild(svgSym((X(a) + X(b)) / 2, yPix - 9, label[0], label.slice(1),
      { 'text-anchor': 'middle', 'font-size': 19, fill: colour, 'font-weight': 700 }));
  };
  const dimY = (a, b, xPix, label, colour = INK, anchor = 'start') => {
    svg.appendChild(svgElement('line', { x1: xPix, x2: xPix, y1: Y(a), y2: Y(b), stroke: colour, 'stroke-width': 2 }));
    [a, b].forEach((v) => svg.appendChild(svgElement('line', {
      x1: xPix - 6, x2: xPix + 6, y1: Y(v), y2: Y(v), stroke: colour, 'stroke-width': 2 })));
    svg.appendChild(svgSym(xPix + (anchor === 'start' ? 10 : -10), (Y(a) + Y(b)) / 2 + 6,
      label[0], label.slice(1),
      { 'text-anchor': anchor, 'font-size': 19, fill: colour, 'font-weight': 700 }));
  };

  const yBot = HH - 26;
  dimX(0, Ls, yBot, 'Ls');
  dimX(Ls, Ls + Lc, yBot, 'Lc', ORANGE);
  dimX(Ls + Lc, L, yBot, 'Lm');
  dimY(0, H, 22, 'H', MUTED, 'start');

  /* amplitudes, measured from each wall */
  const xPeak = X(Ls + Lc / 2);
  svg.appendChild(svgElement('line', {
    x1: xPeak, x2: xPeak, y1: Y(H), y2: Y(H - p.a_strong), stroke: '#c84b4b', 'stroke-width': 3 }));
  const halo = { stroke: '#fffdf8', 'stroke-width': 5, 'paint-order': 'stroke' };
  svg.appendChild(svgSym(xPeak + 12, Y(H - p.a_strong / 2) + 6, 'a', 'strong',
    Object.assign({ 'font-size': 19, fill: '#c84b4b', 'font-weight': 700 }, halo)));
  svg.appendChild(svgElement('line', {
    x1: xPeak, x2: xPeak, y1: Y(0), y2: Y(p.a_weak), stroke: GREEN, 'stroke-width': 3 }));
  svg.appendChild(svgSym(xPeak + 12, Y(p.a_weak / 2) + 6, 'a', 'weak',
    Object.assign({ 'font-size': 19, fill: GREEN, 'font-weight': 700 }, halo)));

  /* baffle thicknesses */
  svg.appendChild(svgSym(X(Ls / 2), Y(H / 2) - 12, 't', 's',
    { 'text-anchor': 'middle', 'font-size': 19, fill: '#ad5c10', 'font-weight': 700 }));
  svg.appendChild(svgSym(X(Ls + Lc + Lm / 2), Y(H / 2) - 12, 't', 'm',
    { 'text-anchor': 'middle', 'font-size': 19, fill: '#ad5c10', 'font-weight': 700 }));

  /* the peak-to-peak gap that the mesh constraint protects */
  if (opts.gap) {
    const gap = H - p.a_weak - p.a_strong;
    const gx = xPeak - 44;
    svg.appendChild(svgElement('line', {
      x1: gx, x2: gx, y1: Y(p.a_weak), y2: Y(H - p.a_strong),
      stroke: '#7a2fa8', 'stroke-width': 3 }));
    [p.a_weak, H - p.a_strong].forEach((v) => svg.appendChild(svgElement('line', {
      x1: gx - 7, x2: gx + 7, y1: Y(v), y2: Y(v), stroke: '#7a2fa8', 'stroke-width': 3 })));
    /* leader up into the top margin, clear of the baffles and the fill */
    svg.appendChild(svgElement('line', {
      x1: gx, x2: gx, y1: Y(H - p.a_strong), y2: 26,
      stroke: '#7a2fa8', 'stroke-width': 1.5, 'stroke-dasharray': '4 4' }));
    svg.appendChild(svgText(gx, 20, `peak gap ${gap.toFixed(2)} H`,
      { 'text-anchor': 'middle', 'font-size': 19, fill: '#7a2fa8', 'font-weight': 700 }));
  }
}

function renderAllCells() {
  document.querySelectorAll('svg[data-cell]').forEach((svg) => {
    if (svg.dataset.rendered) return;
    svg.dataset.rendered = 'true';
    renderCell(svg, JSON.parse(svg.dataset.cell), {
      annotate: svg.dataset.annotate === 'true',
      sections: svg.dataset.sections === 'true',
      gap: svg.dataset.gap === 'true',
    });
  });
}

/* ===================================================================
   4. Forrester function and the kappa dial.
   f(x) = (6x-2)^2 sin(12x-4), negated so the slide can stay a
   maximisation and match the UCB story told earlier.
   =================================================================== */
function forrester(x) {
  const a = 6 * x - 2;
  return -(a * a * Math.sin(12 * x - 4)) / 6;   /* /6 keeps it near unit scale */
}

/* Sample placement chosen so the kappa dial tells a true story:
   kappa=0 re-measures next to the best known point, kappa=3 lands on the
   real optimum at x=0.758.  (Pushing further drifts past it: kappa=15
   picks 0.781.) */
const FX = [0.05, 0.25, 0.45, 0.62, 0.95];

function forresterPosterior(lengthScale = 0.10, signal = 1.4, noise = 1e-5) {
  const k = (a, b) => {
    const d = a - b;
    return signal * signal * Math.exp(-(d * d) / (2 * lengthScale * lengthScale));
  };
  const Y = FX.map(forrester);
  const K = FX.map((xi, i) => FX.map((xj, j) => k(xi, xj) + (i === j ? noise : 0)));
  const alpha = solve(K, Y);
  const grid = [];
  for (let i = 0; i <= 260; i += 1) {
    const xs = i / 260;
    const ks = FX.map((xi) => k(xs, xi));
    const mean = ks.reduce((s2, kk, i2) => s2 + kk * alpha[i2], 0);
    const v = solve(K, ks);
    const sd = Math.sqrt(Math.max(k(xs, xs) - ks.reduce((s2, kk, i2) => s2 + kk * v[i2], 0), 1e-9));
    grid.push({ x: xs, mean, sd, truth: forrester(xs) });
  }
  return grid;
}

function renderForresterPanel(svg, kappa, opts = {}) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const grid = forresterPosterior();
  const pick = grid.reduce((a, b) => ((b.mean + kappa * b.sd) > (a.mean + kappa * a.sd) ? b : a));
  const trueBest = grid.reduce((a, b) => (b.truth > a.truth ? b : a));

  const W = 520, H = 330;
  const m = { left: 26, right: 18, top: 22, bottom: 44 };
  /* the acquisition must be inside the axes too: at kappa=3 it rises well above
     mean+2sigma and was being clipped off the top of the panel */
  const all = grid.flatMap((p) => [
    p.mean + 2.2 * p.sd, p.mean - 2.2 * p.sd, p.truth, p.mean + kappa * p.sd,
  ]);
  const lo = Math.min(...all), hi = Math.max(...all);
  const x = (v) => m.left + v * (W - m.left - m.right);
  const y = (v) => H - m.bottom - (v - lo) / (hi - lo) * (H - m.top - m.bottom);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  if (opts.truth !== false) {
    svg.appendChild(svgElement('polyline', {
      points: grid.map((p) => `${x(p.x)},${y(p.truth)}`).join(' '),
      fill: 'none', stroke: MUTED, 'stroke-width': 2, 'stroke-dasharray': '7 5',
    }));
  }
  const band = grid.map((p) => `${x(p.x)},${y(p.mean + 2 * p.sd)}`)
    .concat(grid.slice().reverse().map((p) => `${x(p.x)},${y(p.mean - 2 * p.sd)}`)).join(' ');
  svg.appendChild(svgElement('polygon', { points: band, fill: BLUE, opacity: 0.16 }));
  svg.appendChild(svgElement('polyline', {
    points: grid.map((p) => `${x(p.x)},${y(p.mean)}`).join(' '),
    fill: 'none', stroke: BLUE, 'stroke-width': 3,
  }));
  /* the acquisition itself, on the same axes */
  svg.appendChild(svgElement('polyline', {
    points: grid.map((p) => `${x(p.x)},${y(p.mean + kappa * p.sd)}`).join(' '),
    fill: 'none', stroke: ORANGE, 'stroke-width': 3.5,
  }));
  FX.forEach((xi) => svg.appendChild(svgElement('circle', { cx: x(xi), cy: y(forrester(xi)), r: 6, fill: INK })));

  svg.appendChild(svgElement('line', {
    x1: x(pick.x), x2: x(pick.x), y1: m.top, y2: H - m.bottom,
    stroke: ORANGE, 'stroke-width': 3, 'stroke-dasharray': '8 5' }));
  svg.appendChild(svgElement('circle', {
    cx: x(pick.x), cy: y(pick.mean + kappa * pick.sd), r: 8,
    fill: ORANGE, stroke: INK, 'stroke-width': 2.5 }));
  svg.appendChild(svgElement('line', {
    x1: x(trueBest.x), x2: x(trueBest.x), y1: H - m.bottom - 12, y2: H - m.bottom,
    stroke: GREEN, 'stroke-width': 4 }));

  svg.appendChild(svgText(W / 2, H - 20, `picks x = ${pick.x.toFixed(2)}`,
    { 'text-anchor': 'middle', 'font-size': 20, fill: INK, 'font-weight': 700 }));
  svg.appendChild(svgText(W / 2, H - 3, `true optimum at x = ${trueBest.x.toFixed(2)}`,
    { 'text-anchor': 'middle', 'font-size': 16, fill: GREEN }));
}

function renderForresterAll() {
  document.querySelectorAll('svg[data-kappa]').forEach((svg) => {
    if (svg.dataset.rendered) return;
    svg.dataset.rendered = 'true';
    renderForresterPanel(svg, parseFloat(svg.dataset.kappa), { truth: svg.dataset.truth !== 'false' });
  });
}


/* ===================================================================
   5. Pareto front and hypervolume, in objective space.
   BOTH objectives are minimised here (pressure ratio, segregation), so
   "better" is down-and-left and the reference point sits up-and-right.
   The shaded staircase IS the hypervolume: the area the front dominates
   with respect to that reference.
   =================================================================== */
const OBJ = [                      /* schematic designs, minimise both */
  { x: 1.4, y: 0.86 }, { x: 2.1, y: 0.62 }, { x: 3.0, y: 0.47 },
  { x: 4.3, y: 0.38 }, { x: 6.2, y: 0.33 },
  { x: 2.6, y: 0.80 }, { x: 3.8, y: 0.70 }, { x: 5.1, y: 0.55 },
  { x: 1.9, y: 0.92 }, { x: 4.9, y: 0.78 },
];
const REF = { x: 7.4, y: 1.0 };
const NEWPT = { x: 2.5, y: 0.41 };

function paretoOf(pts) {                       /* minimise both */
  return pts.filter((p) => !pts.some((o) =>
    o !== p && o.x <= p.x && o.y <= p.y && (o.x < p.x || o.y < p.y)
  )).sort((a, b) => a.x - b.x);
}

/* Staircase polygon of the region dominated by `front`, up to the reference. */
function dominatedPolygon(front, ref, X, Y) {
  const pts = [];
  pts.push(`${X(front[0].x)},${Y(ref.y)}`);
  front.forEach((p, i) => {
    pts.push(`${X(p.x)},${Y(p.y)}`);
    const nextX = (i + 1 < front.length) ? front[i + 1].x : ref.x;
    pts.push(`${X(nextX)},${Y(p.y)}`);
  });
  pts.push(`${X(ref.x)},${Y(ref.y)}`);
  return pts.join(' ');
}

function renderObjectiveSpace(svg, mode) {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = 1340, H = 560;
  const m = { left: 96, right: 300, top: 34, bottom: 78 };
  const xMax = 8.4, yMax = 1.12;
  const X = (v) => m.left + (v / xMax) * (W - m.left - m.right);
  const Y = (v) => H - m.bottom - (v / yMax) * (H - m.top - m.bottom);
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);

  const front = paretoOf(OBJ);

  if (mode === 'hypervolume') {
    svg.appendChild(svgElement('polygon', {
      points: dominatedPolygon(front, REF, X, Y), fill: BLUE, opacity: 0.20,
    }));
    /* what one more design would add */
    const grown = paretoOf(OBJ.concat([NEWPT]));
    svg.appendChild(svgElement('polygon', {
      points: dominatedPolygon(grown, REF, X, Y), fill: ORANGE, opacity: 0.30,
    }));
    svg.appendChild(svgElement('polygon', {
      points: dominatedPolygon(front, REF, X, Y), fill: '#fffdf8', opacity: 1,
    }));
    svg.appendChild(svgElement('polygon', {
      points: dominatedPolygon(front, REF, X, Y), fill: BLUE, opacity: 0.20,
    }));
  }

  /* axes */
  svg.appendChild(svgElement('line', {
    x1: m.left, x2: W - m.right, y1: H - m.bottom, y2: H - m.bottom, stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', {
    x1: m.left, x2: m.left, y1: m.top, y2: H - m.bottom, stroke: INK, 'stroke-width': 2 }));
  svg.appendChild(svgText((m.left + W - m.right) / 2, H - 34,
    'pressure cost  →  worse', { 'text-anchor': 'middle', fill: INK, 'font-size': 21 }));
  const yl = svgText(30, H / 2, 'segregation  →  worse',
    { 'text-anchor': 'middle', fill: INK, 'font-size': 21 });
  yl.setAttribute('transform', `rotate(-90 30 ${H / 2})`);
  svg.appendChild(yl);
  svg.appendChild(svgText(m.left + 12, m.top + 22, 'better designs are down and to the left',
    { fill: MUTED, 'font-size': 19 }));

  /* dominated points, then the front */
  OBJ.forEach((pt) => {
    const on = front.includes(pt);
    svg.appendChild(svgElement('circle', {
      cx: X(pt.x), cy: Y(pt.y), r: on ? 10 : 7,
      fill: on ? ORANGE : MUTED, opacity: on ? 1 : 0.45,
      stroke: on ? INK : 'none', 'stroke-width': on ? 3 : 0,
    }));
  });
  svg.appendChild(svgElement('polyline', {
    points: front.map((p) => `${X(p.x)},${Y(p.y)}`).join(' '),
    fill: 'none', stroke: INK, 'stroke-width': 3, 'stroke-dasharray': '8 5',
  }));

  if (mode === 'hypervolume') {
    svg.appendChild(svgElement('circle', {
      cx: X(REF.x), cy: Y(REF.y), r: 9, fill: '#c84b4b' }));
    svg.appendChild(svgText(X(REF.x), Y(REF.y) - 16, 'reference point',
      { 'text-anchor': 'middle', fill: '#c84b4b', 'font-size': 20, 'font-weight': 700 }));
    svg.appendChild(svgElement('circle', {
      cx: X(NEWPT.x), cy: Y(NEWPT.y), r: 10,
      fill: ORANGE, stroke: INK, 'stroke-width': 3, 'stroke-dasharray': '4 3' }));
    svg.appendChild(svgText(X(NEWPT.x) - 16, Y(NEWPT.y) + 6, 'candidate',
      { 'text-anchor': 'end', fill: INK, 'font-size': 20, 'font-weight': 700,
        stroke: '#fffdf8', 'stroke-width': 5, 'paint-order': 'stroke' }));
  }

  const legend = mode === 'hypervolume'
    ? [[BLUE, 'hypervolume the front already owns'],
       [ORANGE, 'what the candidate would ADD'],
       ['#c84b4b', 'reference point (declared in config)']]
    : [[ORANGE, 'Pareto front \u2014 non-dominated'],
       [MUTED, 'dominated \u2014 some design beats it on both']];
  legend.forEach(([colour, label], i) => {
    const ly = m.top + 20 + i * 36;
    svg.appendChild(svgElement('rect', {
      x: W - m.right + 16, y: ly - 13, width: 26, height: 15, fill: colour,
      opacity: colour === BLUE || colour === ORANGE ? 0.5 : 1 }));
    svg.appendChild(svgText(W - m.right + 52, ly, label, { 'font-size': 18, fill: INK }));
  });
}

function renderObjectiveSpaces() {
  document.querySelectorAll('svg[data-objspace]').forEach((svg) => {
    if (svg.dataset.rendered) return;
    svg.dataset.rendered = 'true';
    renderObjectiveSpace(svg, svg.dataset.objspace);
  });
}


/* Author attribution on every slide.  Injected rather than written out 67
   times, so it cannot drift out of sync -- and injected INTO the section rather
   than positioned against the viewport, so it scales with the slide. */
const AUTHOR = 'Tomislav Maric \u2013 MMA TU Darmstadt \u2013 maric@mma.tu-darmstadt.de';

function addAuthorLine() {
  document.querySelectorAll('.reveal .slides > section').forEach((section) => {
    if (section.querySelector(':scope > .deck-author')) return;
    const el = document.createElement('div');
    el.className = 'deck-author';
    el.textContent = AUTHOR;
    section.appendChild(el);
  });
}

function renderAll() {
  addAuthorLine();
  renderParetoChart();
  renderGpChart();
  renderAllCells();
  renderForresterAll();
  renderObjectiveSpaces();
  renderTheoryFigures();
  renderDidacticFigures();
}

Deck.on(renderAll);
document.addEventListener('DOMContentLoaded', () => {
  Deck.init();
  renderAll();
});
