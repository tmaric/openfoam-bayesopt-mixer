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
const GRID = '#d9e1e7';

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
    svg.appendChild(svgElement('line', { x1: x(t), x2: x(t), y1: m.top, y2: H - m.bottom, stroke: GRID }));
    svg.appendChild(svgText(x(t), H - 46, String(t), { 'text-anchor': 'middle', 'font-size': 19 }));
  });
  [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6].forEach((t) => {
    svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: y(t), y2: y(t), stroke: GRID }));
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
  const yl = svgText(26, H / 2, 'Mixing index  1 − √I s,flux',
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
  const yMin = -2.6, yMax = 2.9;
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
  svg.appendChild(svgElement('line', { x1: m.left, x2: W - m.right, y1: y(0), y2: y(0), stroke: GRID }));

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
    [ORANGE, 'UCB says: spend run 6 here'],
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
    svg.appendChild(svgText((X(a) + X(b)) / 2, yPix - 9, label,
      { 'text-anchor': 'middle', 'font-size': 19, fill: colour, 'font-weight': 700 }));
  };
  const dimY = (a, b, xPix, label, colour = INK, anchor = 'start') => {
    svg.appendChild(svgElement('line', { x1: xPix, x2: xPix, y1: Y(a), y2: Y(b), stroke: colour, 'stroke-width': 2 }));
    [a, b].forEach((v) => svg.appendChild(svgElement('line', {
      x1: xPix - 6, x2: xPix + 6, y1: Y(v), y2: Y(v), stroke: colour, 'stroke-width': 2 })));
    svg.appendChild(svgText(xPix + (anchor === 'start' ? 10 : -10), (Y(a) + Y(b)) / 2 + 6, label,
      { 'text-anchor': anchor, 'font-size': 19, fill: colour, 'font-weight': 700 }));
  };

  const yBot = HH - 26;
  dimX(0, Ls, yBot, 'L s');
  dimX(Ls, Ls + Lc, yBot, 'L c', ORANGE);
  dimX(Ls + Lc, L, yBot, 'L m');
  dimY(0, H, 22, 'H', MUTED, 'start');

  /* amplitudes, measured from each wall */
  const xPeak = X(Ls + Lc / 2);
  svg.appendChild(svgElement('line', {
    x1: xPeak, x2: xPeak, y1: Y(H), y2: Y(H - p.a_strong), stroke: '#c84b4b', 'stroke-width': 3 }));
  const halo = { stroke: '#fffdf8', 'stroke-width': 5, 'paint-order': 'stroke' };
  svg.appendChild(svgText(xPeak + 12, Y(H - p.a_strong / 2) + 6, 'a strong',
    Object.assign({ 'font-size': 19, fill: '#c84b4b', 'font-weight': 700 }, halo)));
  svg.appendChild(svgElement('line', {
    x1: xPeak, x2: xPeak, y1: Y(0), y2: Y(p.a_weak), stroke: GREEN, 'stroke-width': 3 }));
  svg.appendChild(svgText(xPeak + 12, Y(p.a_weak / 2) + 6, 'a weak',
    Object.assign({ 'font-size': 19, fill: GREEN, 'font-weight': 700 }, halo)));

  /* baffle thicknesses */
  svg.appendChild(svgText(X(Ls / 2), Y(H / 2) - 12, 't s',
    { 'text-anchor': 'middle', 'font-size': 19, fill: '#ad5c10', 'font-weight': 700 }));
  svg.appendChild(svgText(X(Ls + Lc + Lm / 2), Y(H / 2) - 12, 't m',
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

Reveal.initialize({
  /* Slides are laid out by CSS flex (title pinned top, body centred in the rest),
     which needs full-height sections.  Reveal's default vertical centring sizes
     each section to its content instead, leaving the dead space at the bottom. */
  center: false,
  hash: true,
  height: 900,
  margin: 0.035,
  navigationMode: 'linear',
  plugins: [RevealMath.KaTeX, RevealNotes],
  slideNumber: 'c/t',
  transition: 'fade',
  width: 1600,
});

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
}

Reveal.on('ready', renderAll);
Reveal.on('slidechanged', renderAll);
