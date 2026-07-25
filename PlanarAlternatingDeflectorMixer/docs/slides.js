const samples = [
  { id: '00000', phase: 'sobol', dp: 0.00139438, segregation: 0.818833 },
  { id: '00001', phase: 'sobol', dp: 0.00262067, segregation: 0.728358 },
  { id: '00002', phase: 'sobol', dp: 0.0015977, segregation: 0.860598 },
  { id: '00003', phase: 'sobol', dp: 0.00192139, segregation: 0.64428 },
  { id: '00004', phase: 'sobol', dp: 0.0018918, segregation: 0.704327 },
  { id: '00005', phase: 'sobol', dp: 0.00212679, segregation: 0.656652 },
  { id: '00006', phase: 'sobol', dp: 0.00254635, segregation: 0.809064 },
  { id: '00007', phase: 'sobol', dp: 0.00199872, segregation: 0.691184 },
  { id: '00008', phase: 'bo', dp: 0.00139317, segregation: 0.760204 },
  { id: '00009', phase: 'bo', dp: 0.00163849, segregation: 0.624195 },
  { id: '00010', phase: 'bo', dp: 0.00186858, segregation: 0.582741 },
  { id: '00011', phase: 'bo', dp: 0.00238017, segregation: 0.535467 },
  { id: '00012', phase: 'bo', dp: 0.00141802, segregation: 0.721376 },
  { id: '00013', phase: 'bo', dp: 0.00105026, segregation: 0.873552 },
  { id: '00014', phase: 'bo', dp: 0.00159905, segregation: 0.643871 },
  { id: '00015', phase: 'bo', dp: 0.001924, segregation: 0.52742 },
  { id: '00016', phase: 'bo', dp: 0.00109318, segregation: 0.871843 },
  { id: '00017', phase: 'bo', dp: 0.00120549, segregation: 0.819882 },
  { id: '00018', phase: 'bo', dp: 0.00193124, segregation: 0.525891 },
  { id: '00019', phase: 'bo', dp: 0.00178059, segregation: 0.583896 },
  { id: '00020', phase: 'bo', dp: 0.00154777, segregation: 0.670629 },
  { id: '00021', phase: 'bo', dp: 0.00129664, segregation: 0.821584 },
  { id: '00022', phase: 'bo', dp: 0.00100487, segregation: 0.880986 },
  { id: '00023', phase: 'bo', dp: 0.00129633, segregation: 0.789391 },
  { id: '00024', phase: 'bo', dp: 0.00244679, segregation: 0.542703 },
  { id: '00025', phase: 'bo', dp: 0.00108038, segregation: 0.883315 },
  { id: '00026', phase: 'bo', dp: 0.00155495, segregation: 0.721264 },
  { id: '00027', phase: 'bo', dp: 0.00198601, segregation: 0.499471 },
];

function nonDominated(data) {
  return data.filter((point) => !data.some((other) =>
    other !== point &&
    other.dp <= point.dp &&
    other.segregation <= point.segregation &&
    (other.dp < point.dp || other.segregation < point.segregation)
  ));
}

function svgElement(name, attributes = {}) {
  const element = document.createElementNS('http://www.w3.org/2000/svg', name);
  Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
  return element;
}

function renderParetoChart() {
  const svg = document.getElementById('pareto-chart');
  if (!svg || svg.dataset.rendered) return;
  svg.dataset.rendered = 'true';

  const width = 1320;
  const height = 555;
  const margin = { left: 115, right: 40, top: 28, bottom: 82 };
  const xMin = Math.log10(0.00095);
  const xMax = Math.log10(0.00275);
  const yMin = 0.08;
  const yMax = 0.53;
  const x = (value) => margin.left + (Math.log10(value) - xMin) / (xMax - xMin) * (width - margin.left - margin.right);
  const y = (quality) => height - margin.bottom - (quality - yMin) / (yMax - yMin) * (height - margin.top - margin.bottom);

  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);

  const xTicks = [0.001, 0.0012, 0.0015, 0.002, 0.0025];
  const yTicks = [0.1, 0.2, 0.3, 0.4, 0.5];

  xTicks.forEach((tick) => {
    const px = x(tick);
    svg.appendChild(svgElement('line', { x1: px, x2: px, y1: margin.top, y2: height - margin.bottom, stroke: '#d9e1e7', 'stroke-width': 1 }));
    const label = svgElement('text', { x: px, y: height - 42, 'text-anchor': 'middle', fill: '#5d7285', 'font-size': 19 });
    label.textContent = `${(tick * 1000).toFixed(tick === 0.001 ? 1 : 2)}`;
    svg.appendChild(label);
  });

  yTicks.forEach((tick) => {
    const py = y(tick);
    svg.appendChild(svgElement('line', { x1: margin.left, x2: width - margin.right, y1: py, y2: py, stroke: '#d9e1e7', 'stroke-width': 1 }));
    const label = svgElement('text', { x: margin.left - 20, y: py + 7, 'text-anchor': 'end', fill: '#5d7285', 'font-size': 19 });
    label.textContent = tick.toFixed(1);
    svg.appendChild(label);
  });

  svg.appendChild(svgElement('line', { x1: margin.left, x2: width - margin.right, y1: height - margin.bottom, y2: height - margin.bottom, stroke: '#102a43', 'stroke-width': 2 }));
  svg.appendChild(svgElement('line', { x1: margin.left, x2: margin.left, y1: margin.top, y2: height - margin.bottom, stroke: '#102a43', 'stroke-width': 2 }));

  const front = nonDominated(samples).sort((a, b) => a.dp - b.dp);
  const points = front.map((point) => `${x(point.dp)},${y(1 - point.segregation)}`).join(' ');
  svg.appendChild(svgElement('polyline', { points, fill: 'none', stroke: '#102a43', 'stroke-width': 3 }));

  samples.forEach((point) => {
    const isFront = front.includes(point);
    const circle = svgElement('circle', {
      cx: x(point.dp),
      cy: y(1 - point.segregation),
      r: isFront ? 9 : 7,
      fill: point.phase === 'sobol' ? '#2f78a8' : '#f28e2b',
      opacity: isFront ? 1 : 0.58,
      stroke: isFront ? '#102a43' : 'none',
      'stroke-width': isFront ? 3 : 0,
    });
    const title = svgElement('title');
    title.textContent = `sample ${point.id}: Jdp=${point.dp}, 1-Is=${(1-point.segregation).toFixed(4)}`;
    circle.appendChild(title);
    svg.appendChild(circle);
  });

  ['00022', '00027'].forEach((id) => {
    const point = samples.find((item) => item.id === id);
    const label = svgElement('text', {
      x: x(point.dp) + (id === '00022' ? 14 : -14),
      y: y(1 - point.segregation) - 14,
      'text-anchor': id === '00022' ? 'start' : 'end',
      fill: '#102a43',
      'font-size': 20,
      'font-weight': 700,
    });
    label.textContent = id;
    svg.appendChild(label);
  });

  const xLabel = svgElement('text', { x: (margin.left + width - margin.right) / 2, y: height - 8, 'text-anchor': 'middle', fill: '#102a43', 'font-size': 22 });
  xLabel.textContent = 'Kinematic pressure drop Jdp (×10⁻³ m²/s², logarithmic axis)';
  svg.appendChild(xLabel);

  const yLabel = svgElement('text', { x: 26, y: height / 2, transform: `rotate(-90 26 ${height / 2})`, 'text-anchor': 'middle', fill: '#102a43', 'font-size': 22 });
  yLabel.textContent = 'Mixing quality 1 − Is';
  svg.appendChild(yLabel);
}

Reveal.initialize({
  hash: true,
  height: 900,
  margin: 0.035,
  navigationMode: 'linear',
  plugins: [RevealMath.KaTeX, RevealNotes],
  slideNumber: 'c/t',
  transition: 'fade',
  width: 1600,
});

Reveal.on('ready', renderParetoChart);
Reveal.on('slidechanged', renderParetoChart);
