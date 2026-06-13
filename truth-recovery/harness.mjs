// ============================================================
// harness.mjs -- Truth-recovery yardstick for LivingNMA.
//
// Wires the app's OWN runNMACore (copied verbatim into engine.mjs) to the
// known-truth network DGP and measures how often its pairwise relative-effect
// CIs cover the TRUE relative effects d_i - d_j.
//
// Candidate fix evaluated: the shipped core is a COMMON-EFFECT NMA (W=1/v); it
// computes a DL tau^2 but never uses it to widen the relative-effect CI. We add
// a RANDOM-EFFECTS refit (same graph-theoretic WLS, W=1/(v+tau^2), reusing the
// engine's own tau^2) and compare coverage.
//
// Truth-first: every number printed comes from seeded simulation here.
// Run:  node truth-recovery/harness.mjs --reps 400
// ============================================================

import { buildNetwork, runNMACore, matCreate, matMul, matT, matInv, normalQuantile }
  from './engine.mjs';
import { generate, makeRng, D_TRUE, SCENARIOS } from './dgp-nma.mjs';

const BASE_SEED = 20260613;
const TRTS = Object.keys(D_TRUE);
// all unordered pairs (i<j) of treatments -> true relative effects
const PAIRS = [];
for (let i = 0; i < TRTS.length; i++)
  for (let j = i + 1; j < TRTS.length; j++) PAIRS.push([TRTS[i], TRTS[j]]);

// Generic graph-theoretic WLS NMA with a supplied weight function -> allPairs CIs.
// getW(v) = 1/v reproduces the shipped fixed-effect model; 1/(v+tau2) is RE.
function nmaWLS(data, net, refIdx, getW) {
  const { trtIdx, T } = net;
  const K = data.length;
  const nonRefIdx = [];
  for (let i = 0; i < T; i++) if (i !== refIdx) nonRefIdx.push(i);
  const X = matCreate(K, T - 1), Y = [], W = [];
  for (let k = 0; k < K; k++) {
    const d = data[k];
    const i = trtIdx[d.t1], j = trtIdx[d.t2];
    Y.push(d.y); W.push(getW(d.v));
    for (let c = 0; c < T - 1; c++) {
      const tIdx = nonRefIdx[c];
      X[k][c] = tIdx === i ? 1 : (tIdx === j ? -1 : 0);
    }
  }
  const Wmat = matCreate(K, K, 0);
  for (let k = 0; k < K; k++) Wmat[k][k] = W[k];
  const Xt = matT(X), XtW = matMul(Xt, Wmat), XtWX = matMul(XtW, X);
  const vcov = matInv(XtWX);
  const beta = matMul(vcov, matMul(XtW, Y.map(y => [y]))).map(r => r[0]);
  const zc = normalQuantile(0.975);
  const eff = (a, b) => {  // effect of treatment a vs b
    const ai = trtIdx[a], bi = trtIdx[b];
    const ic = nonRefIdx.indexOf(ai), jc = nonRefIdx.indexOf(bi);
    let est, se;
    if (ai === refIdx) { est = -beta[jc]; se = Math.sqrt(vcov[jc][jc]); }
    else if (bi === refIdx) { est = beta[ic]; se = Math.sqrt(vcov[ic][ic]); }
    else { est = beta[ic] - beta[jc]; se = Math.sqrt(vcov[ic][ic] + vcov[jc][jc] - 2 * vcov[ic][jc]); }
    return { est, ci: [est - zc * se, est + zc * se] };
  };
  return { eff };
}

// shipped: the verbatim engine (fixed-effect CIs)
function methodShipped(data, net, refIdx) {
  const r = runNMACore(data, net, refIdx);
  return {
    tau2: r.tau2,
    eff: (a, b) => {
      const p = r.allPairs[net.trtIdx[a]][net.trtIdx[b]];
      return { est: p.est, ci: p.ci };
    },
  };
}
// proposed: random-effects refit (reuse engine's DL tau2)
function methodRE(data, net, refIdx) {
  const tau2 = runNMACore(data, net, refIdx).tau2;
  return nmaWLS(data, net, refIdx, v => 1 / (v + tau2));
}

const METHODS = { 'shipped (FE)': methodShipped, 'RE refit': methodRE };
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

export function runCell(studiesPerEdge, scenario, reps, rng) {
  const acc = {};
  for (const n of Object.keys(METHODS)) acc[n] = { cov: 0, m: 0, wSum: 0 };
  for (let r = 0; r < reps; r++) {
    const { data } = generate(studiesPerEdge, scenario, rng);
    const net = buildNetwork(data);
    const refIdx = net.trtIdx['A'];
    for (const [name, build] of Object.entries(METHODS)) {
      let mdl; try { mdl = build(data, net, refIdx); } catch { continue; }
      for (const [a, b] of PAIRS) {
        const trueRel = D_TRUE[a] - D_TRUE[b];
        let o; try { o = mdl.eff(a, b); } catch { continue; }
        if (!o || !isFinite(o.ci[0]) || !isFinite(o.ci[1])) continue;
        acc[name].m++;
        acc[name].wSum += o.ci[1] - o.ci[0];
        if (o.ci[0] <= trueRel && trueRel <= o.ci[1]) acc[name].cov++;
      }
    }
  }
  const res = {};
  for (const [name, a] of Object.entries(acc))
    res[name] = { coverage: a.m ? +(a.cov / a.m).toFixed(4) : null,
                  meanWidth: a.m ? +(a.wSum / a.m).toFixed(4) : null, m: a.m };
  return res;
}

export function runGrid({ reps = 400, perEdge = [1, 2, 4], scenarios = SCENARIOS } = {}) {
  const rng = makeRng(BASE_SEED);
  const grid = [];
  for (const scen of scenarios)
    for (const spe of perEdge) grid.push({ scen, spe, results: runCell(spe, scen, reps, rng) });
  return grid;
}

export function summarize(grid, filter = () => true) {
  const names = Object.keys(grid[0].results);
  const out = {};
  for (const name of names) {
    const cov = [];
    for (const c of grid) if (filter(c) && c.results[name].coverage != null) cov.push(c.results[name].coverage);
    out[name] = { meanCoverage: +mean(cov).toFixed(4) };
  }
  return out;
}

const isMain = process.argv[1]?.endsWith('harness.mjs');
if (isMain) {
  const i = process.argv.indexOf('--reps');
  const reps = i >= 0 ? Number(process.argv[i + 1]) : 400;
  const t0 = Date.now();
  const grid = runGrid({ reps });
  const s = summarize(grid);
  console.log(`\n# Truth-recovery yardstick -- LivingNMA`);
  console.log(`reps=${reps}/cell  true d: A0 B0.4 C0.7 D1.0  seed=${BASE_SEED}\n`);
  console.log('## Mean coverage of TRUE relative effects (all 6 pairs) over all cells\n');
  console.log('method          meanCov');
  for (const [n, v] of Object.entries(s)) console.log(n.padEnd(14), String(v.meanCoverage).padStart(7));
  console.log('\n## Per-cell coverage: shipped (FE) vs RE refit\n');
  console.log('scenario   studies/edge   FE.cov   RE.cov   FE.width  RE.width');
  for (const c of grid) {
    console.log(c.scen.padEnd(10), String(c.spe).padStart(10),
      String(c.results['shipped (FE)'].coverage).padStart(9),
      String(c.results['RE refit'].coverage).padStart(8),
      String(c.results['shipped (FE)'].meanWidth).padStart(10),
      String(c.results['RE refit'].meanWidth).padStart(9));
  }
  console.log(`\n(${(Date.now() - t0) / 1000}s)`);
}
