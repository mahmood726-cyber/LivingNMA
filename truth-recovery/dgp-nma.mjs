// ============================================================
// dgp-nma.mjs -- Known-truth DGP for NETWORK meta-analysis.
//
// LivingNMA runs a graph-theoretic frequentist NMA. Its CIs come from a
// COMMON-EFFECT (fixed-effect) weighted least squares (W = 1/v). It DOES compute
// a DL tau^2, but uses it only for the prediction interval -- never to widen the
// relative-effect CIs. This DGP tests what that costs: it injects a known set of
// true basic parameters d (so every true relative effect d_i - d_j is known) and
// a known common between-study heterogeneity tau^2, then measures how often the
// engine's pairwise CIs cover the truth.
//
// Network (consistency holds by construction): treatments A,B,C,D with edges
//   A-B, A-C, B-C (a closed loop), B-D, C-D. Studies are distributed across edges.
//   For a study on contrast (t1,t2): theta = (d_t1 - d_t2) + sqrt(tau2)*Z,
//   y = theta + se*Z',  se ~ log-uniform.  Data shape matches the app:
//   {study, t1, t2, y, se, v}.
//
// ESTIMAND: the true relative effects d_i - d_j (consistency model target).
// Seeded -> reproducible. Standalone.
// ============================================================

export const SCENARIOS = ['tau0', 'tau_low', 'tau_mod', 'tau_high'];
const TAU2 = { tau0: 0, tau_low: 0.02, tau_mod: 0.06, tau_high: 0.15 };

// true basic parameters (log scale), reference = A = 0
export const D_TRUE = { A: 0, B: 0.40, C: 0.70, D: 1.00 };
const EDGES = [['A', 'B'], ['A', 'C'], ['B', 'C'], ['B', 'D'], ['C', 'D']];

export function makeRng(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function randn(rng) {
  let u1 = rng(), u2 = rng();
  if (u1 < 1e-12) u1 = 1e-12;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}
function drawSe(rng, seLo, seHi) {
  const lo = Math.log(seLo), hi = Math.log(seHi);
  return Math.exp(lo + (hi - lo) * rng());
}

// studiesPerEdge studies on each edge -> k = studiesPerEdge * |EDGES| contrasts.
export function generate(studiesPerEdge, scenario, rng, { seLo = 0.12, seHi = 0.55 } = {}) {
  const tau2 = TAU2[scenario] ?? 0;
  const sd = Math.sqrt(tau2);
  const data = [];
  let sid = 0;
  for (const [t1, t2] of EDGES) {
    const trueRel = D_TRUE[t1] - D_TRUE[t2];
    for (let s = 0; s < studiesPerEdge; s++) {
      const se = drawSe(rng, seLo, seHi);
      const theta = trueRel + sd * randn(rng);
      const y = theta + se * randn(rng);
      data.push({ study: `s${sid++}`, t1, t2, y, se, v: se * se });
    }
  }
  return { data, D_TRUE, tau2, info: { k: data.length } };
}
