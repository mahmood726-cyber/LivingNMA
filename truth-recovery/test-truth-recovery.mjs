// node --test truth-recovery/test-truth-recovery.mjs
// Measured invariants for the LivingNMA truth-recovery yardstick. Seeded; no
// hand-entered numbers.
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { generate, makeRng, D_TRUE, SCENARIOS } from './dgp-nma.mjs';
import { buildNetwork } from './engine.mjs';
import { runCell, runGrid, summarize } from './harness.mjs';

describe('network DGP', () => {
  it('is reproducible for a fixed seed', () => {
    const a = generate(2, 'tau_mod', makeRng(7));
    const b = generate(2, 'tau_mod', makeRng(7));
    assert.deepEqual(a.data, b.data);
  });
  it('produces a connected network for every scenario', () => {
    const rng = makeRng(1);
    for (const scen of SCENARIOS) {
      const { data } = generate(2, scen, rng);
      assert.ok(buildNetwork(data).connected, `disconnected for ${scen}`);
    }
  });
});

describe('Truth-recovery (measured)', () => {
  it('the shipped fixed-effect NMA UNDER-covers true relative effects under heterogeneity', () => {
    const rng = makeRng(20260613);
    const hi = runCell(4, 'tau_high', 400, rng);
    // FE CIs ignore tau^2 -> badly too narrow when between-study heterogeneity exists.
    assert.ok(hi['shipped (FE)'].coverage < 0.75,
      `FE coverage under high het = ${hi['shipped (FE)'].coverage} (expected severe under-coverage)`);
  });

  it('a random-effects refit (engine tau^2) recovers far more truth; never worse than FE', () => {
    const grid = runGrid({ reps: 300 });
    const s = summarize(grid);
    assert.ok(s['RE refit'].meanCoverage > s['shipped (FE)'].meanCoverage + 0.05,
      `RE ${s['RE refit'].meanCoverage} not >> FE ${s['shipped (FE)'].meanCoverage}`);
    // RE wins biggest under high heterogeneity.
    const hi = summarize(grid, c => c.scen === 'tau_high');
    assert.ok(hi['RE refit'].meanCoverage - hi['shipped (FE)'].meanCoverage > 0.15,
      `RE-FE gap under high het only ${hi['RE refit'].meanCoverage - hi['shipped (FE)'].meanCoverage}`);
  });

  it('at tau^2=0 the RE refit does not hurt (both near nominal)', () => {
    const grid = runGrid({ reps: 300 });
    const s0 = summarize(grid, c => c.scen === 'tau0');
    assert.ok(s0['RE refit'].meanCoverage > 0.92, `RE tau0 coverage ${s0['RE refit'].meanCoverage}`);
    assert.ok(s0['shipped (FE)'].meanCoverage > 0.92, `FE tau0 coverage ${s0['shipped (FE)'].meanCoverage}`);
  });
});
