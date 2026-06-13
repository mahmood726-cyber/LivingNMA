# Truth-recovery yardstick — LivingNMA

**Verdict: MAJOR measured improvement — the shipped NMA reports fixed-effect CIs
that severely under-cover the truth under heterogeneity; a random-effects refit
fixes it.**

## The defect (measured)
`runNMACore` is a graph-theoretic frequentist NMA whose relative-effect CIs come
from a **common-effect (fixed-effect) weighted least squares** (`W = 1/v`). It
*does* compute a DerSimonian–Laird `tau^2`, but uses it **only for the prediction
interval** — never to widen the relative-effect confidence intervals. So whenever
real between-study heterogeneity exists (the usual case), the reported CIs are too
narrow and miss the true relative effect far more than 5% of the time.

## Method
- DGP (`dgp-nma.mjs`): a connected 4-treatment network (A,B,C,D; edges A-B, A-C,
  B-C, B-D, C-D) with **known** basic parameters `d = {A:0, B:0.4, C:0.7, D:1.0}`
  (so every true relative effect `d_i-d_j` is known) and a known common `tau^2`.
  Consistency holds by construction. Seeded → reproducible.
- Engine (`engine.mjs`): `runNMACore` + matrix utils + `buildNetwork` copied
  **verbatim** from `index.html`. The SAME math the app ships is measured.
- The candidate fix is a **random-effects refit**: the same WLS, but
  `W = 1/(v + tau^2)`, reusing the engine's own DL `tau^2`. Point estimates barely
  move; only the variance widens correctly.
- 400 reps/cell, true `d` above, `studies/edge ∈ {1,2,4}`, `tau^2 ∈ {0, 0.02,
  0.06, 0.15}`. Coverage averaged over all 6 treatment pairs.

## Results (coverage of the TRUE relative effects)

| scenario  | studies/edge | FE coverage | RE coverage |
|-----------|-------------:|------------:|------------:|
| tau0      | 1            |   0.945     |   0.951     |
| tau0      | 4            |   0.960     |   0.967     |
| tau_low   | 4            |   0.862     |   0.925     |
| tau_mod   | 4            |   0.746     |   0.899     |
| tau_high  | 1            |   0.660     |   0.789     |
| tau_high  | 2            |   0.620     |   0.867     |
| tau_high  | 4            |   **0.586** |   **0.910** |

**Mean over all 12 cells: FE 0.806 → RE 0.905 (+10pp).**

## Findings (all measured)
1. **The shipped fixed-effect NMA under-covers true relative effects badly under
   heterogeneity** — down to 0.59 at `tau^2=0.15`, vs the nominal 0.95. The more
   studies per edge, the *worse* FE gets (more precision concentrated on a CI that
   ignores tau^2).
2. **A random-effects refit recovers the truth far better** — mean +10pp overall,
   **+32pp** in the worst cell (`tau_high`, 4 studies/edge: 0.586 → 0.910). This is
   exactly the netmeta default (random-effects consistency model).
3. **No cost when there is no heterogeneity** — at `tau^2=0` both sit at ~0.95;
   RE never recovers worse than FE in any cell.
4. Honest caveat: RE still slightly under-covers under `tau_high` with only 1
   study/edge (0.789), because DL `tau^2` is poorly estimated with so few studies.
   That is a small-sample tau-estimation limit, not a flaw in the RE principle —
   and still a 13pp improvement over FE's 0.660.

## Recommendation (specific integration point)
Make the random-effects model the default (or a user toggle) for the
relative-effect CIs. The reference implementation is `harness.mjs::nmaWLS` with
`getW = v => 1/(v + tau2)`. Concretely: in `runNMACore`, after `tau2` is computed,
run a second WLS pass with `W_k = 1/(v_k + tau2)` and use that `beta`/`vcov` for
`comparisons` and `allPairs`. (Keep the FE pass for Q, node-splitting and the hat
matrix, which legitimately use within-study weights.)

## What did NOT transfer
NPE/conformal/SBC/PartialID are estimator-of-μ machinery; LivingNMA is a
linear-algebra NMA engine, so only the known-truth harness + a network DGP
transferred. No runtime dependency added.

## Reproduce
```
node truth-recovery/harness.mjs --reps 400
node --test truth-recovery/test-truth-recovery.mjs
```
