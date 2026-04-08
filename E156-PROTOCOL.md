# E156 Protocol — LivingNMA

**Project**: LivingNMA
**Created**: 2026-04-09
**Type**: methodological
**Estimand**: Network treatment effects via graph-theoretic weighted least squares
**Dashboard**: https://mahmood726-cyber.github.io/LivingNMA/

## E156 Body (155 words, 7 sentences)

Can frequentist network meta-analysis — including consistency assessment and treatment ranking — be performed entirely in a web browser? We implemented the graph-theoretic NMA framework (Rücker 2012) using JavaScript matrix algebra: the design matrix X encodes treatment contrasts, weighted least squares yields beta = (X'WX)^(-1)X'WY, and the variance-covariance matrix provides all pairwise comparisons simultaneously. Consistency is assessed globally via the Q statistic decomposition (total Q minus within-design Q) and locally via node-splitting, which separates direct from indirect evidence for each comparison with direct data. Treatment rankings use P-scores — the frequentist analogue of SUCRA — computed as the mean probability of each treatment being superior across all pairwise comparisons. Applied to 12 anti-hypertensive trials across 6 treatments, the tool produces a network graph, league table, forest plots, and rankograms consistent with published NMA results. The living update feature enables sequential re-analysis as new trials emerge, tracking estimate evolution over time. This implementation assumes a common heterogeneity parameter across comparisons and does not handle multi-arm trial covariance adjustment.
