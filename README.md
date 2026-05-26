# LivingNMA

[![ci](https://github.com/mahmood726-cyber/LivingNMA/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/LivingNMA/actions/workflows/ci.yml) [![codeql](https://github.com/mahmood726-cyber/LivingNMA/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/LivingNMA/actions/workflows/codeql.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)

Living Network Meta-Analysis engine — frequentist NMA in the browser.

## Methods
- **Graph-theoretic NMA** (Rücker 2012) — weighted least squares via design matrix
- **Consistency** — global Q decomposition + node-splitting (local)
- **P-scores** — frequentist SUCRA analogue
- **DL heterogeneity** — common tau² across comparisons

## Outputs
- Network graph (circular layout, edge thickness = study count)
- League table (all pairwise, significance-colored)
- Forest plot (vs reference treatment)
- Rankogram (P-score bar chart)
- Consistency assessment (global + node-splitting table)

## Examples
- Anti-hypertensives (6 treatments, 12 contrasts)
- Antidepressants (4 treatments, 8 contrasts)
- DOACs for AF (5 treatments vs Warfarin)

## Live Demo
https://mahmood726-cyber.github.io/LivingNMA/

## Author
Mahmood Ahmad, Tahir Heart Institute
