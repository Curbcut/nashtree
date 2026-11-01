# NASH changes to the fork

The record of each change nashtree makes to stochtree, oldest first. The fork
started as a temporary pin; since the nashtree rename it is maintained long
term (see README.md).

## Sampler fixes

This R-installable branch is generated with `Rscript cran-bootstrap.R 0 0 1`
from the `codex/nash-fixes` source branch at commit
`67cb458b` in Curbcut/stochtree.

It combines these local bug fixes:

- Correct MCMC grow and prune proposal probabilities (5c5f2e8a).
- Fix extrema calculation in MCMC split helpers (1c846cfd).
- Fix NaN ordering in GFR feature presorting (7d364ad2).

NASH pins the generated R-package commit in its DESCRIPTION Remotes field.

## Numeric preprocessing follow-up

The `codex/numeric-preprocessing` branch starts at the exact generated package
commit `06fbcc6bf3bd57c7d5caee475895fa5db32646d6`. It batches ordinary numeric
covariate binding in training and prediction preprocessing while preserving
categorical processing and retaining pairwise behavior for attributed columns.
Native sampler sources are unchanged. The added numeric-binding tests compare
full preprocessing metadata and seeded BART draws and split counts with the
previous numeric loop.

## Bulk sample retention follow-up

`codex/bulk-forest-retention` starts at `4a2ef47d74b671d6d900efa22412581b80186894`.
`ForestSamples$retain_samples()` validates strictly increasing zero-based IDs
before moving surviving forest pointers into a new container. Empty IDs remove
all samples; numerical contents and order are preserved. Model-level variance
arrays and sample maps remain the caller's responsibility. This is post-fit
compaction, not an importance-only sampler or a change to iteration budgets.

Focused forest tests and the complete standard installed suite passed (two
configured skips and 13 existing model-setting warnings). Exact JSON equality
with descending single-sample deletion and invalid-input nonmutation are tested.
An isolated actual 1,000-sample/64-retained NASH candidate measured 6 ms for bulk
retention versus 8–9 ms for the old loop, with identical compact model JSON.
This bounded result is not a full-fit speedup or native RSS measurement.

## Causal sampler correctness (0.4.6.9006)

Based on the complete NASH pin be00e7689f78fa565bab1b6bdbb1a887b8c332cb.
Preserves all earlier NASH fixes and forest-retention features.

* Correct the Gaussian random-effects inverse-gamma conditional to
  shape a + J/2 and scale b + sum(xi^2)/2 (upstream-source fix 55551c88).
* Reconstruct tree partitions in parent-before-child order, including when
  MCMC prune/grow operations recycle node IDs. Preserve deleted-node bookkeeping.

These are correctness fixes, not evidence of posterior convergence.
Regression evidence is included in inst/nash-correctness.

## nashtree identity and native forests (0.4.6.9010)

`126a6320` renames the package, DLL and configure identity to `nashtree`, with
no sampler, prediction or serialization change; seeded BART, random-intercept
BART and BCF fits serialize identically to the stochtree-named package. The
commits after it add the `nashtree_forest_v1` native core, its inference-only
R bindings, the producer-neutral writer and build receipts. Evidence is in
`native/RESULTS.md` and `native/evidence/`. nash pins `ff3eb820` for both
training and its inference bindings.
