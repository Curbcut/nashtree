# Temporary NASH dependency fork

This R-installable branch is generated with `Rscript cran-bootstrap.R 0 0 1`
from the `codex/nash-fixes` source branch at commit
`67cb458b` in Curbcut/stochtree.

It combines these local bug fixes:

- Correct MCMC grow and prune proposal probabilities (5c5f2e8a).
- Fix extrema calculation in MCMC split helpers (1c846cfd).
- Fix NaN ordering in GFR feature presorting (7d364ad2).

NASH pins the generated R-package commit in its DESCRIPTION Remotes field.
After upstream includes all three fixes in a verified release, remove that
Remotes entry and require the fixed upstream version in Imports.

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

## Negative-binomial BART outcome (0.4.6.9011)

Branch `cma-nb`, based on `cma-rfx` at bd5c5c48 (the `nashtree` rename on top
of the NASH pin 0ac0e7d0). Adds a count outcome to `bart()` so nash can fit
negative-binomial models without GPL code.

**Model.** `OutcomeModel(outcome = "count", link = "log")`:
`y_i ~ NB(mean mu_i, size r)`, `log mu_i = o_i + c + f(x_i) [+ u_g(i)]`, with
`o_i` the optional `offset_train` / `offset_test` (log exposure), `c` the
fixed centering constant `log(sum(y) / sum(exp(o)))` (stored as
`outcome_mean`, the role `y_bar_train` plays for probit), `f` the mean forest
and `u` an optional random intercept or random slopes.

**Gibbs scheme**, in the R loops of `bart()` beside the probit branch:

1. `psi_i = log mu_i - log r`; draw `omega_i ~ PG(y_i + r, psi_i)`.
2. The pseudo-outcome `(y_i - r) / (2 omega_i)` is `N(psi_i, 1 / omega_i)`;
   the residual handed to the forest is that pseudo-outcome net of the offset,
   `c`, `-log r` and the current forest and random-effects fits. The forest
   dataset and the random-effects dataset get variance weights `1 / omega_i`
   (times any `observation_weights`), the global error variance stays at 1.
3. Forest, leaf scale and random effects: the existing samplers, unchanged.
4. `r | y, mu`: slice sampling on `log r` (Neal 2003) of the exact full
   conditional under a `Gamma(a0, b0)` prior (defaults 1 and 0.01).

The plan called for the Zhou and Carin (2012) CRT-Gamma update of `r`. That
update is the conditional of `r` with `p_i = mu_i / (mu_i + r)` held fixed,
which is right when the forest models `logit p_i`. Here the forest models
`log mu_i` and `mu` stays fixed during the `r` step, and the CRT-Gamma draw
then does not leave `p(r | y, mu)` invariant: on simulated data (n = 2,000,
r = 2, `mu` at the truth) it settled at a mean of 2.44 against the exact
conditional's 2.05 (sd 0.12), with lag-one autocorrelation 0.92. The slice
step matches the exact conditional (mean 9.449 vs 9.453, sd 0.755 vs 0.746 in
a second check at r = 10) with lag-one autocorrelation 0.01.

**Pólya-Gamma sampler** (`src/include/stochtree/polya_gamma_sampler.h`,
binding `sample_polya_gamma_cpp(b, z, rng)`), written from the papers only;
no third-party code (bartisan, pgdraw, BayesLogit, pg or other) was read:

- Devroye, L. (2009). On exact simulation algorithms for some distributions
  related to Jacobi theta functions. *Statistics & Probability Letters* 79,
  2251-2259.
- Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference for
  logistic models using Pólya-Gamma latent variables. *JASA* 108, 1339-1349.
- Windle, J., Polson, N. G. and Scott, J. G. (2014). Sampling Pólya-Gamma
  random variates: alternate and approximate techniques. arXiv:1405.0506.
- Michael, J. R., Schucany, W. R. and Haas, R. W. (1976). Generating random
  variates using transformations with multiple roots. *The American
  Statistician* 30, 88-90.

Exact PG(1, z) by Devroye's alternating-series sampler; exact sums for integer
`b <= 50`; otherwise a 100-term gamma series with the remainder replaced by its
exact mean (variance error below 2e-8 of the total) up to `b < 10,000`, and a
moment-matched normal above that (skewness below 0.02). `y + r` is rarely an
integer, so the series path does most of the work.

**Random-effects weights.** `src/random_effects.cpp` now weights each row by
`1 / v_i` in the working-parameter and group-parameter conditionals when the
`RandomEffectsDataset` has variance weights. Without weights the code path is
the old one, and seeded Gaussian BART (two chains, custom basis;
intercept-only) and BCF random-effects draws are bitwise identical to 0ac0e7d0
and to bd5c5c48, which reproduce each other exactly
(`tests/testthat/fixtures/rfx-unweighted-reference.rds`, captured at 0ac0e7d0;
the Gaussian, probit and cloglog reference `outcome-paths-reference.rds`
likewise).

Found while doing this and left unchanged:

- A variance forest writes its per-row variances into the forest dataset only;
  the random-effects dataset is created without weights, so with a variance
  forest the random effects are still sampled as if homoskedastic. Fixing it
  means creating the rfx dataset with weights and copying the forest weights
  into it every iteration in `bart()` and `bcf()`, more than a one-line change.
- `observation_weights` are never passed to the random-effects dataset either,
  so Gaussian and probit random effects ignore them.
- `observation_weights` act as variance multipliers in the C++ leaf sufficient
  statistics (`sum_w += 1 / v_i`), the opposite of the roxygen text
  (`N(mu, sigma^2 / w_i)`, "larger weights increase an observation's
  influence"); upstream has the same code and text. A check with weight 100 on
  one half of the rows gave a fit that followed the other half. The count
  outcome follows the C++ convention: user weights multiply `1 / omega`.

**Out of scope here**: zero inflation (ZINB, added below), the count outcome
in `bcf()`, nash integration, and the alternatives in the feasibility note
(conjugate gamma-Poisson leaves, Laplace reversible jump).

## Zero-inflated negative-binomial BART outcome (0.4.6.9011)

Branch `cma-zinb`, on `cma-nb` (da675755). `OutcomeModel("count", "log",
family = "zinb")`: with probability `pi_i = Phi(c0 + g(x_i))` the row is a
structural zero, otherwise `y_i ~ NB(mu_i, r)` as above. `g` is a probit
"zero forest" with its own `ForestDataset` (no variance weights), `Outcome`,
`ForestModelConfig`, `ForestModel` (tracker) and `ForestSamples`, stored like
the variance forest (`zero_forests`, JSON `forest_1` since a count model has no
variance forest).

**Gibbs scheme**: (1) `s_i ~ Bernoulli(pi_i / (pi_i + (1 - pi_i) NB(0 | mu_i, r)))`
for `y_i = 0`, on the log-odds scale; (2) the zero forest on `s` by the
Albert-Chib step (truncated normal by log-scale inversion); (3) the Pólya-Gamma
count step on rows with `s_i = 0`; rows with `s_i = 1` get variance weight
`Inf`, so `1 / w = 0` and `y / w = 0` in every Gaussian leaf and random-effects
sufficient statistic (tests show the forest and rfx draws are bitwise
independent of those rows' residuals, and an all-excluded leaf draws from its
prior); (4) `r` by the slice step over rows with `s_i = 0`; (5) random effects
on the count part only. `min_samples_leaf` counts structural zeros too
(documented). Initialization: a constant-rate ZINB MLE started from the
excess-zero fraction over the NB method-of-moments fit gives `c`, `qnorm(pi)`
(clamped) and `r`; the first `s` draw comes from the root forests, i.e. from
that fit. `observation_weights` are refused for this family.

**Recovery** (n = 2,000, 4 scenarios x 8 seeds, mean pi 0.38-0.40): rate
correlation >= 0.952, log-rate RMSE <= 0.219, pi correlation >= 0.934, pi RMSE
<= 0.083, response RMSE / sd <= 0.365, rfx correlation >= 0.954, r in the 99%
interval 31 / 32. Plain NB on the same data: rate correlation 0.82-0.85,
log-rate RMSE 0.70-0.76. **Timing** (n = 5,000, 200 count trees, ms per
iteration): Gaussian 9.8, NB 30.3, ZINB 28.9 with 50 zero trees (27.2 with 20).

## Count outcomes in bcf() (0.4.6.9011)

Branch `cma-bcf-count`, on `cma-zinb`. `bcf()` accepts
`OutcomeModel("count", "log")` and `family = "zinb"`.

**Count part.** `log mu_i = o_i + c + mu_c(x_i, pi_hat_i) + [tau_0 + tau_c(x_i)] b(z_i) [+ rfx]`,
`b(z)` the treatment basis (`z` or the adaptive coding). The Pólya-Gamma step
(`bcf_negbin_latent_update` in `R/bcf.R`) writes the pseudo-outcome
`(y - r) / (2 omega) + log r - o` into `resid_train` (the role the latent `z`
plays for probit) and variance weights `1 / omega` into the forest and rfx
datasets; the forests, rfx, GFR and chains then run unchanged. `tau_0` and the
adaptive-coding parameters are sampled in R from unweighted sums in the
Gaussian path, so the count path recomputes them with the precision weights
`omega` (the Gaussian and probit paths are untouched). Variance forests and
`sample_sigma2_global = TRUE` are refused. Priors: prognostic leaf scale as in
NB BART; treatment forest `sigma2_leaf = log_rate_ratio_prior_sd^2 / num_trees`
(default sd 0.5 on the log-rate-ratio scale); `tau_0 ~ N(0, 1)`.

**Zero part** (`zinb`): probit prognostic forest (50 trees, leaf scale
`2 / m`) and treatment forest on the raw `z` (25 trees, probit-BCF calibration
with `delta_max = 0.5`), both on the structural-zero indicators shared with
bart(); no rfx, no intercept, no adaptive coding, univariate treatment only.
JSON forests: mu, tau, zero mu, zero tau (`forest_0`..`forest_3`).

**Recovery** (n = 3,000, confounded binary treatment, offset, 8 seeds per
family): average log rate ratio error <= 0.060 (naive 0.18-0.38), 95% coverage
15 / 16, CATE RMSE <= 0.104, r in 99% interval 16 / 16; zinb tau_z error
<= 0.063 (8 / 8 covered), effect on pi error <= 0.020. **Timing** (n = 5,000,
default 250 + 100 trees, ms per iteration): Gaussian 30.4, NB 63.9, ZINB 68.1.

Found and left unchanged (pre-existing, upstream code):

- Probit BCF with `sample_intercept = TRUE` (the default) omits `tau_0 b(z)`
  from both the truncated-normal mean and the residual handed to the forests
  in the latent step, so the forests partly absorb `tau_0`. The count path
  includes it.
- In the GFR loop, sampling the treatment-forest leaf scale updates the
  prognostic forest's config with the prognostic scale
  (`forest_model_config_mu$update_leaf_model_scale(current_leaf_scale_mu)`);
  the MCMC loop is correct.
- The Gaussian `tau_0` and adaptive-coding updates ignore `observation_weights`.
