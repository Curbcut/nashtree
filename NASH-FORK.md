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
adaptive-coding parameters are sampled in R; the count path weights their
sums by the precisions `omega` (Gaussian `observation_weights` were added to
the same sums later, see below). Variance forests and
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

Found while doing this (pre-existing, upstream code): the probit `tau_0`
omission listed under "Known upstream bugs" below, and two bugs fixed in the
next section (the GFR treatment leaf scale, unweighted Gaussian `tau_0` and
adaptive coding).

## bcf() GFR treatment leaf scale and weighted tau_0 / adaptive coding

Two upstream bugs in `R/bcf.R`, fixed on a branch from `cma-rfx` (84f5ccbb).
`tests/testthat/fixtures/make-bcf-fixes-reference.R` captured seeded draws on
84f5ccbb before the fixes; `test-bcf-fixes.R` checks them.

- **GFR treatment leaf scale.** With
  `treatment_effect_forest_params$sample_sigma2_leaf = TRUE`, each GFR
  iteration drew a treatment leaf scale and then passed the prognostic scale
  to the prognostic config, so the treatment forest kept its initial scale for
  all GFR sweeps (the MCMC loop was right). The GFR loop now updates the
  treatment config. The test instruments `createForestModelConfig()` and
  checks the treatment config receives, and holds, every stored
  `sigma2_leaf_tau` draw in order.
- **Weighted Gaussian `tau_0` and adaptive coding.** The C++ leaf statistics
  read `observation_weights` as variance multipliers `v_i`
  (`sum_w += 1 / v_i`), but the R-side `tau_0` and `b_0` / `b_1` updates used
  unweighted sums. They now use precisions `w_i = 1 / v_i` (`Z'W r`, `Z'W Z`,
  weighted `s_tt`, `s_ty`; the adaptive-coding prior `N(0, 1/2)` is
  unchanged), through `bcf_tau_0_posterior()` and `bcf_coding_suff_stats()`,
  which the count path shares. Without weights, and for probit, the
  unweighted arithmetic is kept. Recovery (n = 800, 20% of rows with 400x
  variance passed as weights, constant effect 1, 10 seeds, 10 GFR + 100
  burn-in + 200 draws): without adaptive coding the posterior-mean `tau_0`
  was 8-68 away from 1 (the weighted forest offsetting it) and is now within
  0.72; with adaptive coding the ATE RMSE fell from 0.37 (max error 0.75) to
  0.048 (max 0.08). Without adaptive coding the ATE RMSE was already 0.045
  (now 0.042).

Unchanged, bitwise against the reference: Gaussian fits without weights
(intercept, adaptive coding, both, rfx with two chains, variance forest),
MCMC-only fits sampling the treatment leaf scale, and weighted probit fits.
The count tests' Gaussian reference `gaussian_dose_weighted` carries weights,
so it is no longer compared bitwise. NASH (Gaussian, no weights, no adaptive
coding, treatment leaf scale fixed) takes none of the changed paths.

## Known upstream bugs

Pre-existing in upstream code and not fixed here:

- Probit BCF with `sample_intercept = TRUE` (the default) omits `tau_0 b(z)`
  from both the truncated-normal mean and the residual handed to the forests
  in the latent step, so the forests partly absorb `tau_0`. Fixing it changes
  probit draws. The count path includes the term.
- Probit `tau_0` and adaptive coding still ignore `observation_weights` (the
  probit latent step draws `z` with unit variance whatever the weights).
- The R-side `tau_0` and adaptive-coding updates ignore a variance forest's
  per-row variances (they use only `observation_weights`).
- `observation_weights` of zero are infinite precision under the C++
  convention: any zero weight makes Gaussian `bcf()` return `NaN` draws,
  although the input checks describe all-zero weights as prior sampling.

## Random intercept on the ZINB zero part (0.4.6.9011)

Branch `claude/zinb-zero-rfx`, on `cma-bcf-count` (d8a0cb7b), merged with
`cma-rfx` 2eab02bb (the #12 fixes). With
`rfx_group_ids_train` and `family = "zinb"`, the zero part of `bart()` and
`bcf()` gets its own random effects: `pi_i = Phi(c0 + g(x_i) [+ tau_z(x_i) z_i]
+ w_i' v_k(i))` on the count part's groups, basis and `random_effects_params`
priors, with separate parameters `v`. Default on;
`general_params$zero_random_effects = FALSE` gives the previous model draw for draw.

**Scheme**: a second `RandomEffectsModel` / `RandomEffectsDataset` (no
variance weights) / `RandomEffectsTracker` / `RandomEffectSamples`, sampled
right after the zero forest(s) on the zero part's latent residual
`z - c0 - forests - w'v` with `global_variance = 1`, the way probit `bart()`
samples its rfx. `v` enters the structural-zero draw and the truncated-normal
mean. Chain resets, warm starts and GFR deletion mirror the count rfx. Output:
`zero_rfx_samples`, `zero_rfx_preds_train/test`; `zero_hat_*` are now the full
index. `predict()` term `"zero_rfx"` (probit scale on every scale); `y_hat` on
the response / zero-probability scales includes `v` and needs group ids. JSON:
flag `has_zero_rfx` and `random_effect_container_1`. Unseen groups follow
`unseen_groups` in `bart()`; `bcf()` refuses them on both parts (no policy
there, unchanged).

**Recovery** (60 tracts x 30 rows, intercept SD 0.4 count / 0.7 zero, r = 5,
8 seeds): `bart()` pi correlation 0.882-0.927 with vs 0.535-0.677 without,
response RMSE ratio <= 0.834, intercept correlation 0.870-0.924. `bcf()` with
treatment correlated 0.35-0.49 with the zero intercept: zero-part effect error
0.076-0.296 with vs 0.556-0.824 without; the remainder is leakage through the
partially pooled intercepts (error 0.01-0.2 when treatment is unconfounded).

**Mixing** (`bcf()`, same panel, 4 chains x 1,000 draws, 3 seeds, GFR and root
starts): split R-hat of the response-scale ATE <= 1.008 with vs <= 1.005
without; bulk ESS 468-660 vs 520-851 of 4,000. Zero-part effect: R-hat <= 1.035
vs <= 1.025, ESS 133-192 vs 173-243. **Timing** (n = 5,000, 100 tracts, default
trees, ms per iteration): `bart()` ZINB 32.2 -> 32.6; `bcf()` ZINB 69-71 ->
68.5-68.8 (within run-to-run noise).

The seeded references `tests/testthat/fixtures/zinb-rfx-reference.rds` pin
ZINB without group ids, ZINB with the switch off (`bart()` and `bcf()`) and NB
with group ids (`bart()` and `bcf()`) to d8a0cb7b bitwise.

## A prior of its own for the ZINB zero-part random effects (0.4.6.9011)

Branch `claude/zinb-zero-rfx-prior-3da03e`, on `cma-rfx` (1273f128). `bart()`
and `bcf()` take `zero_random_effects_params`, a list with any of the six
prior keys of `random_effects_params` (`working_parameter_prior_mean`,
`group_parameter_prior_mean`, `working_parameter_prior_cov`,
`group_parameter_prior_cov`, `variance_prior_shape`, `variance_prior_scale`).
Keys left out take the count part's values, so the default `list()` gives the
zero part the count part's prior as before. The groups, basis (`model_spec`)
and `unseen_groups` stay shared and are refused here. A non-empty list without
`family = "zinb"`, group ids and `zero_random_effects = TRUE` is an error.

**Scheme.** R only: `.zinbZeroRfxPrior()` (`R/negbin.R`) resolves the zero
part's starting values and priors the way the count part does. They replace
the count part's at the five places the zero-part `RandomEffectsModel` is set
up or reset: the initial set-up, chain starts from GFR draws, warm starts
(with and without previous zero-part samples) and root starts. The C++ needed
no change: the zero part already has its own model object, whose priors are
set through the existing setters, and the count part's code is untouched.

**Mapping a half-Cauchy prior.** The random effects are parameter-expanded:
effect `alpha * xi_k`, `alpha ~ N(0, working_parameter_prior_cov)`,
`xi_k ~ N(0, s2)`, `s2 ~ IG(variance_prior_shape, variance_prior_scale)` with
`1 / s2 ~ Gamma(shape, rate = scale)`. The prior means are starting values
only (both conditionals have zero prior mean in `src/random_effects.cpp`).
With shape = scale = `nu / 2`, the SD `|alpha| sqrt(s2)` is half-t with `nu`
degrees of freedom and scale `sqrt(working_parameter_prior_cov)` (Gelman 2006).
So `list(working_parameter_prior_cov = A^2, variance_prior_shape = 0.5,
variance_prior_scale = 0.5)` is exactly half-Cauchy(0, A). The default (1, 1,
covariance 1) is half-t_2 with scale 1 (prior median 0.82).

**Tests** (`test-zinb-zero-rfx-prior.R`). Seeded references
(`fixtures/zinb-zero-rfx-prior-reference.rds`, from a build of 1273f128) pin,
bitwise: Gaussian, probit and NB `bart()` / `bcf()` with group ids, cloglog
`bart()`, and default-prior ZINB `bart()` / `bcf()` with zero-part random
effects (GFR-started and root-started chains, warm starts). The count prior
written out as the zero prior, in full or key by key, reproduces the default
draws bitwise, including a non-default count prior that the zero part
inherits. A tight zero prior holds the zero-part effects below 0.01 at every
chain start while the count part's stay free, and the reverse also holds.
A half-Cauchy with scale 0.001 holds strong group effects (true SD about 1),
and with scale 0.1 lets them through, as a heavy tail should.

## Upstream stochtree changes through 2026-10-05 (0.4.6.9012)

Upstream `main` was reviewed on 2026-11-12 up to 737f397c (5 Oct 2026). The
fork point is upstream 5c09ac77 (22 Jul 2026); the fork shares no git history
with upstream because its R package is generated with `cran-bootstrap.R`.
Upstream's 18 commits since then include the four NASH sampler fixes
(StochasticTree/stochtree #425-#428), which are already here.

Ported, with paths mapped from `include/` to `src/include/`:

- aef52833 (#417): remove an unused variable in `cutpoint_candidates.cpp`.
- 5c330e4d (#418): `[[noreturn]]` on `Log::Fatal`.
- b2421090 (#420): `initialize_forest_model_cpp` accepts cloglog forests
  (leaf model 4), and `initialize_forest_model_active_forest_cpp` reads the
  starting value for them. Before, `bart()` with a cloglog link used
  `init_val` uninitialized when it set the active forest's root leaves.
- 1878b0e8 (#421): `init_val` and `leaf_scale` start at 0.0.
- be249664 (#431): `tests/testthat/test-mcmc-exact-posterior.R`, unchanged.
  One tree on four points has an enumerable posterior over leaf counts; a
  100k-draw chain must match it within one percentage point, with the
  covariate in ascending and descending order.

The uninitialized cloglog value was harmless in practice on Apple silicon:
seeded draws from the outcome reference set (Gaussian, weighted, variance
forest, probit with random effects, cloglog binary and ordinal), twelve more
cloglog fits (binary and ordinal, root- and GFR-started) and an NB fit are
bitwise identical before and after.

Not ported:

- The Eigen 3.4.1 to 5.0.1 upgrade (28794fe1, 9ff40f67) and
  `Eigen::placeholders::all` (9c4dcd49, needed only under Eigen 5). It would
  replace the vendored linear algebra for no functional gain and puts
  bitwise-reproducible refits at risk.
- `.gitmodules`, Python tests, and upstream NEWS and contributor updates,
  which do not apply to the generated R package.
- Open upstream PR #419 (`std::lgamma` for `boost::math::lgamma` in
  `leaf_model.cpp`): unmerged, and it could change the last bits of the
  variance-forest marginal likelihood.
