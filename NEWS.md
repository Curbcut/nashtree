# nashtree 0.4.6.9011

* The zero part of zero-inflated count models (`family = "zinb"`) in `bart()`
  and `bcf()` gets random effects of its own whenever `rfx_group_ids_train` is
  given: `pi = Phi(c0 + g(X) [+ tau_z(X) Z] + w' v_group)`, on the count part's
  groups, basis (`model_spec`) and `random_effects_params` priors but with
  separate parameters, sampled after the zero forest(s) on the zero part's
  probit latent residual with unit variance. The intercept enters `pi`
  everywhere: the structural-zero draw, the latent step, `zero_hat_train` /
  `zero_hat_test` (now the full probit index), `predict()` (new term
  `"zero_rfx"`, on the probit scale; `y_hat` on the `"response"` and
  `"zero_probability"` scales now needs `rfx_group_ids`), the posterior
  predictive, contrasts, and JSON (a second random-effects container). Unseen
  groups follow `random_effects_params$unseen_groups` in `bart()`; `bcf()`
  refuses them on both parts, as before. New outputs `zero_rfx_samples`,
  `zero_rfx_preds_train`, `zero_rfx_preds_test`.
  `general_params$zero_random_effects = FALSE` turns it off and reproduces the
  previous draws bitwise; models without group ids, negative-binomial models
  and saved models without the field are unchanged.

* `bcf()` fits count outcomes: `OutcomeModel("count", "log")` (negative
  binomial) and `OutcomeModel("count", "log", family = "zinb")`. The count part
  is `log(mu) = offset + c + mu(X, pi_hat) + [tau_0 + tau(X)] b(Z) [+ rfx]`, so
  the treatment effect is a log rate ratio; the Pólya-Gamma step turns each
  iteration into the weighted Gaussian regression `bcf()` already fits, and the
  prognostic and treatment forests (binary or continuous dose), the propensity
  covariate, adaptive coding, `tau_0`, random effects, GFR and multiple chains
  run on it. `tau_0` and the adaptive-coding parameters use the Pólya-Gamma
  precision weights. The treatment-forest prior is set on the log-rate-ratio
  scale by the new `treatment_effect_forest_params$log_rate_ratio_prior_sd`
  (default 0.5), and `tau_0` has a `N(0, 1)` prior by default. The zero-inflated
  family adds probit zero-part prognostic and treatment forests
  (`zero_prognostic_forest_params`, `zero_treatment_effect_forest_params`). New
  arguments `offset_train` / `offset_test`; `predict.bcfmodel()` gains the
  `"rate"`, `"response"` and `"zero_probability"` scales, the `"zero_prognostic"`,
  `"zero_tau"` and `"zero_forest"` terms and an `offset`;
  `computeContrastBCFModel()` takes `offset_0` / `offset_1` and the count scales;
  `sampleBCFPosteriorPredictive()` draws counts. Variance forests and
  `sample_sigma2_global = TRUE` are refused with counts. Gaussian and probit
  `bcf()` draws are bitwise unchanged.

* `bart()` fits zero-inflated negative-binomial outcomes with
  `OutcomeModel("count", "log", family = "zinb")` (`family = "negbin"` stays the
  default and samples exactly as before). With probability
  `pi = Phi(c0 + g(X))` an observation is a structural zero, otherwise it
  follows the negative-binomial count part. Each iteration draws the
  structural-zero indicators, samples a probit zero forest `g` on them
  (configured by the new `zero_forest_params`, default 50 trees with the leaf
  scale fixed at `2 / num_trees`) and runs the Pólya-Gamma count step on the
  non-structural rows only: structural zeros get an infinite variance weight and
  drop out of the leaf and random-effects sufficient statistics. `r` is updated
  over the non-structural rows; random effects enter the count part. `predict()`
  gains `scale = "zero_probability"` (`pi`), `scale = "response"` becomes
  `(1 - pi) mu` for this family, and the term `"zero_forest"` returns the zero
  part; `sampleBARTPosteriorPredictive()` includes the zero mass. Draws of the
  probit index are returned as `zero_hat_train` / `zero_hat_test`, and the zero
  forest and family are serialized with the model.

* `bart()` fits count outcomes with a negative-binomial likelihood:
  `general_params = list(outcome_model = OutcomeModel("count", "log"))` gives
  `y ~ NB(mean mu, size r)` with `log(mu) = offset + c + f(X)` plus optional
  random effects. Each iteration draws Pólya-Gamma variables
  `omega ~ PG(y + r, log(mu) - log(r))`, which turn the likelihood into a
  Gaussian regression with variance weights `1 / omega` and unit error
  variance, so the existing forest and random-effects samplers run unchanged.
  The dispersion `r` is drawn from its exact full conditional by slice
  sampling on `log(r)` under a `Gamma(negbin_dispersion_prior_shape,
  negbin_dispersion_prior_rate)` prior (defaults 1 and 0.01). New arguments
  `offset_train` and `offset_test` take a log-scale offset (log exposure).
  `predict()` gains `scale = "rate"` (and its alias `"response"`) and an
  `offset` argument; `sampleBARTPosteriorPredictive()` draws negative-binomial
  counts and also takes `offset`. Draws of `r` are returned as
  `negbin_dispersion_samples`, serialized with the model, combined across
  chains and JSON files, and available through
  `extractParameter(model, "dispersion")`. The leaf scale defaults to the
  Chipman, George and McCulloch (2010) calibration on the range of
  `log1p(y)`. Variance forests and `sample_sigma2_global = TRUE` are refused
  with a count outcome.

* New internal Pólya-Gamma sampler (`src/include/stochtree/polya_gamma_sampler.h`),
  written from Devroye (2009), Polson, Scott and Windle (2013) and Windle,
  Polson and Scott (2014): exact for integer `b <= 50`, a 100-term gamma series
  with an exact remainder mean below `b = 10,000`, and a moment-matched normal
  above.

* The random-effects sampler honours variance weights on a
  `RandomEffectsDataset`: row `i` enters the group and working-parameter
  conditionals with precision `1 / (sigma^2 * v_i)`. Datasets without weights
  sample exactly as before; seeded Gaussian BART and BCF random-effects draws
  are bitwise identical to 0.4.6.9010 (0ac0e7d0 and the nashtree rename
  bd5c5c48 give the same draws).

# stochtree 0.4.6.9010

* Outcome summaries use read-only pointers, preserving shared ALTREP outcomes
  during fitting instead of retaining private copies in the input inventory.

# stochtree 0.4.6.9009

* Numeric covariate binding reads ordinary integer and double columns through
  read-only pointers. Shared ALTREP inputs no longer acquire retained private
  copies during training or prediction preprocessing. Matrix values, missing
  values, feature ordering and sampler behavior are unchanged.

# stochtree 0.4.6.9008

* Importance-only BART no longer allocates or returns training prediction draws
  (`y_hat_train` or `rfx_preds_train`). Sampling, split counts, variance draws,
  and retained random-effects parameters are unchanged. This removes
  observation-by-draw buffers and their prediction/transformation copies.

# stochtree 0.4.6.9007

* Random effects work with both retention modes. `forest_retention =
  "importance"` samples them as usual and keeps only the split counts, which are
  therefore conditional on the random effects, and `retained_sample_ids` keeps
  the matching random effects samples alongside the selected forests, through a
  new checked `RandomEffectSamples$retain_samples()`.

* `random_effects_params$unseen_groups` chooses what prediction does with group
  labels that were not sampled. `"error"` (the default) stops, as before.
  `"mean"` gives those rows the posterior mean random effect across the sampled
  groups, draw by draw, which is a population-average prediction for an unseen
  group. The choice is stored with the model and serialized, and `predict()`,
  the posterior predictive samplers and test-set predictions all honour it.

# stochtree 0.4.6.9003

* Add checked `ForestSamples$retain_samples()` for bulk retention of ordered
  zero-based sample IDs without changing retained forests or sampling.

# stochtree 0.4.6.9002

* Specialize the mean-tree residual update for its add/subtract operation, avoiding per-observation type-erased function calls while retaining traversal, arithmetic order and RNG consumption.

# stochtree 0.4.6.9000

# stochtree 0.4.5

## Bug Fixes

* Fix `IndexError` when sampling a mean-only Python BART model with categorical covariates; the excluded-variable weight zero-out is now guarded by the forest-inclusion flags, matching the R implementation [#411](https://github.com/StochasticTree/stochtree/pull/411)
* Fix ordinal class prediction bug in R BART [#399](https://github.com/StochasticTree/stochtree/issues/399)

# stochtree 0.4.4

## Documentation and Other Maintenance

* Switch sampling default to single-threaded execution (`num_threads = 1`) in both R and Python [#396](https://github.com/StochasticTree/stochtree/pull/396). This prevents CRAN pre-checks from failing due to limitations on core usage and also avoid unintentional default interaction with users' own multi-core setups.

# stochtree 0.4.3

## Bug Fixes

* Fix serialization / deserialization roundtrip bugs in BART and BCF [#393](https://github.com/StochasticTree/stochtree/issues/393)
* Safely handle / cast `float32` inputs in python BART / BCF [#389](https://github.com/StochasticTree/stochtree/issues/389)
* Fix ordinal cloglog probability transform for K>=4 categories [#382](https://github.com/StochasticTree/stochtree/pull/382)
* Fixed bug that separated homogeneous parametric treatment effect term (`tau_0`) from treatment effect forest predictions (`tau(x)`), only including the latter in `tau` for BCF predictions when a parametric term is included in the model [#377](https://github.com/StochasticTree/stochtree/pull/377)
* Fixed BCF prediction bug for prognostic function and mu(X) when adaptive coding is used and the tau(x) forest is not also requested by a prediction call [#377](https://github.com/StochasticTree/stochtree/pull/377)
* Fixed R BCF prediction bug when covariates are passed as dataframes and an internal propensity is sampled [#374](https://github.com/StochasticTree/stochtree/issues/374)

# stochtree 0.4.2

## New Features

* Added support for parametric treatment effect term in BCF [#309](https://github.com/StochasticTree/stochtree/pull/309/)
* Added support for observation-level weights passed as data arguments to BART and BCF [#333](https://github.com/StochasticTree/stochtree/pull/333)

## Bug Fixes

* Include boost headers in `sdist` for Python package [#362](https://github.com/StochasticTree/stochtree/pull/362)
* Fixed prediction bug for cloglog BART in R [#361](https://github.com/StochasticTree/stochtree/pull/361)
* Fixed assignment bug in RFX category tracker data structure [#360](https://github.com/StochasticTree/stochtree/pull/360)
* Fixed prediction bugs in binary classification sklearn wrapper [#359](https://github.com/StochasticTree/stochtree/pull/359)
* Fixed several initialization, sampling and prediction bugs for cloglog [#357](https://github.com/StochasticTree/stochtree/pull/357)
* Fixed num_threads pass-through bug in R BART [#339](https://github.com/StochasticTree/stochtree/pull/339)
* Fixed probit offset bug in R and Python [#337](https://github.com/StochasticTree/stochtree/pull/337)
* Standardized multi-chain BCF handling of internal propensity models [#334](https://github.com/StochasticTree/stochtree/pull/334)
* Fixed multi-chain BCF bugs with the parametric intercept term in R and Python [#326](https://github.com/StochasticTree/stochtree/pull/326)
* Fixed indexing bugs for multivariate treatment BCF in Python [#326](https://github.com/StochasticTree/stochtree/pull/326)
* Convert binary factor-valued treatments to 0/1 binary numeric treatment in `bcf()` R function [#332](https://github.com/StochasticTree/stochtree/pull/332)

## Documentation and Other Maintenance

* Add support for python 3.14 and drop 3.8 and 3.9 [#363](https://github.com/StochasticTree/stochtree/pull/363)

# stochtree 0.4.1

## Documentation and Other Maintenance

* Converted all R man page titles to title case [#310](https://github.com/StochasticTree/stochtree/pull/310)
* Updated Python `summary()` methods to print directly to console rather than return a summary string [#310](https://github.com/StochasticTree/stochtree/pull/310)
* Updated vendored pybind11 dependency to version 3.0.2 [#310](https://github.com/StochasticTree/stochtree/pull/310)

# stochtree 0.4.0

## New Features

* Added support for ordinal outcome modeling through complementary log-log link function in the BART models ([#196](https://github.com/StochasticTree/stochtree/pull/196))
* Added `__str__`, `__repr__`, `summary`, and `extract_parameter` methods in Python for the `BARTModel` and `BCFModel` classes ([#298](https://github.com/StochasticTree/stochtree/pull/298))
* Added plotting utility function (`plot_parameter_trace`) in Python that operates on both the `BARTModel` and `BCFModel` classes ([#298](https://github.com/StochasticTree/stochtree/pull/298))
* Added vignettes for summary / model inspection utilities in both R and Python ([#298](https://github.com/StochasticTree/stochtree/pull/298))
* Added `print` methods in R and `__str__` methods in Python for the forest container and random effects container objects ([#298](https://github.com/StochasticTree/stochtree/pull/298))
* Updated R documentation to group related functions into topics ([#302](https://github.com/StochasticTree/stochtree/pull/302)):
  * `BARTSerialization`
  * `BCFSerialization`
  * `ForestSamplesSerialization`
  * `RandomEffectSamplesSerialization`
  * `DataPreprocessing`
  * `ForestKernelComputation`
  * `ForestStateManagement`
  * `RandomEffectStateManagement`
* Converted the following R function names from snake case to camel case ([#302](https://github.com/StochasticTree/stochtree/pull/302)):
  * `compute_bart_posterior_interval` -> `computeBARTPosteriorInterval`
  * `compute_bcf_posterior_interval` -> `computeBCFPosteriorInterval`
  * `compute_contrast_bart_model` -> `computeContrastBARTModel`
  * `compute_contrast_bcf_model` -> `computeContrastBCFModel`
  * `extract_parameter` -> `extractParameter`
  * `sample_bart_posterior_predictive` -> `sampleBARTPosteriorPredictive`
  * `sample_bcf_posterior_predictive` -> `sampleBCFPosteriorPredictive`

## Bug Fixes

* Fixed status logging bugs for multi-chain R MCMC loops ([#298](https://github.com/StochasticTree/stochtree/pull/298))

# stochtree 0.3.1

## New Features

* Replaced C++ standard library distributions (`discrete_distribution`, `uniform_real_distribution`, `normal_distribution`, and `gamma_distribution`) with custom implementations for cross-platform reproducibility.
* Substituted custom implementations for base R `mean()`, `var()`, and `sd()` in the preprocessing logic of the R `bart()` and `bcf()` functions for enhanced numeric stability across platforms.

# stochtree 0.3.0

## New Features

* Added `print`, `summary`, `plot`, and `extract_parameter` generic functions in R for the `bartmodel` and `bcfmodel` classes ([#271](https://github.com/StochasticTree/stochtree/pull/271))
* Added sklearn-compatible estimator wrapper for `BARTModel` in Python ([#270](https://github.com/StochasticTree/stochtree/pull/270))

## Bug Fixes

* Fix R bug where our approach to temporarily modifying users' RNG state failed if `.Random.seed` did not exist (i.e. if the R RNG hadn't yet been accessed by an R session) ([#258](https://github.com/StochasticTree/stochtree/issues/258))
* Fix prediction bug for R BART models with random effects with labels that aren't straightforward `1:num_groups` integers when only `y_hat` is requested ([#256](https://github.com/StochasticTree/stochtree/pull/256))
* Fix issue with C++ standard specification in Windows R package config ([#276](https://github.com/StochasticTree/stochtree/pull/276))

# stochtree 0.2.1

## Bug Fixes

* Fix prediction bug for univariate random effects models in R ([#248](https://github.com/StochasticTree/stochtree/pull/248))
* Fix prediction bug for Python BART and BCF models with random effects with labels that aren't straightforward `0:(num_groups-1)` integers ([#256](https://github.com/StochasticTree/stochtree/pull/256))

## Other Changes

* Encode expectations about which combinations of BART / BCF features work together and ensure warning ([#250](https://github.com/StochasticTree/stochtree/pull/250))

# stochtree 0.2.0

## New Features

* Support for multithreading in various elements of the GFR and MCMC algorithms ([#182](https://github.com/StochasticTree/stochtree/pull/182))
* Support for binary outcomes in BART and BCF with a probit link ([#164](https://github.com/StochasticTree/stochtree/pull/164))
* Enable "restricted sweep" of tree algorithms over a handful of trees ([#173](https://github.com/StochasticTree/stochtree/pull/173))
* Support for multivariate treatment in R ([#183](https://github.com/StochasticTree/stochtree/pull/183))
* Enable modification of dataset variables (weights, etc...) via low-level interface ([#194](https://github.com/StochasticTree/stochtree/pull/194))

## Computational Improvements

* Modified default random effects initialization ([#190](https://github.com/StochasticTree/stochtree/pull/190))
* Avoid double prediction on training set ([#178](https://github.com/StochasticTree/stochtree/pull/178))

## Bug Fixes

* Fixed indexing bug in cleanup of grow-from-root (GFR) samples in BART and BCF models
* Avoid using covariate preprocessor in `computeForestLeafIndices` function when a `ForestSamples` object is provided (rather than a `bartmodel` or `bcfmodel` object)
* Correctly compute feature-specific split counts in R and Python ([#220](https://github.com/StochasticTree/stochtree/issues/220))
* Avoid override of user-specified `num_burnin` parameter in BCF models with an internal propensity score ([#222](https://github.com/StochasticTree/stochtree/issues/222))
* Outcome predictions correctly incorporate adaptive coding of untreated observations in BCF with binary treatment ([#231](https://github.com/StochasticTree/stochtree/issues/231))

## Documentation Improvements

* Clarify structure / layout of samples when users request multiple chains in BART and BCF models ([#220](https://github.com/StochasticTree/stochtree/issues/220))

## Other Changes

* Standardized naming conventions for data elements of BART and BCF models across R and Python interfaces
    * Covariates / features are always referred to as "`X`"
    * Treatment is always referred to as "`Z`"
    * Propensity scores are referred to as "`propensity`" (rather than "`pi`")
    * Outcomes are referred to as "`y`"
    * Basis vectors for leaf-wise regression models in forest terms are referred to as "`leaf_basis`"
    * Group labels for additive random effects models are referred to as "`rfx_group_ids`"
    * Basis vectors for additive random effects models are referred to as "`rfx_basis`"
* Run-time checks for variables that are treated as continuous but have many "ties" (which presents issues with the current GFR algorithm) when only GFR samples are requested ([#243](https://github.com/StochasticTree/stochtree/pull/243))

# stochtree 0.1.1

* Fixed initialization bug in several R package code examples for random effects models

# stochtree 0.1.0

* Initial release on CRAN.
* Support for sampling stochastic tree ensembles using two algorithms: MCMC and Grow-From-Root (GFR)
* High-level model types supported:
    * Supervised learning with constant leaves or user-specified leaf regression models
    * Causal effect estimation with binary or continuous treatments
* Additional high-level modeling features:
    * Forest-based variance function estimation (heteroskedasticity)
    * Additive (univariate or multivariate) group random effects
    * Multi-chain sampling and support for parallelism
    * "Warm-start" initialization of MCMC forest samplers via the Grow-From-Root (GFR) algorithm
    * Automated preprocessing / handling of categorical variables
* Low-level interface:
    * Ability to combine a forest sampler with other (additive) model terms, without using C++
    * Combine and sample an arbitrary number of forests or random effects terms
