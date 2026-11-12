#' BART Serialization Routines
#' @name BARTSerialization
#' @description
#' BART models contains external pointers to C++ objects, which means they cannot
#' be correctly serialized to `.Rds` from an R session in their default state.
#' This group of serialization functions allow us to convert between C++ data structures and a persistent JSON
#' representation. The `CppJson` class wraps a performant C++ JSON API, and the functions
#' `saveBARTModelToJson` and `createBARTModelFromJson` save to and load from this format.
#' This representation, of course, also relies on external C++ pointers, so in order to
#' save and reload BART models across sessions, we provide two other interfaces.
#'
#' `saveBARTModelToJsonString` converts a BART model to an in-memory string containing the model's
#' JSON representation and `createBARTModelFromJsonString` converts this representation back to a BART model object.
#'
#' `saveBARTModelToJsonFile` and `createBARTModelFromJsonFile` save or reload a BART model
#' directly to / from a `.json` file.
#'
#' Finally, for cases in which multiple BART models have been sampled (for instance, multiple processes
#' run via `doParallel`), we offer `createBARTModelFromCombinedJson` and `createBARTModelFromCombinedJsonString` for
#' loading a new combined BART model from a list of BART JSON objects / strings.
#' @returns
#' `saveBARTModelToJson` return an object of type `CppJson`.
#' `saveBARTModelToJsonString` returns a string dump of the BART model's JSON representation.
#' `saveBARTModelToJsonFile` returns nothing, but writes to the provided filename.
#'
#' `createBARTModelFromJson`, `createBARTModelFromJsonFile`, `createBARTModelFromJsonString`,
#' `createBARTModelFromCombinedJson`, and `createBARTModelFromCombinedJsonString` all return
#' objects of type `bartmodel`.
#' @examples
#' # Generate data
#' n <- 100
#' p <- 5
#' X <- matrix(runif(n*p), ncol = p)
#' y <- X[,1] + rnorm(n, 0, 1)
#'
#' # Sample BART model
#' bart_model <- bart(X_train = X, y_train = y, num_gfr = 0,
#'                    num_burnin = 0, num_mcmc = 10)
#'
#' # Save to in-memory JSON
#' bart_json <- saveBARTModelToJson(bart_model)
#' # Save to JSON string
#' bart_json_string <- saveBARTModelToJsonString(bart_model)
#' # Save to JSON file
#' tmpjson <- tempfile(fileext = ".json")
#' saveBARTModelToJsonFile(bart_model, file.path(tmpjson))
#'
#' # Reload BART model from in-memory JSON object
#' bart_model_roundtrip <- createBARTModelFromJson(bart_json)
#' # Reload BART model from JSON string
#' bart_model_roundtrip <- createBARTModelFromJsonString(bart_json_string)
#' # Reload BART model from JSON file
#' bart_model_roundtrip <- createBARTModelFromJsonFile(file.path(tmpjson))
#' unlink(tmpjson)
#' # Reload BART model from list of JSON objects
#' bart_model_roundtrip <- createBARTModelFromCombinedJson(list(bart_json))
#' # Reload BART model from list of JSON strings
#' bart_model_roundtrip <- createBARTModelFromCombinedJsonString(list(bart_json_string))
#'
NULL
#> NULL

#' @title Run BART for Supervised Learning
#' @description
#' Run the BART algorithm for supervised learning.
#'
#' @param X_train Covariates used to split trees in the ensemble. May be provided either as a dataframe or a matrix.
#' Matrix covariates will be assumed to be all numeric. Covariates passed as a dataframe will be
#' preprocessed based on the variable types (e.g. categorical columns stored as unordered factors will be one-hot encoded,
#' categorical columns stored as ordered factors will passed as integers to the core algorithm, along with the metadata
#' that the column is ordered categorical).
#' @param y_train Outcome to be modeled by the ensemble.
#' @param leaf_basis_train (Optional) Bases used to define a regression model `y ~ W` in
#' each leaf of each regression tree. By default, BART assumes constant leaf node
#' parameters, implicitly regressing on a constant basis of ones (i.e. `y ~ 1`).
#' @param rfx_group_ids_train (Optional) Group labels used for an additive random effects model.
#' @param rfx_basis_train (Optional) Basis for "random-slope" regression in an additive random effects model.
#' If `rfx_group_ids_train` is provided with a regression basis, an intercept-only random effects model
#' will be estimated.
#' @param X_test (Optional) Test set of covariates used to define "out of sample" evaluation data.
#' May be provided either as a dataframe or a matrix, but the format of `X_test` must be consistent with
#' that of `X_train`.
#' @param leaf_basis_test (Optional) Test set of bases used to define "out of sample" evaluation data.
#' While a test set is optional, the structure of any provided test set must match that
#' of the training set (i.e. if both `X_train` and `leaf_basis_train` are provided, then a test set must
#' consist of `X_test` and `leaf_basis_test` with the same number of columns).
#' @param rfx_group_ids_test (Optional) Test set group labels used for an additive random effects model.
#' We do not currently support (but plan to in the near future), test set evaluation for group labels
#' that were not in the training set.
#' @param rfx_basis_test (Optional) Test set basis for "random-slope" regression in additive random effects model.
#' @param observation_weights (Optional) Numeric vector of observation weights of length `nrow(X_train)`. Weights are
#'   applied as `y_i | - ~ N(mu(X_i), sigma^2 / w_i)`, so larger weights increase an observation's influence on the fit.
#'   All weights must be non-negative. Default: `NULL` (all observations equally weighted). Compatible with Gaussian
#'   (continuous/identity) and probit outcome models; not compatible with cloglog link functions. Note: these are
#'   referred to internally in the C++ layer as "variance weights" (`var_weights`), since they scale the residual variance.
#' @param num_gfr Number of "warm-start" iterations run using the grow-from-root algorithm (He and Hahn, 2021). Default: 5.
#' @param num_burnin Number of "burn-in" iterations of the MCMC sampler. Default: 0.
#' @param num_mcmc Number of "retained" iterations of the MCMC sampler. Default: 100.
#' @param retained_sample_ids Optional strictly increasing unique zero-based posterior
#' sample IDs. Sampling still runs every iteration. Returns a `bartselected` bundle:
#' `model` is an ordinary BART model containing only these forests and matching
#' training/variance samples; the outer bundle retains full training/variance draws
#' and full test means when `X_test` is supplied. Use
#' [sampleBARTSelectedPosteriorPredictive()] for all-iteration test scoring.
#' Only single-chain Gaussian constant-leaf mean forests without retained warmup,
#' variance forests or previous-model initialization are supported. Random
#' effects are supported: the retained samples of the random effects container
#' are kept alongside the selected forests. Cannot be combined with importance
#' retention. Default `NULL` keeps existing behavior.
#' @param offset_train (Optional) Count outcomes only. Numeric vector of length `nrow(X_train)` added to the
#'   linear predictor on the log scale, `log(mu_i) = offset_train[i] + f(X_i)`, typically the log of an
#'   exposure such as population or time at risk. Default: `NULL` (no offset).
#' @param zero_forest_params (Optional) Zero-inflated count outcomes only (`OutcomeModel("count", "log", family = "zinb")`). A list of zero-forest parameters, each with a default:
#'
#'   - `num_trees` Number of trees in the probit zero forest. Default: `50`.
#'   - `alpha`, `beta`, `min_samples_leaf`, `max_depth` Tree prior and size limits, as for the mean forest. Defaults: `0.95`, `2`, `5`, `10`.
#'   - `sample_sigma2_leaf` Whether to sample the zero-forest leaf scale. Default: `FALSE` (fixed at the probit calibration, as in Chipman, George and McCulloch 2010).
#'   - `sigma2_leaf_init` Zero-forest leaf scale. Default: `2 / num_trees`.
#'   - `sigma2_leaf_shape`, `sigma2_leaf_scale` Inverse-gamma prior of the leaf scale when it is sampled. Defaults: `3` and `1 / num_trees`.
#'   - `keep_vars`, `drop_vars`, `num_features_subsample` Variable selection for the zero forest, as for the mean forest.
#'
#'   Supplying `zero_forest_params` for any other outcome model is an error.
#' @param zero_random_effects_params (Optional) Zero-inflated count outcomes (`family = "zinb"`) with
#'   `rfx_group_ids_train` and `general_params$zero_random_effects = TRUE` only. The prior of the zero
#'   part's random effects, as a list with any of the prior keys of `random_effects_params`:
#'   `working_parameter_prior_mean`, `group_parameter_prior_mean`, `working_parameter_prior_cov`,
#'   `group_parameter_prior_cov`, `variance_prior_shape` and `variance_prior_scale`. A key left out takes the
#'   count part's value from `random_effects_params`, so the default `list()` gives both parts the same prior,
#'   as before this argument existed (draw for draw). The groups, basis (`model_spec`) and unseen-group policy
#'   are always the count part's and cannot be set here. Supplying a non-empty list for any other model is an
#'   error.
#'
#'   The random effects are parameter-expanded: a group's effect is `alpha * xi_k` (per basis component),
#'   with the working parameter `alpha ~ N(0, working_parameter_prior_cov)`, `xi_k ~ N(0, s2)` and
#'   `s2 ~ IG(variance_prior_shape, variance_prior_scale)` (shape, scale; `1 / s2 ~ Gamma(shape, rate =
#'   scale)`). The prior means only set the starting values. With `variance_prior_shape =
#'   variance_prior_scale = nu / 2` the implied prior on the effect's standard deviation `|alpha| sqrt(s2)`
#'   is a half-t with `nu` degrees of freedom and scale `sqrt(working_parameter_prior_cov)` (Gelman 2006);
#'   `nu = 1` (shape and scale `0.5`) is a half-Cauchy with that scale. The defaults (`1`, `1`, covariance 1)
#'   are a half-t with 2 degrees of freedom and scale 1 on the probit scale. For example,
#'   `list(working_parameter_prior_cov = 0.094^2, variance_prior_shape = 0.5, variance_prior_scale = 0.5)`
#'   puts a half-Cauchy(0, 0.094) prior on the SD of a zero-part random intercept.
#' @param offset_test (Optional) Count outcomes only. Offset for `X_test`, of length `nrow(X_test)`. If
#'   `offset_train` is supplied and `offset_test` is not, test-set predictions use an offset of 0 (unit
#'   exposure) and a warning is raised.
#' @param forest_retention Forest output policy: `"all"` (default) retains the usual model.
#' `"importance"` accumulates aggregate mean-forest split counts across exactly
#' the retained iterations, without retaining posterior forests. It returns a
#' `bartimportance` summary with `split_counts`, variance draws and metadata.
#' Training predictions (`y_hat_train` and `rfx_preds_train`) are neither
#' allocated nor returned in this mode. Counts follow the processed covariate ordering.
#' All sampling, burn-in and thinning transitions still run. GFR forests are temporarily kept
#' for chain initialization. This mode currently supports single-chain Gaussian
#' constant-leaf mean forests without test data, variance forests or
#' previous-model initialization. Random effects are supported: they are sampled
#' as usual, so the split counts are conditional on them. It cannot predict new
#' data, serve as a warm-start model or be saved with the BART model JSON
#' functions.
#' @param previous_model_json (Optional) JSON string containing a previous BART model. This can be used to "continue" a sampler interactively after inspecting the samples or to run parallel chains "warm-started" from existing forest samples. Default: `NULL`.
#' @param previous_model_warmstart_sample_num (Optional) Sample number from `previous_model_json` that will be used to warmstart this BART sampler. One-indexed (so that the first sample is used for warm-start by setting `previous_model_warmstart_sample_num = 1`). Default: `NULL`.  If `num_chains` in the `general_params` list is > 1, then each successive chain will be initialized from a different sample, counting backwards from `previous_model_warmstart_sample_num`. That is, if `previous_model_warmstart_sample_num = 10` and `num_chains = 4`, then chain 1 will be initialized from sample 10, chain 2 from sample 9, chain 3 from sample 8, and chain 4 from sample 7. If `previous_model_json` is provided but `previous_model_warmstart_sample_num` is NULL, the last sample in the previous model will be used to initialize the first chain, counting backwards as noted before. If more chains are requested than there are samples in `previous_model_json`, a warning will be raised and only the last sample will be used.
#' @param general_params (Optional) A list of general (non-forest-specific) model parameters, each of which has a default value processed internally, so this argument list is optional.
#'
#'   - `cutpoint_grid_size` Maximum size of the "grid" of potential cutpoints to consider in the GFR algorithm. Default: `100`.
#'   - `standardize` Whether or not to standardize the outcome (and store the offset / scale in the model object). Default: `TRUE`.
#'   - `sample_sigma2_global` Whether or not to update the `sigma^2` global error variance parameter based on `IG(sigma2_global_shape, sigma2_global_scale)`. Default: `TRUE`.
#'   - `sigma2_global_init` Starting value of global error variance parameter. Calibrated internally as `1.0*var(y_train)`, where `y_train` is the possibly standardized outcome, if not set.
#'   - `sigma2_global_shape` Shape parameter in the `IG(sigma2_global_shape, sigma2_global_scale)` global error variance model. Default: `0`.
#'   - `sigma2_global_scale` Scale parameter in the `IG(sigma2_global_shape, sigma2_global_scale)` global error variance model. Default: `0`.
#'   - `variable_weights` Numeric weights reflecting the relative probability of splitting on each variable. Does not need to sum to 1 but cannot be negative. Defaults to `rep(1/ncol(X_train), ncol(X_train))` if not set here. Note that if the propensity score is included as a covariate in either forest, its weight will default to `1/ncol(X_train)`.
#'   - `random_seed` Integer parameterizing the C++ random number generator. If not specified, the C++ random number generator is seeded according to `std::random_device`.
#'   - `keep_burnin` Whether or not "burnin" samples should be included in the stored samples of forests and other parameters. Default `FALSE`. Ignored if `num_mcmc = 0`.
#'   - `keep_gfr` Whether or not "grow-from-root" samples should be included in the stored samples of forests and other parameters. Default `FALSE`. Ignored if `num_mcmc = 0`.
#'   - `keep_every` How many iterations of the burned-in MCMC sampler should be run before forests and parameters are retained. Default `1`. Setting `keep_every <- k` for some `k > 1` will "thin" the MCMC samples by retaining every `k`-th sample, rather than simply every sample. This can reduce the autocorrelation of the MCMC samples.
#'   - `num_chains` How many independent MCMC chains should be sampled. If `num_mcmc = 0`, this is ignored. If `num_gfr = 0`, then each chain is run from root for `num_mcmc * keep_every + num_burnin` iterations, with `num_mcmc` samples retained. If `num_gfr > 0`, each MCMC chain will be initialized from a separate GFR ensemble, with the requirement that `num_gfr >= num_chains`. Default: `1`. Note that if `num_chains > 1`, the returned model object will contain samples from all chains, stored consecutively. That is, if there are 4 chains with 100 samples each, the first 100 samples will be from chain 1, the next 100 samples will be from chain 2, etc... For more detail on working with multi-chain BART models, see [the multi chain vignette](https://stochtree.ai/vignettes/multi-chain.html).
#'   - `verbose` Whether or not to print progress during the sampling loops. Default: `FALSE`.
#'   - `outcome_model` A structured `OutcomeModel` object that specifies the outcome type and desired link function. This argument pre-empts the legacy (deprecated) `probit_outcome_model` option. Default: `OutcomeModel(outcome='continuous', link='identity')`. `OutcomeModel(outcome = 'count', link = 'log')` fits a negative-binomial model, see "Count outcomes" below.
#'   - `negbin_dispersion_prior_shape` Count outcomes only. Shape `a0` of the `Gamma(a0, b0)` (shape, rate) prior on the negative-binomial size (dispersion) parameter `r`. Default: `1`.
#'   - `negbin_dispersion_prior_rate` Count outcomes only. Rate `b0` of the same prior. Default: `0.01`, so the default prior is exponential with mean 100 and allows near-Poisson data.
#'   - `negbin_dispersion_init` Count outcomes only. Starting value of `r`. Defaults to a method-of-moments estimate around the constant-rate fit, clamped to `[0.1, 1000]`.
#'   - `sample_negbin_dispersion` Count outcomes only. Whether to sample `r`. If `FALSE`, `r` stays at `negbin_dispersion_init`. Default: `TRUE`.
#'   - `zero_random_effects` Zero-inflated count outcomes with `rfx_group_ids_train` only. Whether the zero part gets its own random effects on the same groups and basis as the count part (see "Zero-inflated count outcomes"). `FALSE` gives the model without them, draw for draw. Default: `TRUE`.
#'   - `probit_outcome_model` Deprecated in favor of `outcome_model`. Whether or not the outcome should be modeled as explicitly binary via a probit link. If `TRUE`, `y` must only contain the values `0` and `1`. Default: `FALSE`.
#'   - `num_threads` Number of threads to use in the GFR and MCMC algorithms, as well as prediction. Defaults to `1` (single-threaded). Set to `-1` to use the maximum number of available threads, or a positive integer for a specific count. OpenMP must be available for values other than `1`.
#'
#' @param mean_forest_params (Optional) A list of mean forest model parameters, each of which has a default value processed internally, so this argument list is optional.
#'
#'   - `num_trees` Number of trees in the ensemble for the conditional mean model. Default: `200`. If `num_trees = 0`, the conditional mean will not be modeled using a forest, and the function will only proceed if `num_trees > 0` for the variance forest.
#'   - `alpha` Prior probability of splitting for a tree of depth 0 in the mean model. Tree split prior combines `alpha` and `beta` via `alpha*(1+node_depth)^-beta`. Default: `0.95`.
#'   - `beta` Exponent that decreases split probabilities for nodes of depth > 0 in the mean model. Tree split prior combines `alpha` and `beta` via `alpha*(1+node_depth)^-beta`. Default: `2`.
#'   - `min_samples_leaf` Minimum allowable size of a leaf, in terms of training samples, in the mean model. Default: `5`.
#'   - `max_depth` Maximum depth of any tree in the ensemble in the mean model. Default: `10`. Can be overridden with ``-1`` which does not enforce any depth limits on trees.
#'   - `sample_sigma2_leaf` Whether or not to update the leaf scale variance parameter based on `IG(sigma2_leaf_shape, sigma2_leaf_scale)`. Cannot (currently) be set to true if `ncol(leaf_basis_train)>1`. Default: `FALSE`.
#'   - `sigma2_leaf_init` Starting value of leaf node scale parameter. Calibrated internally as `1/num_trees` if not set here.
#'   - `sigma2_leaf_shape` Shape parameter in the `IG(sigma2_leaf_shape, sigma2_leaf_scale)` leaf node parameter variance model. Default: `3`.
#'   - `sigma2_leaf_scale` Scale parameter in the `IG(sigma2_leaf_shape, sigma2_leaf_scale)` leaf node parameter variance model. Calibrated internally as `0.5/num_trees` if not set here.
#'   - `keep_vars` Vector of variable names or column indices denoting variables that should be included in the forest. Default: `NULL`.
#'   - `drop_vars` Vector of variable names or column indices denoting variables that should be excluded from the forest. Default: `NULL`. If both `drop_vars` and `keep_vars` are set, `drop_vars` will be ignored.
#'   - `num_features_subsample` How many features to subsample when growing each tree for the GFR algorithm. Defaults to the number of features in the training dataset.
#'
#' @param variance_forest_params (Optional) A list of variance forest model parameters, each of which has a default value processed internally, so this argument list is optional.
#'
#'   - `num_trees` Number of trees in the ensemble for the conditional variance model. Default: `0`. Variance is only modeled using a tree / forest if `num_trees > 0`.
#'   - `alpha` Prior probability of splitting for a tree of depth 0 in the variance model. Tree split prior combines `alpha` and `beta` via `alpha*(1+node_depth)^-beta`. Default: `0.95`.
#'   - `beta` Exponent that decreases split probabilities for nodes of depth > 0 in the variance model. Tree split prior combines `alpha` and `beta` via `alpha*(1+node_depth)^-beta`. Default: `2`.
#'   - `min_samples_leaf` Minimum allowable size of a leaf, in terms of training samples, in the variance model. Default: `5`.
#'   - `max_depth` Maximum depth of any tree in the ensemble in the variance model. Default: `10`. Can be overridden with ``-1`` which does not enforce any depth limits on trees.
#'   - `leaf_prior_calibration_param` Hyperparameter used to calibrate the `IG(var_forest_prior_shape, var_forest_prior_scale)` conditional error variance model. If `var_forest_prior_shape` and `var_forest_prior_scale` are not set below, this calibration parameter is used to set these values to `num_trees / leaf_prior_calibration_param^2 + 0.5` and `num_trees / leaf_prior_calibration_param^2`, respectively. Default: `1.5`.
#'   - `var_forest_leaf_init` Starting value of root forest prediction in conditional (heteroskedastic) error variance model. Calibrated internally as `log(0.6*var(y_train))/num_trees`, where `y_train` is the possibly standardized outcome, if not set.
#'   - `var_forest_prior_shape` Shape parameter in the `IG(var_forest_prior_shape, var_forest_prior_scale)` conditional error variance model (which is only sampled if `num_trees > 0`). Calibrated internally as `num_trees / leaf_prior_calibration_param^2 + 0.5` if not set.
#'   - `var_forest_prior_scale` Scale parameter in the `IG(var_forest_prior_shape, var_forest_prior_scale)` conditional error variance model (which is only sampled if `num_trees > 0`). Calibrated internally as `num_trees / leaf_prior_calibration_param^2` if not set.
#'   - `keep_vars` Vector of variable names or column indices denoting variables that should be included in the forest. Default: `NULL`.
#'   - `drop_vars` Vector of variable names or column indices denoting variables that should be excluded from the forest. Default: `NULL`. If both `drop_vars` and `keep_vars` are set, `drop_vars` will be ignored.
#'   - `num_features_subsample` How many features to subsample when growing each tree for the GFR algorithm. Defaults to the number of features in the training dataset.
#'
#' @param random_effects_params (Optional) A list of random effects model parameters, each of which has a default value processed internally, so this argument list is optional.
#'
#'   - `model_spec` Specification of the random effects model. Options are "custom" and "intercept_only". If "custom" is specified, then a user-provided basis must be passed through `rfx_basis_train`. If "intercept_only" is specified, a random effects basis of all ones will be dispatched internally at sampling and prediction time. If "intercept_plus_treatment" is specified, a random effects basis that combines an "intercept" basis of all ones with the treatment variable (`Z_train`) will be dispatched internally at sampling and prediction time. Default: "custom". If "intercept_only" is specified, `rfx_basis_train` and `rfx_basis_test` (if provided) will be ignored.
#'   - `working_parameter_prior_mean` Prior mean for the random effects "working parameter". Default: `NULL`. Must be a vector whose dimension matches the number of random effects bases, or a scalar value that will be expanded to a vector.
#'   - `group_parameters_prior_mean` Prior mean for the random effects "group parameters." Default: `NULL`. Must be a vector whose dimension matches the number of random effects bases, or a scalar value that will be expanded to a vector.
#'   - `working_parameter_prior_cov` Prior covariance matrix for the random effects "working parameter." Default: `NULL`. Must be a square matrix whose dimension matches the number of random effects bases, or a scalar value that will be expanded to a diagonal matrix.
#'   - `group_parameter_prior_cov` Prior covariance matrix for the random effects "group parameters." Default: `NULL`. Must be a square matrix whose dimension matches the number of random effects bases, or a scalar value that will be expanded to a diagonal matrix.
#'   - `variance_prior_shape` Shape parameter for the inverse gamma prior on the variance of the random effects "group parameter." Default: `1`.
#'   - `variance_prior_scale` Scale parameter for the inverse gamma prior on the variance of the random effects "group parameter." Default: `1`.
#'   - `unseen_groups` How prediction handles group labels that were not present in `rfx_group_ids_train`. `"error"` (the default) stops, as in previous versions. `"mean"` gives those rows the posterior mean random effect across the sampled groups, draw by draw, which is a population-average prediction for an unseen group. The choice is stored with the model and honoured by `predict()`, the posterior predictive samplers and test-set predictions.
#'
#' @section Count outcomes:
#' With `general_params = list(outcome_model = OutcomeModel("count", "log"))`, `y_train` must hold
#' non-negative integers and the model is
#' `y_i ~ NB(mean mu_i, size r)`, `log(mu_i) = offset_i + c + f(X_i) [+ random effects]`,
#' so that `Var(y_i) = mu_i + mu_i^2 / r`. `c = log(sum(y) / sum(exp(offset)))` is a fixed centering
#' constant (stored as `model_params$outcome_mean`) and `f` is the mean forest. Each iteration draws
#' Pólya-Gamma variables `omega_i ~ PG(y_i + r, log(mu_i) - log(r))` (Polson, Scott and Windle 2013), which turn
#' the likelihood into a Gaussian regression of `(y_i - r) / (2 omega_i)` with variance `1 / omega_i`; the
#' forest and random effects are then sampled by the usual Gaussian samplers with variance weights
#' `1 / omega_i` and the global error variance fixed at 1. `r` is drawn from its full conditional given
#' the current means by slice sampling on `log(r)`. A variance forest cannot be combined with a count
#' outcome, and `sample_sigma2_global = TRUE` is refused. `observation_weights`, if supplied, multiply the
#' per-iteration variance weights.
#'
#' Unless `sigma2_leaf_init` is set, the mean-forest leaf scale is calibrated on the log scale as in
#' Chipman, George and McCulloch (2010): `sigma2_leaf_init = (R / (2 k sqrt(num_trees)))^2` with
#' `k = 2` and `R` the range of `log1p(y_train)`.
#'
#' `y_hat_train` and `y_hat_test` hold draws of the linear predictor `log(mu)` including the offset, and
#' the draws of `r` are returned as `negbin_dispersion_samples`. Use `predict(..., scale = "rate")` for
#' `mu` and [sampleBARTPosteriorPredictive()] for counts.
#'
#' @section Zero-inflated count outcomes:
#' With `OutcomeModel("count", "log", family = "zinb")` the model is a zero-inflated negative binomial:
#' with probability `pi_i = Phi(c0 + g(X_i))` the observation is a structural zero, and otherwise
#' `y_i ~ NB(mu_i, r)` with `mu_i` as above. `g` is a second forest (the "zero forest", configured
#' through `zero_forest_params`) and `c0` a fixed centering constant (`model_params$zero_outcome_mean`).
#' With the structural-zero indicator `s_i` as a latent variable, each iteration
#' (1) draws `s_i ~ Bernoulli(pi_i / (pi_i + (1 - pi_i) NB(0 | mu_i, r)))` for `y_i = 0` (and sets
#' `s_i = 0` for `y_i > 0`); (2) samples the zero forest on `s` by the probit (Albert and Chib 1993)
#' latent-variable step, with its own forest, tracker and model configuration; (3) samples the count
#' forest by the Pólya-Gamma step on the rows with `s_i = 0` only: rows with `s_i = 1` carry an infinite
#' variance weight, so `1 / w_i = 0` and they add nothing to the leaf sufficient statistics, to the
#' leaf marginal likelihoods or to the random-effects conditional; a leaf that holds only such rows
#' has a finite marginal likelihood equal to its prior and its value is drawn from the leaf prior;
#' (4) draws `r` from its conditional over the rows with `s_i = 0`; and (5) samples the count part's
#' random effects.
#'
#' With `rfx_group_ids_train`, the zero part also gets random effects (unless
#' `general_params$zero_random_effects = FALSE`): `pi_i = Phi(c0 + g(X_i) + w_i' v_{k(i)})`, with the
#' basis `w_i` of the count part (a column of ones for `model_spec = "intercept_only"`), the same groups
#' and, unless `zero_random_effects_params` sets its own, the same `random_effects_params` priors, but
#' parameters `v` of their own. They are sampled right
#' after the zero forest, on the zero part's latent residual `z_i - c0 - g(X_i)` with unit variance (no
#' Pólya-Gamma weights), and enter `pi` everywhere: the draw of `s`, the latent step, `zero_hat_train`
#' and `zero_hat_test`, `predict()`, the posterior predictive and the serialized model (a second
#' random-effects container). Groups not seen in training follow `random_effects_params$unseen_groups`
#' on both parts.
#'
#' `min_samples_leaf` of the count forest counts all rows of a leaf, structural zeros included, so a
#' count-forest leaf may hold fewer than `min_samples_leaf` rows that inform it. The bound is a guard
#' against very small leaves, not a likelihood quantity, and this keeps the tree prior independent of
#' the latent `s`; the leaf prior still regularizes a leaf with few informative rows.
#'
#' The sampler is initialised from a constant-rate ZINB fit: the excess-zero fraction over the
#' constant-rate negative-binomial fit (method of moments) starts a maximum-likelihood fit of
#' `(pi, c, r)` without covariates; `c` centres the count forest, `qnorm(pi)` (clamped to
#' `[0.01, 0.99]`) centres the zero forest and `r` is the default `negbin_dispersion_init`. The first
#' draw of `s` is made from the root forests and is therefore the posterior of `s` under that constant
#' fit. The zero-forest leaf scale follows the probit calibration of Chipman, George and McCulloch (2010),
#' `sigma2_leaf_init = 2 / num_trees`, and is held fixed unless `zero_forest_params$sample_sigma2_leaf = TRUE`.
#' `observation_weights` are not supported with this family.
#'
#' The output adds `zero_forests`, `zero_hat_train` and `zero_hat_test` (draws of the full probit index
#' `c0 + g(X)` plus any zero-part random effects, so `pnorm(zero_hat_train)` are draws of `pi`) and, with
#' zero-part random effects, `zero_rfx_samples`, `zero_rfx_preds_train` and `zero_rfx_preds_test`.
#' `predict()` offers the scales `"linear"` (`log(mu)`), `"rate"` (`mu`), `"zero_probability"` (`pi`) and
#' `"response"` (`(1 - pi) mu`, the mean of `y`), and the terms `"zero_forest"` (`c0 + g(X)`) and
#' `"zero_rfx"` (the zero-part random effects, on the probit scale on every prediction scale).
#' [sampleBARTPosteriorPredictive()] includes the zero mass.
#'
#' @return List of sampling outputs and a wrapper around the sampled forests (which can be used for in-memory prediction on new data, or serialized to JSON on disk).
#' @export
#'
#' @examples
#' n <- 100
#' p <- 5
#' X <- matrix(runif(n*p), ncol = p)
#' f_XW <- (
#'     ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) +
#'     ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) +
#'     ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) +
#'     ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
#' )
#' noise_sd <- 1
#' y <- f_XW + rnorm(n, 0, noise_sd)
#' test_set_pct <- 0.2
#' n_test <- round(test_set_pct*n)
#' n_train <- n - n_test
#' test_inds <- sort(sample(1:n, n_test, replace = FALSE))
#' train_inds <- (1:n)[!((1:n) %in% test_inds)]
#' X_test <- X[test_inds,]
#' X_train <- X[train_inds,]
#' y_test <- y[test_inds]
#' y_train <- y[train_inds]
#'
#' bart_model <- bart(X_train = X_train, y_train = y_train, X_test = X_test,
#'                    num_gfr = 10, num_burnin = 0, num_mcmc = 10)
bart <- function(
  X_train,
  y_train,
  leaf_basis_train = NULL,
  rfx_group_ids_train = NULL,
  rfx_basis_train = NULL,
  X_test = NULL,
  leaf_basis_test = NULL,
  rfx_group_ids_test = NULL,
  rfx_basis_test = NULL,
  observation_weights = NULL,
  num_gfr = 5,
  num_burnin = 0,
  num_mcmc = 100,
  previous_model_json = NULL,
  previous_model_warmstart_sample_num = NULL,
  general_params = list(),
  mean_forest_params = list(),
  variance_forest_params = list(),
  random_effects_params = list(),
  forest_retention = c("all", "importance"),
  retained_sample_ids = NULL,
  offset_train = NULL,
  offset_test = NULL,
  zero_forest_params = list(),
  zero_random_effects_params = list()
) {
  forest_retention <- match.arg(forest_retention)
  importance_only <- identical(forest_retention, "importance")
  selected_only <- !is.null(retained_sample_ids)
  # Update general BART parameters
  general_params_default <- list(
    cutpoint_grid_size = 100,
    standardize = TRUE,
    sample_sigma2_global = TRUE,
    sigma2_global_init = NULL,
    sigma2_global_shape = 0,
    sigma2_global_scale = 0,
    variable_weights = NULL,
    random_seed = -1,
    keep_burnin = FALSE,
    keep_gfr = FALSE,
    keep_every = 1,
    num_chains = 1,
    verbose = FALSE,
    outcome_model = OutcomeModel(outcome = "continuous", link = "identity"),
    probit_outcome_model = FALSE,
    num_threads = 1,
    negbin_dispersion_prior_shape = 1,
    negbin_dispersion_prior_rate = 0.01,
    negbin_dispersion_init = NULL,
    sample_negbin_dispersion = TRUE,
    zero_random_effects = TRUE
  )
  general_params_updated <- preprocessParams(
    general_params_default,
    general_params
  )
  # TODO: think about validation and deprecation flow for probit_outcome_model

  # Update mean forest BART parameters
  mean_forest_params_default <- list(
    num_trees = 200,
    alpha = 0.95,
    beta = 2.0,
    min_samples_leaf = 5,
    max_depth = 10,
    sample_sigma2_leaf = TRUE,
    sigma2_leaf_init = NULL,
    sigma2_leaf_shape = 3,
    sigma2_leaf_scale = NULL,
    keep_vars = NULL,
    drop_vars = NULL,
    num_features_subsample = NULL,
    cloglog_leaf_prior_shape = 2.0,
    cloglog_leaf_prior_scale = 2.0
  )
  mean_forest_params_updated <- preprocessParams(
    mean_forest_params_default,
    mean_forest_params
  )

  # Update variance forest BART parameters
  variance_forest_params_default <- list(
    num_trees = 0,
    alpha = 0.95,
    beta = 2.0,
    min_samples_leaf = 5,
    max_depth = 10,
    leaf_prior_calibration_param = 1.5,
    var_forest_leaf_init = NULL,
    var_forest_prior_shape = NULL,
    var_forest_prior_scale = NULL,
    keep_vars = NULL,
    drop_vars = NULL,
    num_features_subsample = NULL
  )
  variance_forest_params_updated <- preprocessParams(
    variance_forest_params_default,
    variance_forest_params
  )

  # Update zero forest parameters (zero-inflated count outcomes only)
  zero_forest_params_default <- list(
    num_trees = 50,
    alpha = 0.95,
    beta = 2.0,
    min_samples_leaf = 5,
    max_depth = 10,
    sample_sigma2_leaf = FALSE,
    sigma2_leaf_init = NULL,
    sigma2_leaf_shape = 3,
    sigma2_leaf_scale = NULL,
    keep_vars = NULL,
    drop_vars = NULL,
    num_features_subsample = NULL
  )
  zero_forest_params_updated <- preprocessParams(
    zero_forest_params_default,
    zero_forest_params
  )

  # Update rfx parameters
  rfx_params_default <- list(
    model_spec = "custom",
    working_parameter_prior_mean = NULL,
    group_parameter_prior_mean = NULL,
    working_parameter_prior_cov = NULL,
    group_parameter_prior_cov = NULL,
    variance_prior_shape = 1,
    variance_prior_scale = 1,
    unseen_groups = "error"
  )
  rfx_params_updated <- preprocessParams(
    rfx_params_default,
    random_effects_params
  )

  ### Unpack all parameter values
  # 1. General parameters
  cutpoint_grid_size <- general_params_updated$cutpoint_grid_size
  standardize <- general_params_updated$standardize
  sample_sigma2_global <- general_params_updated$sample_sigma2_global
  sigma2_init <- general_params_updated$sigma2_global_init
  a_global <- general_params_updated$sigma2_global_shape
  b_global <- general_params_updated$sigma2_global_scale
  variable_weights <- general_params_updated$variable_weights
  random_seed <- general_params_updated$random_seed
  keep_burnin <- general_params_updated$keep_burnin
  keep_gfr <- general_params_updated$keep_gfr
  keep_every <- general_params_updated$keep_every
  num_chains <- general_params_updated$num_chains
  verbose <- general_params_updated$verbose
  outcome_model <- general_params_updated$outcome_model
  probit_outcome_model <- general_params_updated$probit_outcome_model
  num_threads <- general_params_updated$num_threads
  negbin_dispersion_prior_shape <- general_params_updated$negbin_dispersion_prior_shape
  negbin_dispersion_prior_rate <- general_params_updated$negbin_dispersion_prior_rate
  negbin_dispersion_init <- general_params_updated$negbin_dispersion_init
  sample_negbin_dispersion <- general_params_updated$sample_negbin_dispersion
  zero_random_effects <- general_params_updated$zero_random_effects

  # 2. Mean forest parameters
  num_trees_mean <- mean_forest_params_updated$num_trees
  alpha_mean <- mean_forest_params_updated$alpha
  beta_mean <- mean_forest_params_updated$beta
  min_samples_leaf_mean <- mean_forest_params_updated$min_samples_leaf
  max_depth_mean <- mean_forest_params_updated$max_depth
  sample_sigma2_leaf <- mean_forest_params_updated$sample_sigma2_leaf
  sigma2_leaf_init <- mean_forest_params_updated$sigma2_leaf_init
  a_leaf <- mean_forest_params_updated$sigma2_leaf_shape
  b_leaf <- mean_forest_params_updated$sigma2_leaf_scale
  keep_vars_mean <- mean_forest_params_updated$keep_vars
  drop_vars_mean <- mean_forest_params_updated$drop_vars
  num_features_subsample_mean <- mean_forest_params_updated$num_features_subsample
  cloglog_leaf_prior_shape <- mean_forest_params_updated$cloglog_leaf_prior_shape
  cloglog_leaf_prior_scale <- mean_forest_params_updated$cloglog_leaf_prior_scale

  # 3. Variance forest parameters
  num_trees_variance <- variance_forest_params_updated$num_trees
  alpha_variance <- variance_forest_params_updated$alpha
  beta_variance <- variance_forest_params_updated$beta
  min_samples_leaf_variance <- variance_forest_params_updated$min_samples_leaf
  max_depth_variance <- variance_forest_params_updated$max_depth
  a_0 <- variance_forest_params_updated$leaf_prior_calibration_param
  variance_forest_init <- variance_forest_params_updated$var_forest_leaf_init
  a_forest <- variance_forest_params_updated$var_forest_prior_shape
  b_forest <- variance_forest_params_updated$var_forest_prior_scale
  keep_vars_variance <- variance_forest_params_updated$keep_vars
  drop_vars_variance <- variance_forest_params_updated$drop_vars
  num_features_subsample_variance <- variance_forest_params_updated$num_features_subsample

  # 3b. Zero forest parameters
  num_trees_zero <- zero_forest_params_updated$num_trees
  alpha_zero <- zero_forest_params_updated$alpha
  beta_zero <- zero_forest_params_updated$beta
  min_samples_leaf_zero <- zero_forest_params_updated$min_samples_leaf
  max_depth_zero <- zero_forest_params_updated$max_depth
  sample_sigma2_leaf_zero <- zero_forest_params_updated$sample_sigma2_leaf
  sigma2_leaf_zero_init <- zero_forest_params_updated$sigma2_leaf_init
  a_leaf_zero <- zero_forest_params_updated$sigma2_leaf_shape
  b_leaf_zero <- zero_forest_params_updated$sigma2_leaf_scale
  keep_vars_zero <- zero_forest_params_updated$keep_vars
  drop_vars_zero <- zero_forest_params_updated$drop_vars
  num_features_subsample_zero <- zero_forest_params_updated$num_features_subsample

  # 4. RFX parameters
  rfx_model_spec <- rfx_params_updated$model_spec
  rfx_working_parameter_prior_mean <- rfx_params_updated$working_parameter_prior_mean
  rfx_group_parameter_prior_mean <- rfx_params_updated$group_parameter_prior_mean
  rfx_working_parameter_prior_cov <- rfx_params_updated$working_parameter_prior_cov
  rfx_group_parameter_prior_cov <- rfx_params_updated$group_parameter_prior_cov
  rfx_variance_prior_shape <- rfx_params_updated$variance_prior_shape
  rfx_variance_prior_scale <- rfx_params_updated$variance_prior_scale
  rfx_unseen_groups <- validateRfxUnseenGroups(
    rfx_params_updated$unseen_groups
  )
  .zinbCheckZeroRfxParams(zero_random_effects_params)

  # Validate before changing RNG state or allocating sampler objects. Random
  # effects are compatible with both retention modes: they are sampled as usual
  # and only the retained forests or posterior samples are kept.
  if (importance_only && (!identical(outcome_model$outcome, "continuous") ||
      !identical(outcome_model$link, "identity") || isTRUE(probit_outcome_model) ||
      num_chains != 1 || num_trees_mean <= 0 || num_trees_variance != 0 ||
      !is.null(X_test) || !is.null(leaf_basis_train) || !is.null(leaf_basis_test) ||
      !is.null(previous_model_json) || !is.null(previous_model_warmstart_sample_num))) {
    stop(paste("Importance retention requires a single-chain Gaussian constant-leaf",
      "mean forest without test data, variance forests or previous-model initialization."))
  }

  if (selected_only) {
    if (importance_only || !identical(outcome_model$outcome, "continuous") ||
        !identical(outcome_model$link, "identity") || isTRUE(probit_outcome_model) ||
        num_chains != 1 || num_trees_mean <= 0 || num_trees_variance != 0 ||
        isTRUE(keep_gfr) || isTRUE(keep_burnin) || num_mcmc < 1 ||
        !is.null(leaf_basis_train) || !is.null(leaf_basis_test) ||
        !is.null(previous_model_json) || !is.null(previous_model_warmstart_sample_num)) {
      stop("Selected retention requires a single-chain Gaussian constant-leaf mean forest, no retained warmup, and no variance forest or previous model.")
    }
    if (!is.numeric(retained_sample_ids) || !length(retained_sample_ids) ||
        anyNA(retained_sample_ids) || any(!is.finite(retained_sample_ids)) ||
        any(retained_sample_ids != floor(retained_sample_ids)) ||
        any(retained_sample_ids < 0 | retained_sample_ids >= num_mcmc) ||
        is.unsorted(retained_sample_ids, strictly = TRUE)) {
      stop("retained_sample_ids must be strictly increasing, unique, zero-based posterior sample IDs.")
    }
    retained_sample_ids <- as.integer(retained_sample_ids)
    selected_mask <- (seq_len(num_mcmc) - 1L) %in% retained_sample_ids
  }

  # Raise a deprecation warning to use `outcome_model` if `probit_outcome_model = TRUE` is specified
  if (probit_outcome_model) {
    warning(
      "Specifying a probit link through `general_params = list(probit_outcome_model = TRUE)` is deprecated and will be removed in a future version. Please use `general_params = list(outcome_model = OutcomeModel(outcome = 'binary', link = 'probit'))` instead."
    )
  }

  # Unpack outcome model details
  link_is_linear <- FALSE
  link_is_probit <- FALSE
  link_is_cloglog <- FALSE
  link_is_log <- FALSE
  outcome_is_continuous <- FALSE
  outcome_is_binary <- FALSE
  outcome_is_ordinal <- FALSE
  outcome_is_count <- FALSE
  if (
    outcome_model$outcome == "continuous" && outcome_model$link == "identity"
  ) {
    link_is_linear <- TRUE
    outcome_is_continuous <- TRUE
  } else if (
    outcome_model$outcome == "binary" && outcome_model$link == "probit"
  ) {
    link_is_probit <- TRUE
    outcome_is_binary <- TRUE
  } else if (
    outcome_model$outcome == "binary" && outcome_model$link == "cloglog"
  ) {
    link_is_cloglog <- TRUE
    outcome_is_binary <- TRUE
  } else if (
    outcome_model$outcome == "ordinal" && outcome_model$link == "cloglog"
  ) {
    link_is_cloglog <- TRUE
    outcome_is_ordinal <- TRUE
  } else if (outcome_model$outcome == "count" && outcome_model$link == "log") {
    link_is_log <- TRUE
    outcome_is_count <- TRUE
  } else {
    stop(paste0(
      "Invalid outcome model specification, outcome = ",
      outcome_model$outcome,
      ", link = ",
      outcome_model$link
    ))
  }

  # Zero-inflated negative binomial: a count outcome with family "zinb"
  outcome_is_zinb <- link_is_log && identical(outcome_model$family, "zinb")
  if (!outcome_is_zinb && length(zero_forest_params) > 0) {
    stop(
      "zero_forest_params is only used by zero-inflated count outcomes (OutcomeModel('count', 'log', family = 'zinb'))"
    )
  }
  if (
    !is.logical(zero_random_effects) ||
      length(zero_random_effects) != 1 ||
      is.na(zero_random_effects)
  ) {
    stop("general_params$zero_random_effects must be TRUE or FALSE")
  }

  # Count (negative-binomial) outcome: validate before changing RNG state
  if (!link_is_log) {
    if (!is.null(offset_train) || !is.null(offset_test)) {
      stop(
        "offset_train and offset_test are only supported for count outcomes (OutcomeModel(outcome = 'count', link = 'log'))"
      )
    }
  } else {
    if (num_trees_mean <= 0) {
      stop("A count outcome requires a mean forest (mean_forest_params$num_trees > 0)")
    }
    if (num_trees_variance > 0) {
      stop(
        "A variance forest cannot be combined with a count outcome: the negative-binomial variance is mu + mu^2 / r"
      )
    }
    if (isTRUE(general_params[["sample_sigma2_global"]])) {
      stop(
        "sample_sigma2_global = TRUE is not supported with a count outcome: the global error variance of the Polya-Gamma augmented model is fixed at 1"
      )
    }
    sample_sigma2_global <- FALSE
    y_count_check <- as.numeric(y_train)
    if (anyNA(y_count_check) || any(!is.finite(y_count_check))) {
      stop("A count outcome must not contain missing or infinite values")
    }
    if (any(y_count_check < 0) || any(y_count_check != round(y_count_check))) {
      stop("A count outcome must contain non-negative integers only")
    }
    if (all(y_count_check == 0)) {
      stop("A count outcome must contain at least one positive count")
    }
    if (outcome_is_zinb) {
      if (!any(y_count_check == 0)) {
        stop(
          "A zero-inflated count outcome needs at least one zero in y_train; use family = 'negbin'"
        )
      }
      if (!is.null(observation_weights)) {
        stop(
          "observation_weights are not supported with a zero-inflated count outcome"
        )
      }
      if (
        !is.numeric(num_trees_zero) ||
          length(num_trees_zero) != 1 ||
          num_trees_zero < 1
      ) {
        stop(
          "A zero-inflated count outcome requires a zero forest (zero_forest_params$num_trees >= 1)"
        )
      }
      if (
        !is.null(sigma2_leaf_zero_init) &&
          (!is.numeric(sigma2_leaf_zero_init) ||
            length(sigma2_leaf_zero_init) != 1 ||
            !is.finite(sigma2_leaf_zero_init) ||
            sigma2_leaf_zero_init <= 0)
      ) {
        stop("zero_forest_params$sigma2_leaf_init must be a single positive number")
      }
    }
    if (
      !is.numeric(negbin_dispersion_prior_shape) ||
        length(negbin_dispersion_prior_shape) != 1 ||
        !is.finite(negbin_dispersion_prior_shape) ||
        negbin_dispersion_prior_shape <= 0 ||
        !is.numeric(negbin_dispersion_prior_rate) ||
        length(negbin_dispersion_prior_rate) != 1 ||
        !is.finite(negbin_dispersion_prior_rate) ||
        negbin_dispersion_prior_rate < 0
    ) {
      stop(
        "negbin_dispersion_prior_shape must be a positive number and negbin_dispersion_prior_rate a non-negative number"
      )
    }
    if (
      !is.null(negbin_dispersion_init) &&
        (!is.numeric(negbin_dispersion_init) ||
          length(negbin_dispersion_init) != 1 ||
          !is.finite(negbin_dispersion_init) ||
          negbin_dispersion_init <= 0)
    ) {
      stop("negbin_dispersion_init must be a single positive number")
    }
    if (is.null(negbin_dispersion_init) && !isTRUE(sample_negbin_dispersion)) {
      stop("negbin_dispersion_init must be supplied when sample_negbin_dispersion = FALSE")
    }
    negbin_offset_train <- if (is.null(offset_train)) {
      rep(0, NROW(X_train))
    } else {
      .validateNegbinOffset(offset_train, NROW(X_train), "offset_train")
    }
    if (!is.null(offset_test)) {
      if (is.null(X_test)) {
        stop("offset_test was supplied without X_test")
      }
      negbin_offset_test <- .validateNegbinOffset(
        offset_test,
        NROW(X_test),
        "offset_test"
      )
    } else if (!is.null(X_test)) {
      if (!is.null(offset_train)) {
        warning(
          "offset_train was supplied but offset_test was not; test-set predictions use an offset of 0 (unit exposure)"
        )
      }
      negbin_offset_test <- rep(0, NROW(X_test))
    }
  }

  # Set a function-scoped RNG if user provided a random seed
  custom_rng <- random_seed >= 0
  has_existing_random_seed <- F
  if (custom_rng) {
    # Cache original global environment RNG state (if it exists)
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      original_global_seed <- .Random.seed
      has_existing_random_seed <- T
    }
    # Set new seed and store associated RNG state
    set.seed(random_seed)
  }

  # Check if there are enough GFR samples to seed num_chains samplers
  if (num_gfr > 0) {
    if (num_chains > num_gfr) {
      stop(
        "num_chains > num_gfr, meaning we do not have enough GFR samples to seed num_chains distinct MCMC chains"
      )
    }
  }

  # Override keep_gfr if there are no MCMC samples
  if (num_mcmc == 0) {
    keep_gfr <- TRUE
  }

  # Check if previous model JSON is provided and parse it if so
  has_prev_model <- !is.null(previous_model_json)
  has_prev_model_index <- !is.null(previous_model_warmstart_sample_num)
  if (has_prev_model) {
    previous_bart_model <- createBARTModelFromJsonString(
      previous_model_json
    )
    prev_num_samples <- previous_bart_model$model_params$num_samples
    if (!has_prev_model_index) {
      previous_model_warmstart_sample_num <- prev_num_samples
      warning(
        "`previous_model_warmstart_sample_num` was not provided alongside `previous_model_json`, so it will be set to the number of samples available in `previous_model_json`"
      )
    } else {
      if (previous_model_warmstart_sample_num < 1) {
        stop(
          "`previous_model_warmstart_sample_num` must be a positive integer"
        )
      }
      if (previous_model_warmstart_sample_num > prev_num_samples) {
        stop(
          "`previous_model_warmstart_sample_num` exceeds the number of samples in `previous_model_json`"
        )
      }
    }
    previous_model_decrement <- T
    if (num_chains > previous_model_warmstart_sample_num) {
      warning(
        "The number of chains being sampled exceeds the number of previous model samples available from the requested position in `previous_model_json`. All chains will be initialized from the same sample."
      )
      previous_model_decrement <- F
    }
    previous_y_bar <- previous_bart_model$model_params$outcome_mean
    previous_y_scale <- previous_bart_model$model_params$outcome_scale
    if (previous_bart_model$model_params$include_mean_forest) {
      previous_forest_samples_mean <- previous_bart_model$mean_forests
    } else {
      previous_forest_samples_mean <- NULL
    }
    if (previous_bart_model$model_params$include_variance_forest) {
      previous_forest_samples_variance <- previous_bart_model$variance_forests
    } else {
      previous_forest_samples_variance <- NULL
    }
    if (previous_bart_model$model_params$sample_sigma2_global) {
      previous_global_var_samples <- previous_bart_model$sigma2_global_samples /
        (previous_y_scale * previous_y_scale)
    } else {
      previous_global_var_samples <- NULL
    }
    if (previous_bart_model$model_params$sample_sigma2_leaf) {
      previous_leaf_var_samples <- previous_bart_model$sigma2_leaf_samples
    } else {
      previous_leaf_var_samples <- NULL
    }
    if (previous_bart_model$model_params$has_rfx) {
      previous_rfx_samples <- previous_bart_model$rfx_samples
    } else {
      previous_rfx_samples <- NULL
    }
    previous_is_count <- identical(
      previous_bart_model$model_params$outcome_model$outcome,
      "count"
    )
    if (previous_is_count != link_is_log) {
      stop(
        "previous_model_json and the current call must both use a count outcome, or neither"
      )
    }
    previous_negbin_dispersion_samples <- if (previous_is_count) {
      previous_bart_model$negbin_dispersion_samples
    } else {
      NULL
    }
    previous_is_zinb <- identical(
      previous_bart_model$model_params$outcome_model$family,
      "zinb"
    )
    if (previous_is_zinb != outcome_is_zinb) {
      stop(
        "previous_model_json and the current call must both use family = 'zinb', or neither"
      )
    }
    previous_forest_samples_zero <- if (previous_is_zinb) {
      previous_bart_model$zero_forests
    } else {
      NULL
    }
    if (
      previous_is_zinb &&
        previous_bart_model$model_params$num_trees_zero != num_trees_zero
    ) {
      stop(
        "zero_forest_params$num_trees (", num_trees_zero,
        ") must match the zero forest of previous_model_json (",
        previous_bart_model$model_params$num_trees_zero, " trees)"
      )
    }
    previous_leaf_var_zero_samples <- previous_bart_model$sigma2_leaf_zero_samples
    previous_zero_rfx_samples <- if (
      isTRUE(previous_bart_model$model_params$has_zero_rfx)
    ) {
      previous_bart_model$zero_rfx_samples
    } else {
      NULL
    }
    if (previous_bart_model$model_params$outcome_model$link == "cloglog") {
      previous_cloglog_cutpoint_samples <- previous_bart_model$cloglog_cutpoint_samples
      previous_cloglog_num_categories <- previous_bart_model$cloglog_num_categories
    } else {
      previous_cloglog_cutpoint_samples <- NULL
      previous_cloglog_num_categories <- 0
    }
    previous_model_num_samples <- previous_bart_model$model_params$num_samples
    if (previous_model_warmstart_sample_num > previous_model_num_samples) {
      stop(
        "`previous_model_warmstart_sample_num` exceeds the number of samples in `previous_model_json`"
      )
    }
  } else {
    previous_y_bar <- NULL
    previous_y_scale <- NULL
    previous_global_var_samples <- NULL
    previous_leaf_var_samples <- NULL
    previous_rfx_samples <- NULL
    previous_forest_samples_mean <- NULL
    previous_forest_samples_variance <- NULL
    previous_cloglog_cutpoint_samples <- NULL
    previous_cloglog_num_categories <- 0
    previous_model_num_samples <- 0
    previous_negbin_dispersion_samples <- NULL
    previous_forest_samples_zero <- NULL
    previous_leaf_var_zero_samples <- NULL
    previous_zero_rfx_samples <- NULL
  }

  # Determine whether conditional mean, variance, or both will be modeled
  if (num_trees_variance > 0) {
    include_variance_forest = TRUE
  } else {
    include_variance_forest = FALSE
  }
  if (num_trees_mean > 0) {
    include_mean_forest = TRUE
  } else {
    include_mean_forest = FALSE
  }

  # observation_weights compatibility checks
  if (!is.null(observation_weights)) {
    if (link_is_cloglog) {
      stop(
        "observation_weights are not compatible with cloglog link functions."
      )
    }
    if (include_variance_forest) {
      warning(
        "Results may be unreliable when observation_weights are deployed alongside a variance forest model."
      )
    }
  }

  # Set the variance forest priors if not set
  if (include_variance_forest) {
    if (is.null(a_forest)) {
      a_forest <- num_trees_variance / (a_0^2) + 0.5
    }
    if (is.null(b_forest)) b_forest <- num_trees_variance / (a_0^2)
  } else {
    a_forest <- 1.
    b_forest <- 1.
  }

  # Override tau sampling if there is no mean forest
  if (!include_mean_forest) {
    sample_sigma2_leaf <- FALSE
  }

  # Variable weight preprocessing (and initialization if necessary)
  if (is.null(variable_weights)) {
    variable_weights = rep(1 / ncol(X_train), ncol(X_train))
  }
  if (any(variable_weights < 0)) {
    stop("variable_weights cannot have any negative weights")
  }

  # Observation weight validation
  if (!is.null(observation_weights)) {
    if (!is.numeric(observation_weights)) {
      stop("observation_weights must be a numeric vector")
    }
    if (length(observation_weights) != nrow(X_train)) {
      stop("length(observation_weights) must equal nrow(X_train)")
    }
    if (any(observation_weights < 0)) {
      stop("observation_weights cannot have any negative values")
    }
    if (all(observation_weights == 0) && num_gfr > 0) {
      stop(
        "observation_weights are all zero (prior sampling mode) but num_gfr > 0. ",
        "GFR warm-start is data-dependent and ill-defined with zero weights. ",
        "Set num_gfr = 0 when using all-zero observation_weights."
      )
    }
  }

  # Check covariates are matrix or dataframe
  if ((!is.data.frame(X_train)) && (!is.matrix(X_train))) {
    stop("X_train must be a matrix or dataframe")
  }
  if (!is.null(X_test)) {
    if ((!is.data.frame(X_test)) && (!is.matrix(X_test))) {
      stop("X_test must be a matrix or dataframe")
    }
  }
  num_cov_orig <- ncol(X_train)

  # Raise a warning if the data have ties and only GFR is being run
  if ((num_gfr > 0) && (num_mcmc == 0) && (num_burnin == 0)) {
    num_values <- nrow(X_train)
    max_grid_size <- ifelse(
      num_values > cutpoint_grid_size,
      floor(num_values / cutpoint_grid_size),
      1
    )
    x_is_df <- is.data.frame(X_train)
    covs_warning_1 <- NULL
    covs_warning_2 <- NULL
    covs_warning_3 <- NULL
    covs_warning_4 <- NULL
    for (i in 1:num_cov_orig) {
      # Skip check for variables that are treated as categorical
      x_numeric <- T
      if (x_is_df) {
        if (is.factor(X_train[, i])) {
          x_numeric <- F
        }
      }
      if (x_numeric) {
        # Determine the number of unique values
        num_unique_values <- length(unique(X_train[, i]))

        # Determine a "name" for the covariate
        cov_name <- ifelse(
          is.null(colnames(X_train)),
          paste0("X", i),
          colnames(X_train)[i]
        )

        # Check for a small relative number of unique values
        unique_full_ratio <- num_unique_values / num_values
        if (unique_full_ratio < 0.2) {
          covs_warning_1 <- c(covs_warning_1, cov_name)
        }

        # Check for a small absolute number of unique values
        if (num_values > 100) {
          if (num_unique_values < 20) {
            covs_warning_2 <- c(covs_warning_2, cov_name)
          }
        }

        # Check for a large number of duplicates of any individual value
        x_j_hist <- table(X_train[, i])
        if (any(x_j_hist > 2 * max_grid_size)) {
          covs_warning_3 <- c(covs_warning_3, cov_name)
        }

        # Check for binary variables
        if (num_unique_values == 2) {
          covs_warning_4 <- c(covs_warning_4, cov_name)
        }
      }
    }

    if (!is.null(covs_warning_1)) {
      warning(
        paste0(
          "Covariate(s) ",
          paste(covs_warning_1, collapse = ", "),
          " have a ratio of unique to overall observations of less than 0.2. ",
          "This might present some issues with the grow-from-root (GFR) algorithm. ",
          "Consider running with `num_mcmc > 0` and `num_burnin > 0` to improve your model's performance."
        )
      )
    }

    if (!is.null(covs_warning_2)) {
      warning(
        paste0(
          "Covariate(s) ",
          paste(covs_warning_2, collapse = ", "),
          " have fewer than 20 unique values. ",
          "This might present some issues with the grow-from-root (GFR) algorithm. ",
          "Consider running with `num_mcmc > 0` and `num_burnin > 0` to improve your model's performance."
        )
      )
    }

    if (!is.null(covs_warning_3)) {
      warning(
        paste0(
          "Covariates ",
          paste(covs_warning_3, collapse = ", "),
          " have some observed values with more than ",
          2 * max_grid_size,
          " repeated observations. ",
          "This might present some issues with the grow-from-root (GFR) algorithm. ",
          "Consider running with `num_mcmc > 0` and `num_burnin > 0` to improve your model's performance."
        )
      )
    }

    if (!is.null(covs_warning_4)) {
      warning(
        paste0(
          "Covariates ",
          paste(covs_warning_4, collapse = ", "),
          " appear to be binary but are currently treated by stochtree as continuous. ",
          "This might present some issues with the grow-from-root (GFR) algorithm. ",
          "Consider converting binary variables to ordered factor (i.e. `factor(..., ordered = T)`."
        )
      )
    }
  }

  # Standardize the keep variable lists to numeric indices
  if (!is.null(keep_vars_mean)) {
    if (is.character(keep_vars_mean)) {
      if (!all(keep_vars_mean %in% names(X_train))) {
        stop(
          "keep_vars_mean includes some variable names that are not in X_train"
        )
      }
      variable_subset_mean <- unname(which(
        names(X_train) %in% keep_vars_mean
      ))
    } else {
      if (any(keep_vars_mean > ncol(X_train))) {
        stop(
          "keep_vars_mean includes some variable indices that exceed the number of columns in X_train"
        )
      }
      if (any(keep_vars_mean < 0)) {
        stop("keep_vars_mean includes some negative variable indices")
      }
      variable_subset_mean <- keep_vars_mean
    }
  } else if ((is.null(keep_vars_mean)) && (!is.null(drop_vars_mean))) {
    if (is.character(drop_vars_mean)) {
      if (!all(drop_vars_mean %in% names(X_train))) {
        stop(
          "drop_vars_mean includes some variable names that are not in X_train"
        )
      }
      variable_subset_mean <- unname(which(
        !(names(X_train) %in% drop_vars_mean)
      ))
    } else {
      if (any(drop_vars_mean > ncol(X_train))) {
        stop(
          "drop_vars_mean includes some variable indices that exceed the number of columns in X_train"
        )
      }
      if (any(drop_vars_mean < 0)) {
        stop("drop_vars_mean includes some negative variable indices")
      }
      variable_subset_mean <- (1:ncol(X_train))[
        !(1:ncol(X_train) %in% drop_vars_mean)
      ]
    }
  } else {
    variable_subset_mean <- 1:ncol(X_train)
  }
  if (!is.null(keep_vars_variance)) {
    if (is.character(keep_vars_variance)) {
      if (!all(keep_vars_variance %in% names(X_train))) {
        stop(
          "keep_vars_variance includes some variable names that are not in X_train"
        )
      }
      variable_subset_variance <- unname(which(
        names(X_train) %in% keep_vars_variance
      ))
    } else {
      if (any(keep_vars_variance > ncol(X_train))) {
        stop(
          "keep_vars_variance includes some variable indices that exceed the number of columns in X_train"
        )
      }
      if (any(keep_vars_variance < 0)) {
        stop(
          "keep_vars_variance includes some negative variable indices"
        )
      }
      variable_subset_variance <- keep_vars_variance
    }
  } else if ((is.null(keep_vars_variance)) && (!is.null(drop_vars_variance))) {
    if (is.character(drop_vars_variance)) {
      if (!all(drop_vars_variance %in% names(X_train))) {
        stop(
          "drop_vars_variance includes some variable names that are not in X_train"
        )
      }
      variable_subset_variance <- unname(which(
        !(names(X_train) %in% drop_vars_variance)
      ))
    } else {
      if (any(drop_vars_variance > ncol(X_train))) {
        stop(
          "drop_vars_variance includes some variable indices that exceed the number of columns in X_train"
        )
      }
      if (any(drop_vars_variance < 0)) {
        stop(
          "drop_vars_variance includes some negative variable indices"
        )
      }
      variable_subset_variance <- (1:ncol(X_train))[
        !(1:ncol(X_train) %in% drop_vars_variance)
      ]
    }
  } else {
    variable_subset_variance <- 1:ncol(X_train)
  }
  if (outcome_is_zinb) {
    variable_subset_zero <- .resolveForestVariableSubset(
      X_train,
      keep_vars_zero,
      drop_vars_zero,
      "zero_forest_params$keep_vars",
      "zero_forest_params$drop_vars"
    )
  }

  # Preprocess covariates
  if ((!is.data.frame(X_train)) && (!is.matrix(X_train))) {
    stop("X_train must be a matrix or dataframe")
  }
  if (!is.null(X_test)) {
    if ((!is.data.frame(X_test)) && (!is.matrix(X_test))) {
      stop("X_test must be a matrix or dataframe")
    }
  }
  if (ncol(X_train) != length(variable_weights)) {
    stop("length(variable_weights) must equal ncol(X_train)")
  }
  train_cov_preprocess_list <- preprocessTrainData(X_train)
  X_train_metadata <- train_cov_preprocess_list$metadata
  X_train <- train_cov_preprocess_list$data
  original_var_indices <- X_train_metadata$original_var_indices
  feature_types <- X_train_metadata$feature_types
  if (!is.null(X_test)) {
    X_test <- preprocessPredictionData(X_test, X_train_metadata)
  }

  # Update variable weights
  variable_weights_mean <- variable_weights_variance <- variable_weights
  variable_weights_adj <- 1 /
    sapply(original_var_indices, function(x) sum(original_var_indices == x))
  if (include_mean_forest) {
    variable_weights_mean <- variable_weights_mean[original_var_indices] *
      variable_weights_adj
    variable_weights_mean[
      !(original_var_indices %in% variable_subset_mean)
    ] <- 0
  }
  if (include_variance_forest) {
    variable_weights_variance <- variable_weights_variance[
      original_var_indices
    ] *
      variable_weights_adj
    variable_weights_variance[
      !(original_var_indices %in% variable_subset_variance)
    ] <- 0
  }
  if (outcome_is_zinb) {
    variable_weights_zero <- variable_weights[original_var_indices] *
      variable_weights_adj
    variable_weights_zero[
      !(original_var_indices %in% variable_subset_zero)
    ] <- 0
  }

  # Set num_features_subsample to default, ncol(X_train), if not already set
  if (is.null(num_features_subsample_mean)) {
    num_features_subsample_mean <- ncol(X_train)
  }
  if (is.null(num_features_subsample_variance)) {
    num_features_subsample_variance <- ncol(X_train)
  }
  if (is.null(num_features_subsample_zero)) {
    num_features_subsample_zero <- ncol(X_train)
  }

  # Convert all input data to matrices if not already converted
  if ((is.null(dim(leaf_basis_train))) && (!is.null(leaf_basis_train))) {
    leaf_basis_train <- as.matrix(leaf_basis_train)
  }
  if ((is.null(dim(leaf_basis_test))) && (!is.null(leaf_basis_test))) {
    leaf_basis_test <- as.matrix(leaf_basis_test)
  }
  if ((is.null(dim(rfx_basis_train))) && (!is.null(rfx_basis_train))) {
    rfx_basis_train <- as.matrix(rfx_basis_train)
  }
  if ((is.null(dim(rfx_basis_test))) && (!is.null(rfx_basis_test))) {
    rfx_basis_test <- as.matrix(rfx_basis_test)
  }

  # Recode group IDs to integer vector (if passed as, for example, a vector of county names, etc...)
  has_rfx <- FALSE
  has_rfx_test <- FALSE
  if (!is.null(rfx_group_ids_train)) {
    group_ids_factor <- factor(rfx_group_ids_train)
    rfx_group_ids_train <- as.integer(group_ids_factor)
    has_rfx <- TRUE
    if (!is.null(rfx_group_ids_test)) {
      group_ids_factor_test <- factor(
        rfx_group_ids_test,
        levels = levels(group_ids_factor)
      )
      if (
        sum(is.na(group_ids_factor_test)) > 0 &&
          identical(rfx_unseen_groups, "error")
      ) {
        stop(
          "All random effect group labels provided in rfx_group_ids_test must be present in rfx_group_ids_train"
        )
      }
      rfx_group_ids_test <- as.integer(group_ids_factor_test)
      has_rfx_test <- TRUE
    }
  }

  # Data consistency checks
  if ((!is.null(X_test)) && (ncol(X_test) != ncol(X_train))) {
    stop("X_train and X_test must have the same number of columns")
  }
  if (
    (!is.null(leaf_basis_test)) &&
      (ncol(leaf_basis_test) != ncol(leaf_basis_train))
  ) {
    stop(
      "leaf_basis_train and leaf_basis_test must have the same number of columns"
    )
  }
  if (
    (!is.null(leaf_basis_train)) &&
      (nrow(leaf_basis_train) != nrow(X_train))
  ) {
    stop("leaf_basis_train and X_train must have the same number of rows")
  }
  if ((!is.null(leaf_basis_test)) && (nrow(leaf_basis_test) != nrow(X_test))) {
    stop("leaf_basis_test and X_test must have the same number of rows")
  }
  if (nrow(X_train) != length(y_train)) {
    stop("X_train and y_train must have the same number of observations")
  }
  if (
    (!is.null(rfx_basis_test)) &&
      (ncol(rfx_basis_test) != ncol(rfx_basis_train))
  ) {
    stop(
      "rfx_basis_train and rfx_basis_test must have the same number of columns"
    )
  }
  if (!is.null(rfx_group_ids_train)) {
    if (!is.null(rfx_group_ids_test)) {
      if ((!is.null(rfx_basis_train)) && (is.null(rfx_basis_test))) {
        stop(
          "rfx_basis_train is provided but rfx_basis_test is not provided"
        )
      }
    }
  }

  # Handle the rfx basis matrices
  has_basis_rfx <- FALSE
  num_basis_rfx <- 0
  if (has_rfx) {
    if (rfx_model_spec == "custom") {
      if (is.null(rfx_basis_train)) {
        stop(
          "A user-provided basis (`rfx_basis_train`) must be provided when the random effects model spec is 'custom'"
        )
      }
      has_basis_rfx <- TRUE
      num_basis_rfx <- ncol(rfx_basis_train)
    } else if (rfx_model_spec == "intercept_only") {
      rfx_basis_train <- matrix(
        rep(1, nrow(X_train)),
        nrow = nrow(X_train),
        ncol = 1
      )
      has_basis_rfx <- TRUE
      num_basis_rfx <- 1
    }
    num_rfx_groups <- length(unique(rfx_group_ids_train))
    num_rfx_components <- ncol(rfx_basis_train)
    if (num_rfx_groups == 1) {
      warning(
        "Only one group was provided for random effect sampling, so the random effects model is likely overkill"
      )
    }
  }
  # Zero-inflated family: the zero part gets its own random effects on the
  # same groups and basis unless general_params$zero_random_effects = FALSE
  has_zero_rfx <- outcome_is_zinb && has_rfx && zero_random_effects
  if (length(zero_random_effects_params) > 0 && !has_zero_rfx) {
    stop(
      "zero_random_effects_params applies only to family = 'zinb' with rfx_group_ids_train and general_params$zero_random_effects = TRUE"
    )
  }
  if (has_rfx_test) {
    if (rfx_model_spec == "custom") {
      if (is.null(rfx_basis_test)) {
        stop(
          "A user-provided basis (`rfx_basis_test`) must be provided when the random effects model spec is 'custom'"
        )
      }
    } else if (rfx_model_spec == "intercept_only") {
      rfx_basis_test <- matrix(
        rep(1, nrow(X_test)),
        nrow = nrow(X_test),
        ncol = 1
      )
    }
  }

  # Convert y_train to numeric vector if not already converted
  if (!is.null(dim(y_train))) {
    y_train <- as.matrix(y_train)
  }

  # Determine whether a basis vector is provided
  has_basis = !is.null(leaf_basis_train)

  # Determine whether a test set is provided
  has_test = !is.null(X_test)

  # Preliminary runtime checks for probit link
  if (!include_mean_forest) {
    link_is_probit <- FALSE
    # TODO: think about allowing binary models with probit link for homoskedastic RFX-only models?
  }
  if (link_is_probit) {
    if (!(length(unique(y_train)) == 2)) {
      stop(
        "You specified a probit link, but supplied an outcome with more than 2 unique values. Probit is only currently supported for binary outcomes."
      )
    }
    unique_outcomes <- sort(unique(y_train))
    if (!(all(unique_outcomes == c(0, 1)))) {
      stop(
        "You specified a probit link, but supplied an outcome with 2 unique values other than 0 and 1"
      )
    }
    if (include_variance_forest) {
      stop("We do not support heteroskedasticity with a probit link")
    }
    if (sample_sigma2_global) {
      warning(
        "Global error variance will not be sampled with a probit link as it is fixed at 1"
      )
      sample_sigma2_global <- F
    }
  }

  # Preliminary runtime checks for cloglog link
  if (!include_mean_forest) {
    link_is_cloglog <- FALSE
    # TODO: think about allowing binary models with cloglog link for homoskedastic RFX-only models?
  }
  if (link_is_cloglog) {
    if (!all(as.integer(y_train) == y_train)) {
      stop(
        "You specified a cloglog link, but supplied an outcome with non-integer values. Cloglog is only currently supported for integer outcomes."
      )
    }
    unique_outcomes <- sort(unique(y_train))
    if (!(min(unique_outcomes) %in% c(0, 1))) {
      stop(
        "You specified a cloglog link, but supplied an integer outcome that does not start with 0 or 1. Please remap / shift the outcomes so that the smallest category label is either 0 or 1."
      )
    }
    if (!all(diff(unique_outcomes) == 1)) {
      stop(
        "You specified a cloglog link, but supplied an integer outcome that is not a sequence of consecutive integers"
      )
    }
    if (include_variance_forest) {
      stop("We do not support heteroskedasticity with a cloglog link")
    }
    if (has_basis) {
      stop("We do not support leaf basis regression with a cloglog link")
    }
    if (sample_sigma2_global) {
      warning(
        "Global error variance will not be sampled with a cloglog link"
      )
      sample_sigma2_global <- F
    }
    if (sample_sigma2_leaf) {
      warning(
        "Leaf scale parameter will not be sampled with a cloglog link"
      )
      sample_sigma2_leaf <- F
    }
  }

  # Runtime checks for variance forest
  if (include_variance_forest) {
    if (sample_sigma2_global) {
      warning(
        "Global error variance will not be sampled with a heteroskedasticity forest"
      )
      sample_sigma2_global <- F
    }
  }

  # Handle standardization, prior calibration, and initialization of forest
  # differently for binary and continuous outcomes
  if (link_is_probit) {
    # Probit-scale intercept: center the forest on the population-average latent mean.
    # The forest predicts mu(X) and y_bar_train is added back at prediction time.
    # The latent z sampling uses y_bar_train to set the correct truncated normal mean and to center z before the residual update.
    y_bar_train <- qnorm(mean_cpp(as.numeric(y_train)))
    y_std_train <- 1
    standardize <- FALSE

    # Set a pseudo outcome by subtracting mean_cpp(y_train) from y_train
    resid_train <- y_train - mean_cpp(as.numeric(y_train))

    # Set initial values of root nodes to 0.0 (in probit scale)
    init_val_mean <- 0.0

    # Calibrate priors for sigma^2 and tau
    # Set sigma2_init to 1, ignoring default provided
    sigma2_init <- 1.0
    # Skip variance_forest_init, since variance forests are not supported with probit link
    if (is.null(b_leaf)) {
      b_leaf <- 1 / (num_trees_mean)
    }
    if (has_basis) {
      if (ncol(leaf_basis_train) > 1) {
        if (is.null(sigma2_leaf_init)) {
          sigma2_leaf_init <- diag(
            2 / (num_trees_mean),
            ncol(leaf_basis_train)
          )
        }
        if (!is.matrix(sigma2_leaf_init)) {
          current_leaf_scale <- as.matrix(diag(
            sigma2_leaf_init,
            ncol(leaf_basis_train)
          ))
        } else {
          current_leaf_scale <- sigma2_leaf_init
        }
      } else {
        if (is.null(sigma2_leaf_init)) {
          sigma2_leaf_init <- as.matrix(2 / (num_trees_mean))
        }
        if (!is.matrix(sigma2_leaf_init)) {
          current_leaf_scale <- as.matrix(diag(sigma2_leaf_init, 1))
        } else {
          current_leaf_scale <- sigma2_leaf_init
        }
      }
    } else {
      if (is.null(sigma2_leaf_init)) {
        sigma2_leaf_init <- as.matrix(2 / (num_trees_mean))
      }
      if (!is.matrix(sigma2_leaf_init)) {
        current_leaf_scale <- as.matrix(diag(sigma2_leaf_init, 1))
      } else {
        current_leaf_scale <- sigma2_leaf_init
      }
    }
    current_sigma2 <- sigma2_init
  } else if (link_is_cloglog) {
    # Fix offset to 0 and scale to 1
    y_bar_train <- 0
    y_std_train <- 1
    standardize <- FALSE

    # Remap outcomes to start from 0
    resid_train <- as.numeric(y_train - min(unique_outcomes))
    cloglog_num_categories <- max(resid_train) + 1

    # Set initial values of root nodes to 0.0 (in linear scale)
    init_val_mean <- 0.0

    # Calibrate priors for sigma^2 and tau
    # Set sigma2_init to 1, ignoring default provided
    sigma2_init <- 1.0
    if (is.null(sigma2_leaf_init)) {
      sigma2_leaf_init <- as.matrix(2 / (num_trees_mean))
    }
    current_sigma2 <- sigma2_init
    current_leaf_scale <- sigma2_leaf_init

    # Set first cutpoint to 0 for identifiability
    cloglog_cutpoint_0 <- 0

    # Set shape and rate parameters for conditional gamma model
    cloglog_forest_shape <- 2.0
    cloglog_forest_rate <- 2.0
  } else if (link_is_log) {
    # Negative-binomial outcome on the log scale. The centering constant c plays
    # the role of y_bar_train for probit: the forest models log(mu) - offset - c
    # and c is added back at prediction time. c is the maximum-likelihood
    # constant log rate, log(sum(y) / sum(exp(offset))).
    negbin_y <- as.numeric(y_train)
    y_bar_train <- log(sum(negbin_y) / sum(exp(negbin_offset_train)))
    y_std_train <- 1
    standardize <- FALSE

    # Zero-inflated family: centre the count forest on the log rate of the
    # constant-rate ZINB fit and the zero forest on qnorm of its zero
    # probability. The first structural-zero draw is made from the root
    # forests, so it is the posterior of s_i under that constant fit: this is
    # how s is initialised from the excess-zero fraction.
    if (outcome_is_zinb) {
      zinb_init <- .zinbConstantFit(negbin_y, negbin_offset_train)
      y_bar_train <- zinb_init$log_rate
      zinb_zero_probability_init <- zinb_init$zero_probability
      zero_outcome_mean <- qnorm(min(max(zinb_zero_probability_init, 0.01), 0.99))
      # Probit calibration of the zero forest, as for a probit mean forest
      if (is.null(sigma2_leaf_zero_init)) {
        sigma2_leaf_zero_init <- 2 / num_trees_zero
      }
      if (is.null(b_leaf_zero)) {
        b_leaf_zero <- 1 / num_trees_zero
      }
      current_leaf_scale_zero <- as.matrix(sigma2_leaf_zero_init)
      current_structural_zero <- integer(length(negbin_y))
    }

    # Placeholder residual; the first Polya-Gamma step overwrites it
    resid_train <- rep(0, length(negbin_y))
    init_val_mean <- 0.0

    # The augmented Gaussian regression has its error variance fixed at 1
    sigma2_init <- 1.0
    current_sigma2 <- sigma2_init

    # Starting value of the dispersion r
    if (is.null(negbin_dispersion_init)) {
      negbin_dispersion_init <- if (outcome_is_zinb) {
        zinb_init$dispersion
      } else {
        .negbinDispersionInit(
          negbin_y,
          exp(negbin_offset_train + y_bar_train)
        )
      }
    }
    current_dispersion <- negbin_dispersion_init

    # Leaf scale calibrated on the log scale (Chipman, George and McCulloch
    # 2010, with k = 2): leaf sd = R / (2 k sqrt(num_trees)), R the range of log1p(y)
    log_count_range <- diff(range(log1p(negbin_y)))
    if (!is.finite(log_count_range) || log_count_range <= 0) {
      log_count_range <- 1
    }
    negbin_leaf_var <- (log_count_range / (2 * 2 * sqrt(num_trees_mean)))^2
    if (is.null(b_leaf)) {
      b_leaf <- negbin_leaf_var / 2
    }
    leaf_dim_init <- if (has_basis) ncol(leaf_basis_train) else 1
    if (is.null(sigma2_leaf_init)) {
      sigma2_leaf_init <- if (leaf_dim_init > 1) {
        diag(negbin_leaf_var, leaf_dim_init)
      } else {
        as.matrix(negbin_leaf_var)
      }
    }
    if (!is.matrix(sigma2_leaf_init)) {
      current_leaf_scale <- as.matrix(diag(sigma2_leaf_init, leaf_dim_init))
    } else {
      current_leaf_scale <- sigma2_leaf_init
    }

    # Variance weights: user observation weights (if any) times 1 / omega
    negbin_base_variance_weights <- if (is.null(observation_weights)) {
      rep(1, length(negbin_y))
    } else {
      as.numeric(observation_weights)
    }
  } else {
    # Only standardize if user requested
    if (standardize) {
      y_bar_train <- mean_cpp(as.numeric(y_train))
      y_std_train <- sd_cpp(as.numeric(y_train))
    } else {
      y_bar_train <- 0
      y_std_train <- 1
    }

    # Compute standardized outcome
    resid_train <- (y_train - y_bar_train) / y_std_train

    # Compute initial value of root nodes in mean forest
    init_val_mean <- mean_cpp(as.numeric(resid_train))

    # Calibrate priors for sigma^2 and tau
    if (is.null(sigma2_init)) {
      sigma2_init <- 1.0 * var_cpp(as.numeric(resid_train))
    }
    if (is.null(variance_forest_init)) {
      variance_forest_init <- 1.0 * var_cpp(as.numeric(resid_train))
    }
    if (is.null(b_leaf)) {
      b_leaf <- var_cpp(as.numeric(resid_train)) / (2 * num_trees_mean)
    }
    if (has_basis) {
      if (ncol(leaf_basis_train) > 1) {
        if (is.null(sigma2_leaf_init)) {
          sigma2_leaf_init <- diag(
            2 * var_cpp(as.numeric(resid_train)) / (num_trees_mean),
            ncol(leaf_basis_train)
          )
        }
        if (!is.matrix(sigma2_leaf_init)) {
          current_leaf_scale <- as.matrix(diag(
            sigma2_leaf_init,
            ncol(leaf_basis_train)
          ))
        } else {
          current_leaf_scale <- sigma2_leaf_init
        }
      } else {
        if (is.null(sigma2_leaf_init)) {
          sigma2_leaf_init <- as.matrix(
            2 * var_cpp(as.numeric(resid_train)) / (num_trees_mean)
          )
        }
        if (!is.matrix(sigma2_leaf_init)) {
          current_leaf_scale <- as.matrix(diag(sigma2_leaf_init, 1))
        } else {
          current_leaf_scale <- sigma2_leaf_init
        }
      }
    } else {
      if (is.null(sigma2_leaf_init)) {
        sigma2_leaf_init <- as.matrix(
          2 * var_cpp(as.numeric(resid_train)) / (num_trees_mean)
        )
      }
      if (!is.matrix(sigma2_leaf_init)) {
        current_leaf_scale <- as.matrix(diag(sigma2_leaf_init, 1))
      } else {
        current_leaf_scale <- sigma2_leaf_init
      }
    }
    current_sigma2 <- sigma2_init
  }

  # Determine leaf model type
  if ((!has_basis) && (!link_is_cloglog)) {
    leaf_model_mean_forest <- 0
  } else if ((!has_basis) && (link_is_cloglog)) {
    leaf_model_mean_forest <- 4
  } else if (ncol(leaf_basis_train) == 1) {
    leaf_model_mean_forest <- 1
  } else if (ncol(leaf_basis_train) > 1) {
    leaf_model_mean_forest <- 2
  } else {
    stop("leaf_basis_train passed must be a matrix with at least 1 column")
  }

  # Set variance leaf model type (currently only one option)
  leaf_model_variance_forest <- 3

  # Unpack model type info
  if (leaf_model_mean_forest == 0) {
    leaf_dimension = 1
    is_leaf_constant = TRUE
    leaf_regression = FALSE
  } else if (leaf_model_mean_forest == 1) {
    stopifnot(has_basis)
    stopifnot(ncol(leaf_basis_train) == 1)
    leaf_dimension = 1
    is_leaf_constant = FALSE
    leaf_regression = TRUE
  } else if (leaf_model_mean_forest == 2) {
    stopifnot(has_basis)
    stopifnot(ncol(leaf_basis_train) > 1)
    leaf_dimension = ncol(leaf_basis_train)
    is_leaf_constant = FALSE
    leaf_regression = TRUE
    if (sample_sigma2_leaf) {
      warning(
        "Sampling leaf scale not yet supported for multivariate leaf models, so the leaf scale parameter will not be sampled in this model."
      )
      sample_sigma2_leaf <- FALSE
    }
  } else if (leaf_model_mean_forest == 4) {
    leaf_dimension = 1
    is_leaf_constant = TRUE
    leaf_regression = FALSE
  }

  # Data. A count outcome always carries variance weights (1 / omega), which
  # the Polya-Gamma step overwrites every iteration.
  forest_variance_weights <- if (link_is_log) {
    negbin_base_variance_weights
  } else {
    observation_weights
  }
  if (leaf_regression) {
    forest_dataset_train <- createForestDataset(
      X_train,
      leaf_basis_train,
      forest_variance_weights
    )
    if (has_test) {
      forest_dataset_test <- createForestDataset(X_test, leaf_basis_test)
    }
    requires_basis <- TRUE
  } else {
    forest_dataset_train <- createForestDataset(
      X_train,
      variance_weights = forest_variance_weights
    )
    if (has_test) {
      forest_dataset_test <- createForestDataset(X_test)
    }
    requires_basis <- FALSE
  }
  outcome_train <- createOutcome(resid_train)
  # The zero forest of a zero-inflated count outcome sees the covariates
  # without the Polya-Gamma variance weights, and its own latent outcome
  if (outcome_is_zinb) {
    forest_dataset_zero <- createForestDataset(X_train)
    outcome_zero <- createOutcome(rep(0, nrow(X_train)))
  }

  # Random number generator (std::mt19937)
  if (is.null(random_seed)) {
    random_seed = sample(1:10000, 1, FALSE)
  }
  rng <- createCppRNG(random_seed)

  # Separate ordinal sampler object for cloglog
  if (link_is_cloglog) {
    ordinal_sampler <- ordinal_sampler_cpp()
  }

  # Sampling data structures
  feature_types <- as.integer(feature_types)
  global_model_config <- createGlobalModelConfig(
    global_error_variance = current_sigma2
  )
  if (include_mean_forest) {
    forest_model_config_mean <- createForestModelConfig(
      feature_types = feature_types,
      num_trees = num_trees_mean,
      num_features = ncol(X_train),
      num_observations = nrow(X_train),
      variable_weights = variable_weights_mean,
      leaf_dimension = leaf_dimension,
      alpha = alpha_mean,
      beta = beta_mean,
      min_samples_leaf = min_samples_leaf_mean,
      max_depth = max_depth_mean,
      leaf_model_type = leaf_model_mean_forest,
      leaf_model_scale = current_leaf_scale,
      cutpoint_grid_size = cutpoint_grid_size,
      num_features_subsample = num_features_subsample_mean
    )
    if (link_is_cloglog) {
      forest_model_config_mean$update_cloglog_forest_shape(cloglog_forest_shape)
      forest_model_config_mean$update_cloglog_forest_rate(cloglog_forest_rate)
    }
    forest_model_mean <- createForestModel(
      forest_dataset_train,
      forest_model_config_mean,
      global_model_config
    )
  }
  if (include_variance_forest) {
    forest_model_config_variance <- createForestModelConfig(
      feature_types = feature_types,
      num_trees = num_trees_variance,
      num_features = ncol(X_train),
      num_observations = nrow(X_train),
      variable_weights = variable_weights_variance,
      leaf_dimension = 1,
      alpha = alpha_variance,
      beta = beta_variance,
      min_samples_leaf = min_samples_leaf_variance,
      max_depth = max_depth_variance,
      leaf_model_type = leaf_model_variance_forest,
      variance_forest_shape = a_forest,
      variance_forest_scale = b_forest,
      cutpoint_grid_size = cutpoint_grid_size,
      num_features_subsample = num_features_subsample_variance
    )
    forest_model_variance <- createForestModel(
      forest_dataset_train,
      forest_model_config_variance,
      global_model_config
    )
  }
  if (outcome_is_zinb) {
    # Probit zero forest: latent error variance fixed at 1
    global_model_config_zero <- createGlobalModelConfig(
      global_error_variance = 1.0
    )
    forest_model_config_zero <- createForestModelConfig(
      feature_types = feature_types,
      num_trees = num_trees_zero,
      num_features = ncol(X_train),
      num_observations = nrow(X_train),
      variable_weights = variable_weights_zero,
      leaf_dimension = 1,
      alpha = alpha_zero,
      beta = beta_zero,
      min_samples_leaf = min_samples_leaf_zero,
      max_depth = max_depth_zero,
      leaf_model_type = 0,
      leaf_model_scale = current_leaf_scale_zero,
      cutpoint_grid_size = cutpoint_grid_size,
      num_features_subsample = num_features_subsample_zero
    )
    forest_model_zero <- createForestModel(
      forest_dataset_zero,
      forest_model_config_zero,
      global_model_config_zero
    )
    forest_samples_zero <- createForestSamples(num_trees_zero, 1, TRUE, FALSE)
    active_forest_zero <- createForest(num_trees_zero, 1, TRUE, FALSE)
  }

  # Container of forest samples
  if (include_mean_forest) {
    forest_samples_mean <- createForestSamples(
      num_trees_mean,
      leaf_dimension,
      is_leaf_constant,
      FALSE
    )
    active_forest_mean <- createForest(
      num_trees_mean,
      leaf_dimension,
      is_leaf_constant,
      FALSE
    )
  }
  if (include_variance_forest) {
    forest_samples_variance <- createForestSamples(
      num_trees_variance,
      1,
      TRUE,
      TRUE
    )
    active_forest_variance <- createForest(
      num_trees_variance,
      1,
      TRUE,
      TRUE
    )
  }

  # Random effects initialization
  if (has_rfx) {
    # Prior parameters
    if (is.null(rfx_working_parameter_prior_mean)) {
      if (num_rfx_components == 1) {
        alpha_init <- c(0)
      } else if (num_rfx_components > 1) {
        alpha_init <- rep(0, num_rfx_components)
      } else {
        stop("There must be at least 1 random effect component")
      }
    } else {
      alpha_init <- expand_dims_1d(
        rfx_working_parameter_prior_mean,
        num_rfx_components
      )
    }

    if (is.null(rfx_group_parameter_prior_mean)) {
      xi_init <- matrix(
        rep(alpha_init, num_rfx_groups),
        num_rfx_components,
        num_rfx_groups
      )
    } else {
      xi_init <- expand_dims_2d(
        rfx_group_parameter_prior_mean,
        num_rfx_components,
        num_rfx_groups
      )
    }

    if (is.null(rfx_working_parameter_prior_cov)) {
      sigma_alpha_init <- diag(1, num_rfx_components, num_rfx_components)
    } else {
      sigma_alpha_init <- expand_dims_2d_diag(
        rfx_working_parameter_prior_cov,
        num_rfx_components
      )
    }

    if (is.null(rfx_group_parameter_prior_cov)) {
      sigma_xi_init <- diag(1, num_rfx_components, num_rfx_components)
    } else {
      sigma_xi_init <- expand_dims_2d_diag(
        rfx_group_parameter_prior_cov,
        num_rfx_components
      )
    }

    sigma_xi_shape <- rfx_variance_prior_shape
    sigma_xi_scale <- rfx_variance_prior_scale

    # Random effects data structure and storage container
    rfx_dataset_train <- createRandomEffectsDataset(
      rfx_group_ids_train,
      rfx_basis_train,
      if (link_is_log) negbin_base_variance_weights else NULL
    )
    rfx_tracker_train <- createRandomEffectsTracker(rfx_group_ids_train)
    rfx_model <- createRandomEffectsModel(
      num_rfx_components,
      num_rfx_groups
    )
    rfx_model$set_working_parameter(alpha_init)
    rfx_model$set_group_parameters(xi_init)
    rfx_model$set_working_parameter_cov(sigma_alpha_init)
    rfx_model$set_group_parameter_cov(sigma_xi_init)
    rfx_model$set_variance_prior_shape(sigma_xi_shape)
    rfx_model$set_variance_prior_scale(sigma_xi_scale)
    rfx_samples <- createRandomEffectSamples(
      num_rfx_components,
      num_rfx_groups,
      rfx_tracker_train
    )
  }

  # Zero-part random effects: same groups and basis as the count part, priors
  # from zero_random_effects_params (default: the count part's), sampled on
  # the probit latent residual of the zero part with unit variance (no
  # Polya-Gamma weights)
  if (has_zero_rfx) {
    zero_rfx_prior <- .zinbZeroRfxPrior(
      rfx_params_updated,
      zero_random_effects_params,
      num_rfx_components,
      num_rfx_groups
    )
    zero_rfx_dataset_train <- createRandomEffectsDataset(
      rfx_group_ids_train,
      rfx_basis_train
    )
    zero_rfx_tracker_train <- createRandomEffectsTracker(rfx_group_ids_train)
    zero_rfx_model <- createRandomEffectsModel(
      num_rfx_components,
      num_rfx_groups
    )
    rootResetRandomEffectsModel(
      zero_rfx_model,
      zero_rfx_prior$alpha_init,
      zero_rfx_prior$xi_init,
      zero_rfx_prior$sigma_alpha_init,
      zero_rfx_prior$sigma_xi_init,
      zero_rfx_prior$sigma_xi_shape,
      zero_rfx_prior$sigma_xi_scale
    )
    zero_rfx_samples <- createRandomEffectSamples(
      num_rfx_components,
      num_rfx_groups,
      zero_rfx_tracker_train
    )
  }

  # Container of parameter samples
  num_actual_mcmc_iter <- num_mcmc * keep_every
  num_samples <- num_gfr + num_burnin + num_actual_mcmc_iter
  # Delete GFR samples from these containers after the fact if desired
  # num_retained_samples <- ifelse(keep_gfr, num_gfr, 0) + ifelse(keep_burnin, num_burnin, 0) + num_mcmc
  num_retained_samples <- num_gfr +
    ifelse(keep_burnin, num_burnin, 0) +
    num_mcmc * num_chains
  if (sample_sigma2_global) {
    global_var_samples <- rep(NA, num_retained_samples)
  }
  if (sample_sigma2_leaf) {
    leaf_scale_samples <- rep(NA, num_retained_samples)
  }
  if (link_is_log) {
    negbin_dispersion_samples <- rep(NA_real_, num_retained_samples)
  }
  if (outcome_is_zinb) {
    zero_forest_pred_train <- matrix(
      NA_real_,
      nrow(X_train),
      num_retained_samples
    )
    if (sample_sigma2_leaf_zero) {
      leaf_scale_zero_samples <- rep(NA_real_, num_retained_samples)
    }
  }
  if (link_is_cloglog) {
    cloglog_cutpoint_samples <- matrix(
      NA_real_,
      cloglog_num_categories - 1,
      num_retained_samples
    )
  }
  # Importance needs sampler state, but no observation-by-draw output buffer.
  if (include_mean_forest && !importance_only) {
    mean_forest_pred_train <- matrix(
      NA_real_,
      nrow(X_train),
      num_retained_samples
    )
  }
  if (include_variance_forest) {
    variance_forest_pred_train <- matrix(
      NA_real_,
      nrow(X_train),
      num_retained_samples
    )
  }
  if (selected_only && has_test) {
    selected_test_pred <- matrix(NA_real_, nrow(X_test), num_retained_samples)
  }
  sample_counter <- 0
  if (importance_only) importance_split_counts <- numeric(ncol(X_train))

  # Initialize the leaves of each tree in the mean forest
  if (include_mean_forest) {
    if (requires_basis) {
      # Handle the case in which we must initialize root values in a leaf basis regression
      # when init_val_mean != 0. To do this, we regress rep(init_val_mean, nrow(y_train))
      # on leaf_basis_train and use (coefs / num_trees_mean) as initial values
      if (abs(init_val_mean) > 0.00001) {
        init_val_y <- rep(init_val_mean, nrow(y_train))
        init_val_model <- lm(init_val_y ~ 0 + leaf_basis_train)
        init_values_mean_forest <- coef(init_val_model)
        if (any(is.na(init_values_mean_forest))) {
          init_values_mean_forest[which(is.na(init_values_mean_forest))] <- 0.
        }
      } else {
        init_values_mean_forest <- rep(init_val_mean, ncol(leaf_basis_train))
      }
    } else {
      init_values_mean_forest <- init_val_mean
    }
    active_forest_mean$prepare_for_sampler(
      forest_dataset_train,
      outcome_train,
      forest_model_mean,
      leaf_model_mean_forest,
      init_values_mean_forest
    )
  }

  # Initialize the leaves of each tree in the variance forest
  if (include_variance_forest) {
    active_forest_variance$prepare_for_sampler(
      forest_dataset_train,
      outcome_train,
      forest_model_variance,
      leaf_model_variance_forest,
      variance_forest_init
    )
  }

  # Initialize the leaves of each tree in the zero forest at 0 (probit scale,
  # around zero_outcome_mean)
  if (outcome_is_zinb) {
    active_forest_zero$prepare_for_sampler(
      forest_dataset_zero,
      outcome_zero,
      forest_model_zero,
      0,
      0.0
    )
  }

  # Initialize auxiliary data for cloglog
  if (link_is_cloglog) {
    ## Allocate auxiliary data
    train_size <- nrow(X_train)
    # Latent variable (Z in Alam et al (2025) notation)
    forest_dataset_train$add_auxiliary_dimension(train_size)
    # Forest predictions (eta in Alam et al (2025) notation)
    forest_dataset_train$add_auxiliary_dimension(train_size)
    # Log-scale non-cumulative cutpoint (gamma in Alam et al (2025) notation)
    forest_dataset_train$add_auxiliary_dimension(cloglog_num_categories - 1)
    # Exponentiated cumulative cutpoints (exp(c_k) in Alam et al (2025) notation)
    # This auxiliary series is designed so that the element stored at position `i`
    # corresponds to the sum of all exponentiated gamma_j values for j < i.
    # It has cloglog_num_categories elements instead of cloglog_num_categories - 1 because
    # even the largest categorical index has a valid value of sum_{j < i} exp(gamma_j)
    forest_dataset_train$add_auxiliary_dimension(cloglog_num_categories)

    ## Set initial values for auxiliary data
    # Initialize latent variables to zero (slot 0)
    for (i in 1:train_size) {
      forest_dataset_train$set_auxiliary_data_value(0, i - 1, 0.0)
    }
    # Initialize forest predictions to zero (slot 1)
    for (i in 1:train_size) {
      forest_dataset_train$set_auxiliary_data_value(1, i - 1, 0.0)
    }
    # Initialize log-scale cutpoints to 0
    initial_gamma <- rep(0.0, cloglog_num_categories - 1)
    for (i in seq_along(initial_gamma)) {
      forest_dataset_train$set_auxiliary_data_value(2, i - 1, initial_gamma[i])
    }
    # Convert to cumulative exponentiated cutpoints directly in C++
    ordinal_sampler_update_cumsum_exp_cpp(
      ordinal_sampler,
      forest_dataset_train$data_ptr
    )
  }

  # Gibbs steps of the count (negative-binomial) outcome, shared by the GFR
  # and MCMC loops. See R/negbin.R for the derivation.
  if (link_is_log) {
    # Polya-Gamma step: omega | eta, r, then the pseudo-outcome and the
    # variance weights seen by the forest and random-effects samplers
    negbin_latent_update <- function(forest_pred = NULL, rfx_pred = NULL) {
      if (is.null(forest_pred)) {
        forest_pred <- active_forest_mean$predict(forest_dataset_train)
      }
      if (is.null(rfx_pred)) {
        rfx_pred <- if (has_rfx) {
          rfx_model$predict(rfx_dataset_train, rfx_tracker_train)
        } else {
          0
        }
      }
      log_dispersion <- log(current_dispersion)
      psi <- negbin_offset_train + y_bar_train + forest_pred + rfx_pred -
        log_dispersion
      if (outcome_is_zinb) {
        # Only the rows with s_i = 0 follow the negative binomial. The others
        # get an infinite variance weight (1 / w = 0), so they drop out of the
        # forest and random-effects sufficient statistics; their residual is
        # set to 0 to keep it finite.
        active <- which(current_structural_zero == 0L)
        omega <- sample_polya_gamma_cpp(
          negbin_y[active] + current_dispersion,
          psi[active],
          rng$rng_ptr
        )
        residual <- numeric(length(negbin_y))
        residual[active] <- (negbin_y[active] - current_dispersion) /
          (2 * omega) -
          psi[active]
        outcome_train$update_data(residual)
        variance_weights <- rep(Inf, length(negbin_y))
        variance_weights[active] <- 1 / omega
        forest_dataset_train$update_variance_weights(
          variance_weights,
          exponentiate = FALSE
        )
        if (has_rfx) {
          rfx_dataset_train$update_variance_weights(
            variance_weights,
            exponentiate = FALSE
          )
        }
        return(invisible(NULL))
      }
      omega <- sample_polya_gamma_cpp(
        negbin_y + current_dispersion,
        psi,
        rng$rng_ptr
      )
      # (y - r) / (2 omega) ~ N(psi, 1 / omega); forest and rfx fit psi net of
      # the offset, the centering constant and -log(r)
      pseudo_outcome <- (negbin_y - current_dispersion) / (2 * omega) -
        negbin_offset_train -
        y_bar_train +
        log_dispersion
      outcome_train$update_data(pseudo_outcome - forest_pred - rfx_pred)
      variance_weights <- negbin_base_variance_weights / omega
      forest_dataset_train$update_variance_weights(
        variance_weights,
        exponentiate = FALSE
      )
      if (has_rfx) {
        rfx_dataset_train$update_variance_weights(
          variance_weights,
          exponentiate = FALSE
        )
      }
    }
    # Dispersion step: r | y, mu
    negbin_dispersion_update <- function() {
      if (!isTRUE(sample_negbin_dispersion)) {
        return(current_dispersion)
      }
      forest_pred <- forest_model_mean$get_cached_forest_predictions()
      rfx_pred <- if (has_rfx) {
        rfx_model$predict(rfx_dataset_train, rfx_tracker_train)
      } else {
        0
      }
      mu <- exp(negbin_offset_train + y_bar_train + forest_pred + rfx_pred)
      if (outcome_is_zinb) {
        # r | y, mu over the rows that follow the negative binomial
        active <- current_structural_zero == 0L
        return(.negbinSampleDispersion(
          current_dispersion,
          negbin_y[active],
          mu[active],
          negbin_dispersion_prior_shape,
          negbin_dispersion_prior_rate
        ))
      }
      .negbinSampleDispersion(
        current_dispersion,
        negbin_y,
        mu,
        negbin_dispersion_prior_shape,
        negbin_dispersion_prior_rate
      )
    }
  }

  # Zero-inflated family: structural-zero indicators and the probit zero
  # forest (steps 1 and 2 of the scheme in R/negbin.R). Returns the current
  # count-forest and random-effects predictions for the Polya-Gamma step,
  # which these steps do not change.
  if (outcome_is_zinb) {
    zinb_zero_update <- function(gfr, keep_sample) {
      forest_pred <- active_forest_mean$predict(forest_dataset_train)
      rfx_pred <- if (has_rfx) {
        rfx_model$predict(rfx_dataset_train, rfx_tracker_train)
      } else {
        0
      }
      eta_count <- negbin_offset_train + y_bar_train + forest_pred + rfx_pred
      zero_pred <- active_forest_zero$predict(forest_dataset_zero)
      eta_zero <- zero_outcome_mean + zero_pred
      if (has_zero_rfx) {
        zero_rfx_pred <- zero_rfx_model$predict(
          zero_rfx_dataset_train,
          zero_rfx_tracker_train
        )
        eta_zero <- eta_zero + zero_rfx_pred
      }
      s <- .zinbSampleStructuralZeros(
        negbin_y,
        eta_count,
        eta_zero,
        current_dispersion
      )
      # Albert and Chib latent z_i ~ N(eta_zero_i, 1) truncated to the side
      # given by s_i, as the probit outcome does
      z <- .probitLatentDraw(eta_zero, s)
      if (has_zero_rfx) {
        outcome_zero$update_data(z - zero_outcome_mean - zero_pred - zero_rfx_pred)
      } else {
        outcome_zero$update_data(z - zero_outcome_mean - zero_pred)
      }
      forest_model_zero$sample_one_iteration(
        forest_dataset = forest_dataset_zero,
        residual = outcome_zero,
        forest_samples = forest_samples_zero,
        active_forest = active_forest_zero,
        rng = rng,
        forest_model_config = forest_model_config_zero,
        global_model_config = global_model_config_zero,
        num_threads = num_threads,
        keep_forest = keep_sample,
        gfr = gfr
      )
      if (keep_sample) {
        zero_forest_pred_train[,
          sample_counter
        ] <<- forest_model_zero$get_cached_forest_predictions()
      }
      if (sample_sigma2_leaf_zero) {
        leaf_scale_zero <- sampleLeafVarianceOneIteration(
          active_forest_zero,
          rng,
          a_leaf_zero,
          b_leaf_zero
        )
        if (keep_sample) {
          leaf_scale_zero_samples[sample_counter] <<- leaf_scale_zero
        }
        forest_model_config_zero$update_leaf_model_scale(
          as.matrix(leaf_scale_zero)
        )
      }
      # Zero-part random effects on the latent residual z - c0 - g(X) - u,
      # with the latent error variance fixed at 1
      if (has_zero_rfx) {
        zero_rfx_model$sample_random_effect(
          zero_rfx_dataset_train,
          outcome_zero,
          zero_rfx_tracker_train,
          zero_rfx_samples,
          keep_sample,
          1.0,
          rng
        )
      }
      current_structural_zero <<- s
      list(forest_pred = forest_pred, rfx_pred = rfx_pred)
    }
  }

  # Run GFR (warm start) if specified
  if (num_gfr > 0) {
    for (i in 1:num_gfr) {
      # Keep all GFR samples at this stage -- remove from ForestSamples after MCMC
      # keep_sample <- ifelse(keep_gfr, TRUE, FALSE)
      keep_sample <- TRUE
      if (keep_sample) {
        sample_counter <- sample_counter + 1
      }
      # Print progress
      if (verbose) {
        if ((i %% 10 == 0) || (i == num_gfr)) {
          cat(
            "Sampling",
            i,
            "out of",
            num_gfr,
            "XBART (grow-from-root) draws\n"
          )
        }
      }

      if (include_mean_forest) {
        if (link_is_probit) {
          # Sample latent probit variable, z | -
          # outcome_pred is the centered forest prediction (not including y_bar_train).
          # The truncated normal mean is outcome_pred + y_bar_train (the full eta on the probit scale).
          # The residual stored is z - y_bar_train - outcome_pred so the forest sees a
          # zero-centered signal and the prior shrinkage toward 0 is well-calibrated.
          outcome_pred <- active_forest_mean$predict(
            forest_dataset_train
          )
          if (has_rfx) {
            rfx_pred <- rfx_model$predict(
              rfx_dataset_train,
              rfx_tracker_train
            )
            outcome_pred <- outcome_pred + rfx_pred
          }
          eta_pred <- outcome_pred + y_bar_train
          mu0 <- eta_pred[y_train == 0]
          mu1 <- eta_pred[y_train == 1]
          u0 <- runif(sum(y_train == 0), 0, pnorm(0 - mu0))
          u1 <- runif(sum(y_train == 1), pnorm(0 - mu1), 1)
          resid_train[y_train == 0] <- mu0 + qnorm(u0)
          resid_train[y_train == 1] <- mu1 + qnorm(u1)

          # Update outcome: center z by y_bar_train before passing to forest
          outcome_train$update_data(resid_train - y_bar_train - outcome_pred)
        }
        if (outcome_is_zinb) {
          zinb_state <- zinb_zero_update(gfr = TRUE, keep_sample = keep_sample)
          negbin_latent_update(zinb_state$forest_pred, zinb_state$rfx_pred)
        } else if (link_is_log) {
          negbin_latent_update()
        }

        # Sample mean forest
        forest_model_mean$sample_one_iteration(
          forest_dataset = forest_dataset_train,
          residual = outcome_train,
          forest_samples = forest_samples_mean,
          active_forest = active_forest_mean,
          rng = rng,
          forest_model_config = forest_model_config_mean,
          global_model_config = global_model_config,
          num_threads = num_threads,
          keep_forest = keep_sample,
          gfr = TRUE
        )

        if (importance_only && keep_gfr) {
          importance_split_counts <- importance_split_counts +
            active_forest_mean$get_forest_split_counts(ncol(X_train))
        }

        # Cache train set predictions since they are already computed during sampling
        if (keep_sample && !importance_only) {
          mean_forest_pred_train[,
            sample_counter
          ] <- forest_model_mean$get_cached_forest_predictions()
        }

        # Additional Gibbs updates needed for the cloglog model
        if (link_is_cloglog) {
          # Update auxiliary data to current forest predictions
          forest_pred_current <- forest_model_mean$get_cached_forest_predictions()
          for (i in 1:train_size) {
            forest_dataset_train$set_auxiliary_data_value(
              1,
              i - 1,
              forest_pred_current[i]
            )
          }

          # Sample latent z_i's using truncated exponential
          ordinal_sampler_update_latent_variables_cpp(
            ordinal_sampler,
            forest_dataset_train$data_ptr,
            outcome_train$data_ptr,
            rng$rng_ptr
          )

          # Sample gamma parameters (cutpoints)
          ordinal_sampler_update_gamma_params_cpp(
            ordinal_sampler,
            forest_dataset_train$data_ptr,
            outcome_train$data_ptr,
            cloglog_forest_shape,
            cloglog_forest_rate,
            cloglog_cutpoint_0,
            rng$rng_ptr
          )

          # Update cumulative sum of exp(gamma) values
          ordinal_sampler_update_cumsum_exp_cpp(
            ordinal_sampler,
            forest_dataset_train$data_ptr
          )

          # Retain cutpoint draw
          if (keep_sample) {
            cloglog_cutpoints <- forest_dataset_train$get_auxiliary_data_vector(
              2
            )
            cloglog_cutpoint_samples[, sample_counter] <- cloglog_cutpoints
          }
        }
      }
      if (include_variance_forest) {
        forest_model_variance$sample_one_iteration(
          forest_dataset = forest_dataset_train,
          residual = outcome_train,
          forest_samples = forest_samples_variance,
          active_forest = active_forest_variance,
          rng = rng,
          forest_model_config = forest_model_config_variance,
          global_model_config = global_model_config,
          num_threads = num_threads,
          keep_forest = keep_sample,
          gfr = TRUE
        )

        # Cache train set predictions since they are already computed during sampling
        if (keep_sample) {
          variance_forest_pred_train[,
            sample_counter
          ] <- forest_model_variance$get_cached_forest_predictions()
        }
      }
      if (sample_sigma2_global) {
        current_sigma2 <- sampleGlobalErrorVarianceOneIteration(
          outcome_train,
          forest_dataset_train,
          rng,
          a_global,
          b_global
        )
        if (keep_sample) {
          global_var_samples[sample_counter] <- current_sigma2
        }
        global_model_config$update_global_error_variance(current_sigma2)
      }
      if (sample_sigma2_leaf) {
        leaf_scale_double <- sampleLeafVarianceOneIteration(
          active_forest_mean,
          rng,
          a_leaf,
          b_leaf
        )
        current_leaf_scale <- as.matrix(leaf_scale_double)
        if (keep_sample) {
          leaf_scale_samples[sample_counter] <- leaf_scale_double
        }
        forest_model_config_mean$update_leaf_model_scale(
          current_leaf_scale
        )
      }
      if (has_rfx) {
        rfx_model$sample_random_effect(
          rfx_dataset_train,
          outcome_train,
          rfx_tracker_train,
          rfx_samples,
          keep_sample,
          current_sigma2,
          rng
        )
      }
      if (link_is_log) {
        current_dispersion <- negbin_dispersion_update()
        if (keep_sample) {
          negbin_dispersion_samples[sample_counter] <- current_dispersion
        }
      }
    }
  }

  # Run MCMC
  if (num_burnin + num_mcmc > 0) {
    for (chain_num in 1:num_chains) {
      if (verbose) {
        cat("Sampling chain", chain_num, "of", num_chains, "\n")
      }
      if (num_gfr > 0) {
        # Reset state of active_forest and forest_model based on a previous GFR sample
        forest_ind <- num_gfr - chain_num
        if (include_mean_forest) {
          resetActiveForest(
            active_forest_mean,
            forest_samples_mean,
            forest_ind
          )
          resetForestModel(
            forest_model_mean,
            active_forest_mean,
            forest_dataset_train,
            outcome_train,
            TRUE
          )
          if (sample_sigma2_leaf) {
            leaf_scale_double <- leaf_scale_samples[forest_ind + 1]
            current_leaf_scale <- as.matrix(leaf_scale_double)
            forest_model_config_mean$update_leaf_model_scale(
              current_leaf_scale
            )
          }
          if (link_is_cloglog) {
            # Restore ordinal labels corrupted by resetForestModel's
            # residual adjustment (outcome stores category labels, not residuals)
            outcome_train$update_data(resid_train)
            # We can reset cutpoints from warm-start since cutpoints are retained
            current_cutpoints <- cloglog_cutpoint_samples[, forest_ind + 1]
            for (i in seq_along(current_cutpoints)) {
              forest_dataset_train$set_auxiliary_data_value(
                2,
                i - 1,
                current_cutpoints[i]
              )
            }
            ordinal_sampler_update_cumsum_exp_cpp(
              ordinal_sampler,
              forest_dataset_train$data_ptr
            )
            # Re-predict from the reconstituted active forest
            active_forest_preds <- active_forest_mean$predict(
              forest_dataset_train
            )
            for (i in 1:train_size) {
              forest_dataset_train$set_auxiliary_data_value(
                1,
                i - 1,
                active_forest_preds[i]
              )
              # Latent variables must be reset to 0 and burnt in
              forest_dataset_train$set_auxiliary_data_value(0, i - 1, 0.0)
            }
          }
        }
        if (include_variance_forest) {
          resetActiveForest(
            active_forest_variance,
            forest_samples_variance,
            forest_ind
          )
          resetForestModel(
            forest_model_variance,
            active_forest_variance,
            forest_dataset_train,
            outcome_train,
            FALSE
          )
        }
        if (has_rfx) {
          resetRandomEffectsModel(
            rfx_model,
            rfx_samples,
            forest_ind,
            sigma_alpha_init
          )
          resetRandomEffectsTracker(
            rfx_tracker_train,
            rfx_model,
            rfx_dataset_train,
            outcome_train,
            rfx_samples
          )
        }
        if (sample_sigma2_global) {
          current_sigma2 <- global_var_samples[forest_ind + 1]
          global_model_config$update_global_error_variance(
            current_sigma2
          )
        }
        if (link_is_log) {
          current_dispersion <- negbin_dispersion_samples[forest_ind + 1]
        }
        if (outcome_is_zinb) {
          resetActiveForest(
            active_forest_zero,
            forest_samples_zero,
            forest_ind
          )
          resetForestModel(
            forest_model_zero,
            active_forest_zero,
            forest_dataset_zero,
            outcome_zero,
            TRUE
          )
          if (sample_sigma2_leaf_zero) {
            forest_model_config_zero$update_leaf_model_scale(
              as.matrix(leaf_scale_zero_samples[forest_ind + 1])
            )
          }
          if (has_zero_rfx) {
            resetRandomEffectsModel(
              zero_rfx_model,
              zero_rfx_samples,
              forest_ind,
              zero_rfx_prior$sigma_alpha_init
            )
            resetRandomEffectsTracker(
              zero_rfx_tracker_train,
              zero_rfx_model,
              zero_rfx_dataset_train,
              outcome_zero,
              zero_rfx_samples
            )
          }
        }
      } else if (has_prev_model) {
        warmstart_index <- ifelse(
          previous_model_decrement,
          previous_model_warmstart_sample_num - chain_num + 1,
          previous_model_warmstart_sample_num
        )
        if (include_mean_forest) {
          resetActiveForest(
            active_forest_mean,
            previous_forest_samples_mean,
            warmstart_index - 1
          )
          resetForestModel(
            forest_model_mean,
            active_forest_mean,
            forest_dataset_train,
            outcome_train,
            TRUE
          )
          if (
            sample_sigma2_leaf &&
              (!is.null(previous_leaf_var_samples))
          ) {
            leaf_scale_double <- previous_leaf_var_samples[
              warmstart_index
            ]
            current_leaf_scale <- as.matrix(leaf_scale_double)
            forest_model_config_mean$update_leaf_model_scale(
              current_leaf_scale
            )
          }
          if (link_is_cloglog) {
            # Restore ordinal labels corrupted by resetForestModel's
            # residual adjustment (outcome stores category labels, not residuals)
            outcome_train$update_data(resid_train)
            # We can reset cutpoints from warm-start since cutpoints are retained
            current_cutpoints <- previous_cloglog_cutpoint_samples[,
              warmstart_index
            ]
            for (i in seq_along(current_cutpoints)) {
              forest_dataset_train$set_auxiliary_data_value(
                2,
                i - 1,
                current_cutpoints[i]
              )
            }
            ordinal_sampler_update_cumsum_exp_cpp(
              ordinal_sampler,
              forest_dataset_train$data_ptr
            )
            # Re-predict from the reconstituted active forest
            active_forest_preds <- active_forest_mean$predict(
              forest_dataset_train
            )
            for (i in 1:train_size) {
              forest_dataset_train$set_auxiliary_data_value(
                1,
                i - 1,
                active_forest_preds[i]
              )
              # Latent variables must be reset to 0 and burnt in
              forest_dataset_train$set_auxiliary_data_value(0, i - 1, 0.0)
            }
          }
        }
        if (include_variance_forest) {
          resetActiveForest(
            active_forest_variance,
            previous_forest_samples_variance,
            warmstart_index - 1
          )
          resetForestModel(
            forest_model_variance,
            active_forest_variance,
            forest_dataset_train,
            outcome_train,
            FALSE
          )
        }
        if (has_rfx) {
          if (is.null(previous_rfx_samples)) {
            warning(
              "`previous_model_json` did not have any random effects samples, so the RFX sampler will be run from scratch while the forests and any other parameters are warm started"
            )
            rootResetRandomEffectsModel(
              rfx_model,
              alpha_init,
              xi_init,
              sigma_alpha_init,
              sigma_xi_init,
              sigma_xi_shape,
              sigma_xi_scale
            )
            rootResetRandomEffectsTracker(
              rfx_tracker_train,
              rfx_model,
              rfx_dataset_train,
              outcome_train
            )
          } else {
            resetRandomEffectsModel(
              rfx_model,
              previous_rfx_samples,
              warmstart_index - 1,
              sigma_alpha_init
            )
            resetRandomEffectsTracker(
              rfx_tracker_train,
              rfx_model,
              rfx_dataset_train,
              outcome_train,
              rfx_samples
            )
          }
        }
        if (sample_sigma2_global) {
          if (!is.null(previous_global_var_samples)) {
            current_sigma2 <- previous_global_var_samples[
              warmstart_index
            ]
            global_model_config$update_global_error_variance(
              current_sigma2
            )
          }
        }
        if (link_is_log && !is.null(previous_negbin_dispersion_samples)) {
          current_dispersion <- previous_negbin_dispersion_samples[
            warmstart_index
          ]
        }
        if (outcome_is_zinb) {
          resetActiveForest(
            active_forest_zero,
            previous_forest_samples_zero,
            warmstart_index - 1
          )
          resetForestModel(
            forest_model_zero,
            active_forest_zero,
            forest_dataset_zero,
            outcome_zero,
            TRUE
          )
          if (
            sample_sigma2_leaf_zero && !is.null(previous_leaf_var_zero_samples)
          ) {
            forest_model_config_zero$update_leaf_model_scale(
              as.matrix(previous_leaf_var_zero_samples[warmstart_index])
            )
          }
          if (has_zero_rfx) {
            if (is.null(previous_zero_rfx_samples)) {
              warning(
                "`previous_model_json` did not have any zero-part random effects samples, so the zero-part RFX sampler will be run from scratch while the forests and any other parameters are warm started"
              )
              rootResetRandomEffectsModel(
                zero_rfx_model,
                zero_rfx_prior$alpha_init,
                zero_rfx_prior$xi_init,
                zero_rfx_prior$sigma_alpha_init,
                zero_rfx_prior$sigma_xi_init,
                zero_rfx_prior$sigma_xi_shape,
                zero_rfx_prior$sigma_xi_scale
              )
              rootResetRandomEffectsTracker(
                zero_rfx_tracker_train,
                zero_rfx_model,
                zero_rfx_dataset_train,
                outcome_zero
              )
            } else {
              resetRandomEffectsModel(
                zero_rfx_model,
                previous_zero_rfx_samples,
                warmstart_index - 1,
                zero_rfx_prior$sigma_alpha_init
              )
              resetRandomEffectsTracker(
                zero_rfx_tracker_train,
                zero_rfx_model,
                zero_rfx_dataset_train,
                outcome_zero,
                zero_rfx_samples
              )
            }
          }
        }
      } else {
        if (include_mean_forest) {
          resetActiveForest(active_forest_mean)
          active_forest_mean$set_root_leaves(
            init_values_mean_forest / num_trees_mean
          )
          resetForestModel(
            forest_model_mean,
            active_forest_mean,
            forest_dataset_train,
            outcome_train,
            TRUE
          )
          if (sample_sigma2_leaf) {
            current_leaf_scale <- as.matrix(sigma2_leaf_init)
            forest_model_config_mean$update_leaf_model_scale(
              current_leaf_scale
            )
          }
          if (link_is_cloglog) {
            # Restore ordinal labels corrupted by resetForestModel's
            # residual adjustment (outcome stores category labels, not residuals)
            outcome_train$update_data(resid_train)
            # Reset all cloglog parameters to default values
            for (i in 1:train_size) {
              forest_dataset_train$set_auxiliary_data_value(0, i - 1, 0.0)
              forest_dataset_train$set_auxiliary_data_value(1, i - 1, 0.0)
            }
            # Initialize log-scale cutpoints to 0
            initial_gamma <- rep(0.0, cloglog_num_categories - 1)
            for (i in seq_along(initial_gamma)) {
              forest_dataset_train$set_auxiliary_data_value(
                2,
                i - 1,
                initial_gamma[i]
              )
            }
            # Convert to cumulative exponentiated cutpoints directly in C++
            ordinal_sampler_update_cumsum_exp_cpp(
              ordinal_sampler,
              forest_dataset_train$data_ptr
            )
          }
        }
        if (include_variance_forest) {
          resetActiveForest(active_forest_variance)
          active_forest_variance$set_root_leaves(
            log(variance_forest_init) / num_trees_variance
          )
          resetForestModel(
            forest_model_variance,
            active_forest_variance,
            forest_dataset_train,
            outcome_train,
            FALSE
          )
        }
        if (has_rfx) {
          rootResetRandomEffectsModel(
            rfx_model,
            alpha_init,
            xi_init,
            sigma_alpha_init,
            sigma_xi_init,
            sigma_xi_shape,
            sigma_xi_scale
          )
          rootResetRandomEffectsTracker(
            rfx_tracker_train,
            rfx_model,
            rfx_dataset_train,
            outcome_train
          )
        }
        if (sample_sigma2_global) {
          current_sigma2 <- sigma2_init
          global_model_config$update_global_error_variance(
            current_sigma2
          )
        }
        if (link_is_log) {
          current_dispersion <- negbin_dispersion_init
        }
        if (outcome_is_zinb) {
          resetActiveForest(active_forest_zero)
          active_forest_zero$set_root_leaves(0)
          resetForestModel(
            forest_model_zero,
            active_forest_zero,
            forest_dataset_zero,
            outcome_zero,
            TRUE
          )
          if (sample_sigma2_leaf_zero) {
            forest_model_config_zero$update_leaf_model_scale(
              as.matrix(sigma2_leaf_zero_init)
            )
          }
          if (has_zero_rfx) {
            rootResetRandomEffectsModel(
              zero_rfx_model,
              zero_rfx_prior$alpha_init,
              zero_rfx_prior$xi_init,
              zero_rfx_prior$sigma_alpha_init,
              zero_rfx_prior$sigma_xi_init,
              zero_rfx_prior$sigma_xi_shape,
              zero_rfx_prior$sigma_xi_scale
            )
            rootResetRandomEffectsTracker(
              zero_rfx_tracker_train,
              zero_rfx_model,
              zero_rfx_dataset_train,
              outcome_zero
            )
          }
        }
      }
      for (i in (num_gfr + 1):num_samples) {
        is_mcmc <- i > (num_gfr + num_burnin)
        if (is_mcmc) {
          mcmc_counter <- i - (num_gfr + num_burnin)
          if (mcmc_counter %% keep_every == 0) {
            keep_sample <- TRUE
          } else {
            keep_sample <- FALSE
          }
        } else {
          if (keep_burnin) {
            keep_sample <- TRUE
          } else {
            keep_sample <- FALSE
          }
        }
        if (keep_sample) {
          sample_counter <- sample_counter + 1
        }
        # Print progress
        if (verbose) {
          if (num_burnin > 0 && !is_mcmc) {
            if (
              ((i - num_gfr) %% 100 == 0) ||
                ((i - num_gfr) == num_burnin)
            ) {
              cat(
                "Sampling",
                i - num_gfr,
                "out of",
                num_burnin,
                "BART burn-in draws; Chain number ",
                chain_num,
                "\n"
              )
            }
          }
          if (num_mcmc > 0 && is_mcmc) {
            raw_iter <- i - num_gfr - num_burnin
            if ((raw_iter %% 100 == 0) || (i == num_samples)) {
              if (keep_every == 1) {
                cat(
                  "Sampling",
                  raw_iter,
                  "out of",
                  num_mcmc,
                  "BART MCMC draws; Chain number ",
                  chain_num,
                  "\n"
                )
              } else {
                cat(
                  "Sampling raw draw",
                  raw_iter,
                  "of",
                  num_actual_mcmc_iter,
                  "BART MCMC draws (thinning by",
                  keep_every,
                  ":",
                  raw_iter %/% keep_every,
                  "of",
                  num_mcmc,
                  "retained); Chain number ",
                  chain_num,
                  "\n"
                )
              }
            }
          }
        }

        if (include_mean_forest) {
          if (link_is_probit) {
            # Sample latent probit variable, z | -
            outcome_pred <- active_forest_mean$predict(
              forest_dataset_train
            )
            if (has_rfx) {
              rfx_pred <- rfx_model$predict(
                rfx_dataset_train,
                rfx_tracker_train
              )
              outcome_pred <- outcome_pred + rfx_pred
            }
            eta_pred <- outcome_pred + y_bar_train
            mu0 <- eta_pred[y_train == 0]
            mu1 <- eta_pred[y_train == 1]
            u0 <- runif(sum(y_train == 0), 0, pnorm(0 - mu0))
            u1 <- runif(sum(y_train == 1), pnorm(0 - mu1), 1)
            resid_train[y_train == 0] <- mu0 + qnorm(u0)
            resid_train[y_train == 1] <- mu1 + qnorm(u1)

            # Update outcome: center z by y_bar_train before passing to forest
            outcome_train$update_data(
              resid_train - y_bar_train - outcome_pred
            )
          }
          if (outcome_is_zinb) {
            zinb_state <- zinb_zero_update(
              gfr = FALSE,
              keep_sample = keep_sample
            )
            negbin_latent_update(zinb_state$forest_pred, zinb_state$rfx_pred)
          } else if (link_is_log) {
            negbin_latent_update()
          }

          forest_model_mean$sample_one_iteration(
            forest_dataset = forest_dataset_train,
            residual = outcome_train,
            forest_samples = forest_samples_mean,
            active_forest = active_forest_mean,
            rng = rng,
            forest_model_config = forest_model_config_mean,
            global_model_config = global_model_config,
            num_threads = num_threads,
            keep_forest = keep_sample && !importance_only &&
              (!selected_only || selected_mask[[sample_counter - num_gfr]]),
            gfr = FALSE
          )
          if (selected_only && has_test && keep_sample) {
            selected_test_pred[, sample_counter] <- active_forest_mean$predict(forest_dataset_test)
          }
          if (importance_only && keep_sample) {
            importance_split_counts <- importance_split_counts +
              active_forest_mean$get_forest_split_counts(ncol(X_train))
          }

          # Cache train set predictions since they are already computed during sampling
          if (keep_sample && !importance_only) {
            mean_forest_pred_train[,
              sample_counter
            ] <- forest_model_mean$get_cached_forest_predictions()
          }

          # Additional Gibbs updates needed for the cloglog model
          if (link_is_cloglog) {
            # Update auxiliary data to current forest predictions
            forest_pred_current <- forest_model_mean$get_cached_forest_predictions()
            for (i in 1:train_size) {
              forest_dataset_train$set_auxiliary_data_value(
                1,
                i - 1,
                forest_pred_current[i]
              )
            }

            # Sample latent z_i's using truncated exponential
            ordinal_sampler_update_latent_variables_cpp(
              ordinal_sampler,
              forest_dataset_train$data_ptr,
              outcome_train$data_ptr,
              rng$rng_ptr
            )

            # Sample gamma parameters (cutpoints)
            ordinal_sampler_update_gamma_params_cpp(
              ordinal_sampler,
              forest_dataset_train$data_ptr,
              outcome_train$data_ptr,
              cloglog_forest_shape,
              cloglog_forest_rate,
              cloglog_cutpoint_0,
              rng$rng_ptr
            )

            # Update cumulative sum of exp(gamma) values
            ordinal_sampler_update_cumsum_exp_cpp(
              ordinal_sampler,
              forest_dataset_train$data_ptr
            )

            # Retain cutpoint draw
            if (keep_sample) {
              cloglog_cutpoints <- forest_dataset_train$get_auxiliary_data_vector(
                2
              )
              cloglog_cutpoint_samples[, sample_counter] <- cloglog_cutpoints
            }
          }
        }
        if (include_variance_forest) {
          forest_model_variance$sample_one_iteration(
            forest_dataset = forest_dataset_train,
            residual = outcome_train,
            forest_samples = forest_samples_variance,
            active_forest = active_forest_variance,
            rng = rng,
            forest_model_config = forest_model_config_variance,
            global_model_config = global_model_config,
            num_threads = num_threads,
            keep_forest = keep_sample,
            gfr = FALSE
          )

          # Cache train set predictions since they are already computed during sampling
          if (keep_sample) {
            variance_forest_pred_train[,
              sample_counter
            ] <- forest_model_variance$get_cached_forest_predictions()
          }
        }
        if (sample_sigma2_global) {
          current_sigma2 <- sampleGlobalErrorVarianceOneIteration(
            outcome_train,
            forest_dataset_train,
            rng,
            a_global,
            b_global
          )
          if (keep_sample) {
            global_var_samples[sample_counter] <- current_sigma2
          }
          global_model_config$update_global_error_variance(
            current_sigma2
          )
        }
        if (sample_sigma2_leaf) {
          leaf_scale_double <- sampleLeafVarianceOneIteration(
            active_forest_mean,
            rng,
            a_leaf,
            b_leaf
          )
          current_leaf_scale <- as.matrix(leaf_scale_double)
          if (keep_sample) {
            leaf_scale_samples[sample_counter] <- leaf_scale_double
          }
          forest_model_config_mean$update_leaf_model_scale(
            current_leaf_scale
          )
        }
        if (has_rfx) {
          rfx_model$sample_random_effect(
            rfx_dataset_train,
            outcome_train,
            rfx_tracker_train,
            rfx_samples,
            keep_sample,
            current_sigma2,
            rng
          )
        }
        if (link_is_log) {
          current_dispersion <- negbin_dispersion_update()
          if (keep_sample) {
            negbin_dispersion_samples[sample_counter] <- current_dispersion
          }
        }
      }
    }
  }

  # Remove GFR samples if they are not to be retained
  if ((!keep_gfr) && (num_gfr > 0)) {
    for (i in 1:num_gfr) {
      if (include_mean_forest) {
        forest_samples_mean$delete_sample(0)
      }
      if (include_variance_forest) {
        forest_samples_variance$delete_sample(0)
      }
      if (outcome_is_zinb) {
        forest_samples_zero$delete_sample(0)
      }
      if (has_rfx) {
        rfx_samples$delete_sample(0)
      }
      if (has_zero_rfx) {
        zero_rfx_samples$delete_sample(0)
      }
    }
    if (include_mean_forest && !importance_only) {
      mean_forest_pred_train <- mean_forest_pred_train[,
        (num_gfr + 1):ncol(mean_forest_pred_train)
      ]
      if (link_is_cloglog) {
        cloglog_cutpoint_samples <- cloglog_cutpoint_samples[,
          (num_gfr + 1):ncol(cloglog_cutpoint_samples),
          drop = FALSE
        ]
      }
    }
    if (include_variance_forest) {
      variance_forest_pred_train <- variance_forest_pred_train[,
        (num_gfr + 1):ncol(variance_forest_pred_train)
      ]
    }
    if (sample_sigma2_global) {
      global_var_samples <- global_var_samples[
        (num_gfr + 1):length(global_var_samples)
      ]
    }
    if (sample_sigma2_leaf) {
      leaf_scale_samples <- leaf_scale_samples[
        (num_gfr + 1):length(leaf_scale_samples)
      ]
    }
    if (link_is_log) {
      negbin_dispersion_samples <- negbin_dispersion_samples[
        (num_gfr + 1):length(negbin_dispersion_samples)
      ]
    }
    if (outcome_is_zinb) {
      zero_forest_pred_train <- zero_forest_pred_train[,
        (num_gfr + 1):ncol(zero_forest_pred_train),
        drop = FALSE
      ]
      if (sample_sigma2_leaf_zero) {
        leaf_scale_zero_samples <- leaf_scale_zero_samples[
          (num_gfr + 1):length(leaf_scale_zero_samples)
        ]
      }
    }
    num_retained_samples <- num_retained_samples - num_gfr
  }

  # Mean forest predictions
  if (include_mean_forest && !importance_only) {
    # y_hat_train <- forest_samples_mean$predict(forest_dataset_train)*y_std_train + y_bar_train
    y_hat_train <- mean_forest_pred_train * y_std_train + y_bar_train
    if (has_test) {
      y_hat_test <- if (selected_only) {
        selected_test_pred[, num_gfr + seq_len(num_mcmc), drop = FALSE] * y_std_train + y_bar_train
      } else {
        forest_samples_mean$predict(forest_dataset_test) * y_std_train + y_bar_train
      }
    }
  }

  # Variance forest predictions
  if (include_variance_forest) {
    # sigma2_x_hat_train <- forest_samples_variance$predict(forest_dataset_train)
    sigma2_x_hat_train <- exp(variance_forest_pred_train)
    if (has_test) {
      sigma2_x_hat_test <- forest_samples_variance$predict(
        forest_dataset_test
      )
    }
  }

  # Random effects predictions
  if (has_rfx && !importance_only) {
    rfx_preds_train <- rfx_samples$predict(
      rfx_group_ids_train,
      rfx_basis_train
    ) *
      y_std_train
    y_hat_train <- y_hat_train + rfx_preds_train
  }
  if ((has_rfx_test) && (has_test)) {
    rfx_preds_test <- predictRandomEffectsForLabels(
      rfx_samples,
      group_ids_factor_test,
      rfx_basis_test,
      unseen_groups = rfx_unseen_groups,
      label = "rfx_group_ids_test"
    ) *
      y_std_train
    y_hat_test <- y_hat_test + rfx_preds_test
  }

  # Count outcome: y_hat holds the full linear predictor log(mu), offset included
  if (link_is_log && !importance_only) {
    y_hat_train <- y_hat_train + negbin_offset_train
    if (has_test) {
      y_hat_test <- y_hat_test + negbin_offset_test
    }
  }

  # Zero-inflated family: probit index of the zero probability, pi = Phi(.)
  if (outcome_is_zinb) {
    zero_hat_train <- zero_forest_pred_train + zero_outcome_mean
    if (has_test) {
      zero_hat_test <- forest_samples_zero$predict(forest_dataset_test) +
        zero_outcome_mean
    }
  }
  # Zero-part random effects enter the probit index on both sets
  if (has_zero_rfx) {
    zero_rfx_preds_train <- zero_rfx_samples$predict(
      rfx_group_ids_train,
      rfx_basis_train
    )
    zero_hat_train <- zero_hat_train + zero_rfx_preds_train
    if (has_rfx_test && has_test) {
      zero_rfx_preds_test <- predictRandomEffectsForLabels(
        zero_rfx_samples,
        group_ids_factor_test,
        rfx_basis_test,
        unseen_groups = rfx_unseen_groups,
        label = "rfx_group_ids_test"
      )
      zero_hat_test <- zero_hat_test + zero_rfx_preds_test
    }
  }

  # Global error variance
  if (sample_sigma2_global) {
    sigma2_global_samples <- global_var_samples * (y_std_train^2)
  }

  # Leaf parameter variance
  if (sample_sigma2_leaf) {
    tau_samples <- leaf_scale_samples
  }

  # Rescale variance forest prediction by global sigma2 (sampled or constant)
  if (include_variance_forest) {
    if (sample_sigma2_global) {
      sigma2_x_hat_train <- sapply(1:num_retained_samples, function(i) {
        sigma2_x_hat_train[, i] * sigma2_global_samples[i]
      })
      if (has_test) {
        sigma2_x_hat_test <- sapply(
          1:num_retained_samples,
          function(i) {
            sigma2_x_hat_test[, i] * sigma2_global_samples[i]
          }
        )
      }
    } else {
      sigma2_x_hat_train <- sigma2_x_hat_train *
        sigma2_init *
        y_std_train *
        y_std_train
      if (has_test) {
        sigma2_x_hat_test <- sigma2_x_hat_test *
          sigma2_init *
          y_std_train *
          y_std_train
      }
    }
  }

  # Return results as a list
  model_params <- list(
    "sigma2_init" = sigma2_init,
    "sigma2_leaf_init" = sigma2_leaf_init,
    "a_global" = a_global,
    "b_global" = b_global,
    "a_leaf" = a_leaf,
    "b_leaf" = b_leaf,
    "a_forest" = a_forest,
    "b_forest" = b_forest,
    "outcome_mean" = y_bar_train,
    "outcome_scale" = y_std_train,
    "standardize" = standardize,
    "leaf_dimension" = leaf_dimension,
    "is_leaf_constant" = is_leaf_constant,
    "leaf_regression" = leaf_regression,
    "requires_basis" = requires_basis,
    "num_covariates" = num_cov_orig,
    "num_basis" = ifelse(
      is.null(leaf_basis_train),
      0,
      ncol(leaf_basis_train)
    ),
    "num_samples" = num_retained_samples,
    "num_gfr" = num_gfr,
    "num_burnin" = num_burnin,
    "num_mcmc" = num_mcmc,
    "keep_every" = keep_every,
    "num_chains" = num_chains,
    "has_basis" = !is.null(leaf_basis_train),
    "has_rfx" = has_rfx,
    "has_rfx_basis" = has_basis_rfx,
    "num_rfx_basis" = num_basis_rfx,
    "sample_sigma2_global" = sample_sigma2_global,
    "sample_sigma2_leaf" = sample_sigma2_leaf,
    "include_mean_forest" = include_mean_forest,
    "include_variance_forest" = include_variance_forest,
    "outcome_model" = outcome_model,
    "probit_outcome_model" = probit_outcome_model,
    "cloglog_num_categories" = ifelse(
      link_is_cloglog,
      cloglog_num_categories,
      0
    ),
    "rfx_model_spec" = rfx_model_spec,
    "rfx_unseen_groups" = rfx_unseen_groups
  )
  if (link_is_log) {
    model_params[["has_offset"]] <- !is.null(offset_train)
    model_params[["negbin_dispersion_prior_shape"]] <- negbin_dispersion_prior_shape
    model_params[["negbin_dispersion_prior_rate"]] <- negbin_dispersion_prior_rate
    model_params[["negbin_dispersion_init"]] <- negbin_dispersion_init
    model_params[["sample_negbin_dispersion"]] <- isTRUE(sample_negbin_dispersion)
  }
  if (outcome_is_zinb) {
    model_params[["zero_outcome_mean"]] <- zero_outcome_mean
    model_params[["zinb_zero_probability_init"]] <- zinb_zero_probability_init
    model_params[["num_trees_zero"]] <- num_trees_zero
    model_params[["sample_sigma2_leaf_zero"]] <- isTRUE(sample_sigma2_leaf_zero)
    model_params[["sigma2_leaf_zero_init"]] <- as.numeric(sigma2_leaf_zero_init)
    model_params[["has_zero_rfx"]] <- has_zero_rfx
  }
  result <- list(
    "model_params" = model_params,
    "train_set_metadata" = X_train_metadata
  )
  if (include_mean_forest) {
    result[["mean_forests"]] = forest_samples_mean
    if (!importance_only) result[["y_hat_train"]] = y_hat_train
    if (has_test) {
      result[["y_hat_test"]] = y_hat_test
    }
    if (link_is_cloglog && !outcome_is_binary) {
      result[["cloglog_cutpoint_samples"]] = cloglog_cutpoint_samples
    }
  }
  if (include_variance_forest) {
    result[["variance_forests"]] = forest_samples_variance
    result[["sigma2_x_hat_train"]] = sigma2_x_hat_train
    if (has_test) result[["sigma2_x_hat_test"]] = sigma2_x_hat_test
  }
  if (sample_sigma2_global) {
    result[["sigma2_global_samples"]] = sigma2_global_samples
  }
  if (sample_sigma2_leaf) {
    result[["sigma2_leaf_samples"]] = tau_samples
  }
  if (link_is_log) {
    result[["negbin_dispersion_samples"]] = negbin_dispersion_samples
  }
  if (outcome_is_zinb) {
    result[["zero_forests"]] = forest_samples_zero
    result[["zero_hat_train"]] = zero_hat_train
    if (has_test) {
      result[["zero_hat_test"]] = zero_hat_test
    }
    if (sample_sigma2_leaf_zero) {
      result[["sigma2_leaf_zero_samples"]] = leaf_scale_zero_samples
    }
    if (has_zero_rfx) {
      result[["zero_rfx_samples"]] = zero_rfx_samples
      result[["zero_rfx_preds_train"]] = zero_rfx_preds_train
      if (has_rfx_test && has_test) {
        result[["zero_rfx_preds_test"]] = zero_rfx_preds_test
      }
    }
  }
  if (has_rfx) {
    result[["rfx_samples"]] = rfx_samples
    if (!importance_only) result[["rfx_preds_train"]] = rfx_preds_train
    result[["rfx_unique_group_ids"]] = levels(group_ids_factor)
  }
  if ((has_rfx_test) && (has_test)) {
    result[["rfx_preds_test"]] = rfx_preds_test
  }
  if (importance_only) {
    result$mean_forests <- NULL
    if (any(importance_split_counts > .Machine$integer.max)) {
      stop("Importance split counts exceed the supported integer range.")
    }
    result$split_counts <- as.integer(importance_split_counts)
    class(result) <- "bartimportance"
  } else {
    class(result) <- "bartmodel"
  }

  if (selected_only) {
    compact_model <- result
    keep <- retained_sample_ids + 1L
    compact_model$y_hat_train <- as.matrix(result$y_hat_train)[, keep, drop = FALSE]
    compact_model$y_hat_test <- NULL
    if (sample_sigma2_global) compact_model$sigma2_global_samples <- result$sigma2_global_samples[keep]
    if (sample_sigma2_leaf) compact_model$sigma2_leaf_samples <- result$sigma2_leaf_samples[keep]
    if (has_rfx) {
      # The sample container is shared with `result`, whose train and test
      # predictions are already materialized as matrices above, so the retained
      # draws are all the compact model needs.
      compact_model$rfx_preds_train <- as.matrix(result$rfx_preds_train)[, keep, drop = FALSE]
      compact_model$rfx_samples$retain_samples(retained_sample_ids)
    }
    compact_model$model_params$num_samples <- length(keep)
    result$mean_forests <- NULL
    result$model <- compact_model
    result$retained_sample_ids <- retained_sample_ids
    class(result) <- "bartselected"
  }

  # Clean up classes with external pointers to C++ data structures
  if (include_mean_forest) {
    rm(forest_model_mean)
  }
  if (include_variance_forest) {
    rm(forest_model_variance)
  }
  if (outcome_is_zinb) {
    rm(forest_model_zero, forest_dataset_zero, outcome_zero)
  }
  rm(forest_dataset_train)
  if (has_test) {
    rm(forest_dataset_test)
  }
  if (has_rfx) {
    rm(rfx_dataset_train, rfx_tracker_train, rfx_model)
  }
  if (has_zero_rfx) {
    rm(zero_rfx_dataset_train, zero_rfx_tracker_train, zero_rfx_model)
  }
  rm(outcome_train)
  rm(rng)

  # Restore global RNG state if user provided a random seed
  if (custom_rng) {
    if (has_existing_random_seed) {
      .Random.seed <- original_global_seed
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  return(result)
}

#' @title Predict from a BART Model
#' @description
#' Predict from a sampled BART model on new data
#'
#' @param object Object of type `bart` containing draws of a regression forest and associated sampling outputs.
#' @param X Covariates used to determine tree leaf predictions for each observation. Must be passed as a matrix or dataframe.
#' @param leaf_basis (Optional) Bases used for prediction (by e.g. dot product with leaf values). Default: `NULL`.
#' @param rfx_group_ids (Optional) Test set group labels used for an additive random effects model.
#' We do not currently support (but plan to in the near future), test set evaluation for group labels
#' that were not in the training set.
#' @param rfx_basis (Optional) Test set basis for "random-slope" regression in additive random effects model.
#' @param type (Optional) Type of prediction to return. Options are "mean", which averages the predictions from every draw of a BART model, and "posterior", which returns the entire matrix of posterior predictions. Default: "posterior".
#' @param terms (Optional) Which model terms to include in the prediction. This can be a single term or a list of model terms. Options include "y_hat", "mean_forest", "rfx", "variance_forest", "zero_forest" and "zero_rfx" (zero-inflated count models), or "all". If a model doesn't have mean forest, random effects, or variance forest predictions, but one of those terms is request, the request will simply be ignored. If none of the requested terms are present in a model, this function will return `NULL` along with a warning. Default: "all".
#' @param scale (Optional) Scale of mean function predictions. Options are "linear", which returns predictions on the original scale of the mean forest / RFX terms, "probability", which transforms predictions into class probabilities for models with discrete outcomes, and "class", which returns predicted outcome categories for discrete outcome models. "probability" is only valid for outcome models with `outcome == 'binary'` or `outcome == 'ordinal'`. For binary outcomes, this will return the probability that `y == 1`, and for ordinal outcomes, this will return probabilities for each outcome label. For count (negative-binomial) models, "linear" returns the log rate `eta = offset + c + f(X) [+ rfx]`, "rate" returns `mu = exp(eta)`, and "response" is an alias of "rate"; on the "rate" scale the mean-forest term is `exp(c + f(X))` and the random-effects term is the multiplicative rate ratio `exp(rfx)`. "rate" and "response" are only valid for count models, and "probability" and "class" are not valid for them. For zero-inflated count models (`family = "zinb"`), "linear" and "rate" describe the count part as above, "response" returns the mean `(1 - pi) * mu`, and "zero_probability" returns the structural-zero probability `pi = Phi(c0 + g(X) [+ zero-part rfx])` (only the "y_hat", "zero_forest" and "zero_rfx" terms are defined on that scale, and "all" is narrowed to them); the "zero_forest" term is the forest's probit index `c0 + g(X)` on the "linear" scale and `Phi(c0 + g(X))` on the other scales, and the "zero_rfx" term (models with zero-part random effects) is the zero-part random effect on the probit scale on every scale. With zero-part random effects, `rfx_group_ids` are needed for "y_hat" on the "response" and "zero_probability" scales. With `type = "mean"`, transformed draws are averaged, so "rate" gives the posterior mean of `mu`. Default: "linear".
#' @param offset (Optional) Count models only. Log-scale offset of length `nrow(X)` added to the linear predictor, as `offset_train` in [bart()]. If the model was fit with an offset and none is given here, an offset of 0 (unit exposure) is used with a warning. Default: `NULL`.
#' @param ... (Optional) Other prediction parameters.
#'
#' @return List of prediction matrices or single prediction matrix / vector, depending on the terms requested.
#' @export
#'
#' @examples
#' n <- 100
#' p <- 5
#' X <- matrix(runif(n*p), ncol = p)
#' f_XW <- (
#'     ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) +
#'     ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) +
#'     ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) +
#'     ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
#' )
#' noise_sd <- 1
#' y <- f_XW + rnorm(n, 0, noise_sd)
#' test_set_pct <- 0.2
#' n_test <- round(test_set_pct*n)
#' n_train <- n - n_test
#' test_inds <- sort(sample(1:n, n_test, replace = FALSE))
#' train_inds <- (1:n)[!((1:n) %in% test_inds)]
#' X_test <- X[test_inds,]
#' X_train <- X[train_inds,]
#' y_test <- y[test_inds]
#' y_train <- y[train_inds]
#' bart_model <- bart(X_train = X_train, y_train = y_train,
#'                    num_gfr = 10, num_burnin = 0, num_mcmc = 10)
#' y_hat_test <- predict(bart_model, X=X_test)$y_hat
predict.bartmodel <- function(
  object,
  X,
  leaf_basis = NULL,
  rfx_group_ids = NULL,
  rfx_basis = NULL,
  type = "posterior",
  terms = "all",
  scale = "linear",
  offset = NULL,
  ...
) {
  # Handle mean function scale
  if (!is.character(scale)) {
    stop("scale must be a string or character vector")
  }
  if (
    !(scale %in%
      c("linear", "probability", "class", "rate", "response", "zero_probability"))
  ) {
    stop(
      "scale must be 'linear', 'probability', or 'class' ('rate' or 'response' for count models, 'zero_probability' for zero-inflated count models)"
    )
  }
  outcome_model <- object$model_params$outcome_model
  is_count <- identical(outcome_model$outcome, "count") &&
    identical(outcome_model$link, "log")
  is_zinb <- is_count && identical(outcome_model$family, "zinb")
  if (!is_zinb && scale == "zero_probability") {
    stop(
      "scale can only be 'zero_probability' for a zero-inflated count model (family = 'zinb')"
    )
  }
  zero_probability_scale <- scale == "zero_probability"
  if (is_count && (scale %in% c("probability", "class"))) {
    stop(
      "scale cannot be 'probability' or 'class' for a count model; use 'linear', 'rate' or 'response'"
    )
  }
  if (!is_count && (scale %in% c("rate", "response"))) {
    stop("scale can only be 'rate' or 'response' for a count (negative-binomial) model")
  }
  if (!is_count && !is.null(offset)) {
    stop("offset is only supported for count (negative-binomial) models")
  }
  rate_scale <- is_count && (scale %in% c("rate", "response"))
  is_probit <- (outcome_model$link == "probit" &&
    outcome_model$outcome == "binary")
  is_binary_cloglog <- (outcome_model$link == "cloglog" &&
    outcome_model$outcome == "binary")
  is_ordinal_cloglog <- (outcome_model$link == "cloglog" &&
    outcome_model$outcome == "ordinal")
  is_cloglog <- is_binary_cloglog || is_ordinal_cloglog
  if ((scale == "probability") && (!(is_probit || is_cloglog))) {
    stop(
      "scale cannot be 'probability' for models not fit with a probit or cloglog outcome model"
    )
  }
  if ((scale == "class") && (!(is_probit || is_cloglog))) {
    stop(
      "scale cannot be 'class' for models not fit with a probit or cloglog outcome model"
    )
  }
  probability_scale <- scale == "probability"
  class_scale <- scale == "class"

  # Handle prediction type
  if (!is.character(type)) {
    stop("type must be a string or character vector")
  }
  if (!(type %in% c("mean", "posterior"))) {
    stop("type must either be 'mean' or 'posterior'")
  }
  predict_mean <- type == "mean"
  if (predict_mean && class_scale) {
    stop("Posterior mean predictions are not supported for scale = 'class'")
  }

  # Handle prediction terms
  rfx_model_spec <- object$model_params$rfx_model_spec
  rfx_intercept <- rfx_model_spec == "intercept_only"
  if (!is.character(terms)) {
    stop("type must be a string or character vector")
  }
  num_terms <- length(terms)
  has_mean_forest <- object$model_params$include_mean_forest
  has_variance_forest <- object$model_params$include_variance_forest
  has_rfx <- object$model_params$has_rfx
  has_zero_forest <- is_zinb
  has_zero_rfx <- is_zinb && isTRUE(object$model_params$has_zero_rfx)
  has_y_hat <- has_mean_forest || has_rfx
  predict_y_hat <- (((has_y_hat) && ("y_hat" %in% terms)) ||
    ((has_y_hat) && ("all" %in% terms)))
  predict_mean_forest <- (((has_mean_forest) && ("mean_forest" %in% terms)) ||
    ((has_mean_forest) && ("all" %in% terms)))
  predict_rfx <- (((has_rfx) && ("rfx" %in% terms)) ||
    ((has_rfx) && ("all" %in% terms)))
  predict_variance_forest <- (((has_variance_forest) &&
    ("variance_forest" %in% terms)) ||
    ((has_variance_forest) && ("all" %in% terms)))
  predict_zero_forest <- has_zero_forest &&
    (("zero_forest" %in% terms) || ("all" %in% terms))
  predict_zero_rfx <- has_zero_rfx &&
    (("zero_rfx" %in% terms) || ("all" %in% terms))
  if (zero_probability_scale) {
    # On the zero-probability scale only y_hat (= pi), the zero forest and the
    # zero-part random effects are defined; "all" is narrowed to those terms
    if (any(c("mean_forest", "rfx") %in% terms)) {
      stop(
        "scale = 'zero_probability' is only defined for the 'y_hat', 'zero_forest' and 'zero_rfx' terms"
      )
    }
    predict_mean_forest <- FALSE
    predict_rfx <- FALSE
  }
  predict_count <- sum(c(
    predict_y_hat,
    predict_mean_forest,
    predict_rfx,
    predict_variance_forest,
    predict_zero_forest,
    predict_zero_rfx
  ))
  if (predict_count == 0) {
    warning(paste0(
      "None of the requested model terms, ",
      paste(terms, collapse = ", "),
      ", were fit in this model"
    ))
    return(NULL)
  }
  predict_rfx_intermediate <- (predict_y_hat && has_rfx &&
    !zero_probability_scale)
  predict_mean_forest_intermediate <- (predict_y_hat && has_mean_forest &&
    !zero_probability_scale)
  # The zero part enters y_hat on the "response" and "zero_probability" scales
  predict_zero_forest_intermediate <- predict_y_hat && is_zinb &&
    (scale %in% c("response", "zero_probability"))
  predict_zero_rfx_intermediate <- predict_zero_forest_intermediate &&
    has_zero_rfx
  if (class_scale) {
    if (!((predict_count == 1) && (predict_y_hat))) {
      stop("Class scale can only be used with y_hat predictions")
    }
  }

  # Check that we have at least one term to predict on probability scale
  if (
    probability_scale &&
      !predict_y_hat &&
      !predict_mean_forest &&
      !predict_rfx
  ) {
    stop(
      "scale can only be 'probability' if at least one mean term is requested"
    )
  }

  # Check that covariates are matrix or data frame
  if ((!is.data.frame(X)) && (!is.matrix(X))) {
    stop("X must be a matrix or dataframe")
  }

  # Convert all input data to matrices if not already converted
  if ((is.null(dim(leaf_basis))) && (!is.null(leaf_basis))) {
    leaf_basis <- as.matrix(leaf_basis)
  }
  if ((is.null(dim(rfx_basis))) && (!is.null(rfx_basis))) {
    if (predict_rfx || predict_zero_rfx || predict_zero_rfx_intermediate) {
      rfx_basis <- as.matrix(rfx_basis)
    }
  }

  # Data checks
  if ((object$model_params$requires_basis) && (is.null(leaf_basis))) {
    stop("Basis (leaf_basis) must be provided for this model")
  }
  if ((!is.null(leaf_basis)) && (nrow(X) != nrow(leaf_basis))) {
    stop("X and leaf_basis must have the same number of rows")
  }
  if (object$model_params$num_covariates != ncol(X)) {
    stop(
      "X must contain the same number of columns as the BART model's training dataset"
    )
  }
  if ((predict_rfx) && (is.null(rfx_group_ids))) {
    stop(
      "Random effect group labels (rfx_group_ids) must be provided for this model"
    )
  }
  if (
    (predict_zero_rfx || predict_zero_rfx_intermediate) &&
      is.null(rfx_group_ids)
  ) {
    stop(
      "Random effect group labels (rfx_group_ids) must be provided for this model: its zero part has random effects"
    )
  }
  if (
    (predict_zero_rfx || predict_zero_rfx_intermediate) &&
      is.null(rfx_basis) &&
      !rfx_intercept
  ) {
    stop("Random effects basis (rfx_basis) must be provided for this model")
  }
  if ((predict_rfx) && (is.null(rfx_basis)) && (!rfx_intercept)) {
    stop("Random effects basis (rfx_basis) must be provided for this model")
  }
  if ((object$model_params$num_rfx_basis > 0) && (!rfx_intercept)) {
    if (ncol(rfx_basis) != object$model_params$num_rfx_basis) {
      stop(
        "Random effects basis has a different dimension than the basis used to train this model"
      )
    }
  }

  # Count models: offset on the log scale, added to y_hat only
  if (is_count) {
    if (is.null(offset)) {
      if (
        isTRUE(object$model_params$has_offset) &&
          predict_y_hat &&
          !zero_probability_scale
      ) {
        warning(
          "This model was fit with an offset but none was supplied; predictions use an offset of 0 (unit exposure)"
        )
      }
      offset <- rep(0, nrow(X))
    } else {
      offset <- .validateNegbinOffset(offset, nrow(X), "offset")
    }
  }

  # Preprocess covariates
  train_set_metadata <- object$train_set_metadata
  X <- preprocessPredictionData(X, train_set_metadata)

  # Unseen group labels either stop the prediction or take the posterior mean
  # random effect across the sampled groups, depending on how the model was
  # sampled. Models serialized before this option default to "error".
  rfx_unseen_groups <- validateRfxUnseenGroups(
    object$model_params[["rfx_unseen_groups"]]
  )

  # Recode group IDs to integer vector (if passed as, for example, a vector of county names, etc...)
  group_ids_factor <- NULL
  if (!is.null(rfx_group_ids)) {
    rfx_unique_group_ids <- object$rfx_unique_group_ids
    group_ids_factor <- factor(rfx_group_ids, levels = rfx_unique_group_ids)
    if (
      sum(is.na(group_ids_factor)) > 0 &&
        identical(rfx_unseen_groups, "error")
    ) {
      stop(
        "All random effect group labels provided in rfx_group_ids must have been present in rfx_group_ids_train"
      )
    }
    rfx_group_ids <- as.integer(group_ids_factor)
  }

  # Handle RFX model specification
  if (has_rfx) {
    if (object$model_params$rfx_model_spec == "custom") {
      if (is.null(rfx_basis)) {
        stop(
          "A user-provided basis (`rfx_basis`) must be provided when the model was sampled with a random effects model spec set to 'custom'"
        )
      }
    } else if (object$model_params$rfx_model_spec == "intercept_only") {
      # Only construct a basis if user-provided basis missing
      if (is.null(rfx_basis)) {
        rfx_basis <- matrix(
          rep(1, nrow(X)),
          nrow = nrow(X),
          ncol = 1
        )
      }
    }
  }

  # Create prediction dataset
  if (!is.null(leaf_basis)) {
    prediction_dataset <- createForestDataset(X, leaf_basis)
  } else {
    prediction_dataset <- createForestDataset(X)
  }

  # Compute variance forest predictions
  if (predict_variance_forest) {
    s_x_raw <- object$variance_forests$predict(prediction_dataset)
  }

  # Scale variance forest predictions
  num_samples <- object$model_params$num_samples
  y_std <- object$model_params$outcome_scale
  y_bar <- object$model_params$outcome_mean
  sigma2_init <- object$model_params$sigma2_init
  if (predict_variance_forest) {
    if (object$model_params$sample_sigma2_global) {
      sigma2_global_samples <- object$sigma2_global_samples
      variance_forest_predictions <- sapply(1:num_samples, function(i) {
        s_x_raw[, i] * sigma2_global_samples[i]
      })
    } else {
      variance_forest_predictions <- s_x_raw * sigma2_init * y_std * y_std
    }
    if (predict_mean) {
      variance_forest_predictions <- rowMeans(variance_forest_predictions)
    }
  }

  # Compute mean forest predictions
  if (predict_mean_forest || predict_mean_forest_intermediate) {
    mean_forest_predictions <- object$mean_forests$predict(
      prediction_dataset
    ) *
      y_std +
      y_bar
  }

  # Compute zero forest predictions (probit index of the zero probability)
  if (predict_zero_forest || predict_zero_forest_intermediate) {
    zero_forest_linear <- object$zero_forests$predict(prediction_dataset) +
      object$model_params$zero_outcome_mean
    if (predict_zero_forest) {
      zero_forest_predictions <- if (scale == "linear") {
        zero_forest_linear
      } else {
        pnorm(zero_forest_linear)
      }
    }
  }

  # Zero-part random effects, on the probit scale on every prediction scale;
  # they shift the probit index of pi (unseen groups follow the model's
  # unseen-group policy, as the count part does)
  if (predict_zero_rfx || predict_zero_rfx_intermediate) {
    zero_rfx_predictions <- predictRandomEffectsForLabels(
      object$zero_rfx_samples,
      group_ids_factor,
      rfx_basis,
      unseen_groups = rfx_unseen_groups
    )
    if (predict_zero_rfx_intermediate) {
      zero_forest_linear <- zero_forest_linear + zero_rfx_predictions
    }
  }

  # Compute rfx predictions (if needed)
  if (predict_rfx || predict_rfx_intermediate) {
    if (!is.null(rfx_basis)) {
      rfx_predictions <- predictRandomEffectsForLabels(
        object$rfx_samples,
        group_ids_factor,
        rfx_basis,
        unseen_groups = rfx_unseen_groups
      ) *
        y_std
    } else {
      # Sanity check -- this branch should only occur if rfx_model_spec == "intercept_only"
      if (!rfx_intercept) {
        stop(
          "rfx_basis must be provided for random effects models with random slopes"
        )
      }

      # Extract the raw RFX samples and scale by train set outcome standard deviation
      rfx_param_list <- object$rfx_samples$extract_parameter_samples()
      rfx_beta_draws <- rfx_param_list$beta_samples * y_std

      # Promote to an array with consistent dimensions when there's one rfx term
      if (length(dim(rfx_beta_draws)) == 2) {
        dim(rfx_beta_draws) <- c(1, dim(rfx_beta_draws))
      }

      # Construct a matrix with the appropriate group random effects arranged for each observation
      rfx_predictions_raw <- array(
        NA,
        dim = c(
          nrow(X),
          ncol(rfx_basis),
          object$model_params$num_samples
        )
      )
      for (i in 1:nrow(X)) {
        rfx_predictions_raw[i, , ] <-
          rfx_beta_draws[, rfx_group_ids[i], ]
      }

      # Intercept-only model, so the random effect prediction is simply the
      # value of the respective group's intercept coefficient for each observation
      rfx_predictions = rfx_predictions_raw[, 1, ]
    }
  }

  # Combine into y hat predictions
  if (probability_scale || class_scale) {
    if (is_probit) {
      if (predict_y_hat) {
        if (has_mean_forest && has_rfx) {
          y_hat <- pnorm(mean_forest_predictions + rfx_predictions)
          mean_forest_predictions <- pnorm(mean_forest_predictions)
          rfx_predictions <- pnorm(rfx_predictions)
        } else if (has_mean_forest) {
          y_hat <- pnorm(mean_forest_predictions)
          mean_forest_predictions <- pnorm(mean_forest_predictions)
        } else if (has_rfx) {
          y_hat <- pnorm(rfx_predictions)
          rfx_predictions <- pnorm(rfx_predictions)
        }
      } else {
        if (has_mean_forest && has_rfx) {
          mean_forest_predictions <- pnorm(mean_forest_predictions)
          rfx_predictions <- pnorm(rfx_predictions)
        } else if (has_mean_forest) {
          mean_forest_predictions <- pnorm(mean_forest_predictions)
        } else if (has_rfx) {
          rfx_predictions <- pnorm(rfx_predictions)
        }
      }
    } else if (is_binary_cloglog) {
      mean_forest_predictions <- exp(-exp(mean_forest_predictions))
      if (predict_y_hat) {
        y_hat <- mean_forest_predictions
      }
    } else if (is_ordinal_cloglog) {
      cloglog_num_categories <- object$model_params$cloglog_num_categories
      cloglog_cutpoint_samples <- object$cloglog_cutpoint_samples
      n_obs_pred <- nrow(X)
      n_samp_pred <- object$model_params$num_samples
      mean_forest_probabilities <- array(
        NA_real_,
        dim = c(n_obs_pred, cloglog_num_categories, n_samp_pred)
      )
      # Sequential ordinal cloglog: P(Y=k) = prod_{j<k} S_j * (1 - S_k)
      # S_k = exp(-exp(gamma_k + f)), running survival product across k.
      survival_product <- matrix(1.0, nrow = n_obs_pred, ncol = n_samp_pred)
      for (k in seq_len(cloglog_num_categories - 1)) {
        S_k <- exp(
          -exp(sweep(
            mean_forest_predictions,
            2,
            cloglog_cutpoint_samples[k, ],
            "+"
          ))
        )
        mean_forest_probabilities[, k, ] <- survival_product * (1 - S_k)
        survival_product <- survival_product * S_k
      }
      mean_forest_probabilities[, cloglog_num_categories, ] <- survival_product
      if (predict_y_hat) {
        y_hat <- mean_forest_probabilities
      }
      mean_forest_predictions <- mean_forest_probabilities
    }
  } else if (zero_probability_scale) {
    if (predict_y_hat) {
      y_hat <- pnorm(zero_forest_linear)
    }
  } else {
    if (predict_y_hat && has_mean_forest && has_rfx) {
      y_hat <- mean_forest_predictions + rfx_predictions
    } else if (predict_y_hat && has_mean_forest) {
      y_hat <- mean_forest_predictions
    } else if (predict_y_hat && has_rfx) {
      y_hat <- rfx_predictions
    }
    if (is_count) {
      if (predict_y_hat) {
        y_hat <- y_hat + offset
      }
      if (rate_scale) {
        if (predict_y_hat) {
          y_hat <- exp(y_hat)
        }
        if (predict_mean_forest) {
          mean_forest_predictions <- exp(mean_forest_predictions)
        }
        if (predict_rfx) {
          rfx_predictions <- exp(rfx_predictions)
        }
      }
      # Zero-inflated mean E[y | x] = (1 - pi) mu
      if (is_zinb && scale == "response" && predict_y_hat) {
        y_hat <- (1 - pnorm(zero_forest_linear)) * y_hat
      }
    }
  }

  # Collapse to posterior mean predictions if requested
  if (predict_mean) {
    if (predict_mean_forest) {
      if (is_ordinal_cloglog && probability_scale) {
        mean_forest_predictions <- apply(mean_forest_predictions, c(1, 2), mean)
      } else {
        mean_forest_predictions <- rowMeans(mean_forest_predictions)
      }
    }
    if (predict_rfx) {
      # Note: random effects not supported for cloglog, so we don't have to handle the ordinal / cloglog case here
      rfx_predictions <- rowMeans(rfx_predictions)
    }
    if (predict_zero_forest) {
      zero_forest_predictions <- rowMeans(zero_forest_predictions)
    }
    if (predict_zero_rfx) {
      zero_rfx_predictions <- rowMeans(zero_rfx_predictions)
    }
    if (predict_y_hat) {
      if (is_ordinal_cloglog && probability_scale) {
        y_hat <- apply(y_hat, c(1, 2), mean)
      } else {
        y_hat <- rowMeans(y_hat)
      }
    }
  }

  # Convert probabilities to classes if requested
  if (class_scale) {
    if (is_ordinal_cloglog) {
      y_hat <- apply(y_hat, c(1, 3), which.max)
    } else {
      y_hat <- ifelse(y_hat < 0.5, 0, 1)
    }
  }

  if (predict_count == 1) {
    if (predict_y_hat) {
      return(y_hat)
    } else if (predict_mean_forest) {
      return(mean_forest_predictions)
    } else if (predict_rfx) {
      return(rfx_predictions)
    } else if (predict_variance_forest) {
      return(variance_forest_predictions)
    } else if (predict_zero_forest) {
      return(zero_forest_predictions)
    } else if (predict_zero_rfx) {
      return(zero_rfx_predictions)
    }
  } else {
    result <- list()
    if (predict_y_hat) {
      result[["y_hat"]] = y_hat
    } else {
      result[["y_hat"]] <- NULL
    }
    if (predict_mean_forest) {
      result[["mean_forest_predictions"]] = mean_forest_predictions
    } else {
      result[["mean_forest_predictions"]] <- NULL
    }
    if (predict_rfx) {
      result[["rfx_predictions"]] = rfx_predictions
    } else {
      result[["rfx_predictions"]] <- NULL
    }
    if (predict_variance_forest) {
      result[["variance_forest_predictions"]] = variance_forest_predictions
    } else {
      result[["variance_forest_predictions"]] <- NULL
    }
    if (predict_zero_forest) {
      result[["zero_forest_predictions"]] = zero_forest_predictions
    }
    if (predict_zero_rfx) {
      result[["zero_rfx_predictions"]] = zero_rfx_predictions
    }
    return(result)
  }
}

#' @title Print Summary of BART Model
#' @description Prints a summary of the BART model, including the model terms and their specifications.
#' @param x The BART model object
#' @param ... Additional arguments
#' @export
#' @return BART model object unchanged after printing summary
print.bartmodel <- function(x, ...) {
  # What type of model was run
  model_terms <- c()
  if (x$model_params$include_mean_forest) {
    model_terms <- c(model_terms, "mean forest")
  }
  if (x$model_params$include_variance_forest) {
    model_terms <- c(model_terms, "variance forest")
  }
  if (x$model_params$has_rfx) {
    model_terms <- c(model_terms, "additive random effects")
  }
  if (x$model_params$sample_sigma2_global) {
    model_terms <- c(model_terms, "global error variance model")
  }
  if (x$model_params$sample_sigma2_leaf) {
    model_terms <- c(model_terms, "mean forest leaf scale model")
  }
  if (identical(x$model_params$outcome_model$family, "zinb")) {
    model_terms <- c(model_terms, "zero forest")
    if (isTRUE(x$model_params$has_zero_rfx)) {
      model_terms <- c(model_terms, "zero-part random effects")
    }
  }
  if (isTRUE(x$model_params$sample_negbin_dispersion)) {
    model_terms <- c(model_terms, "negative-binomial dispersion model")
  }
  if (length(model_terms) > 2) {
    summary_message <- paste0(
      "nashtree::bart() run with ",
      paste0(
        paste0(model_terms[1:(length(model_terms) - 1)], collapse = ", "),
        ", and ",
        model_terms[length(model_terms)]
      )
    )
  } else if (length(model_terms) == 2) {
    summary_message <- paste0(
      "nashtree::bart() run with ",
      paste0(model_terms, collapse = " and ")
    )
  } else {
    summary_message <- paste0("nashtree::bart() run with ", model_terms)
  }

  # Outcome and leaf model details
  outcome_model <- x$model_params$outcome_model
  is_probit <- (outcome_model$link == "probit" &&
    outcome_model$outcome == "binary")
  is_binary_cloglog <- (outcome_model$link == "cloglog" &&
    outcome_model$outcome == "binary")
  is_ordinal_cloglog <- (outcome_model$link == "cloglog" &&
    outcome_model$outcome == "ordinal")
  if (is_ordinal_cloglog) {
    num_categories <- x$model_params$cloglog_num_categories
    outcome_model_summary <- paste0(
      "Ordinal outcome with ",
      num_categories,
      " categories was modeled with a complementary log-log (cloglog) link function"
    )
  } else if (is_binary_cloglog) {
    outcome_model_summary <- paste0(
      "Binary outcome was modeled with a complementary log-log (cloglog) link function"
    )
  } else if (is_probit) {
    outcome_model_summary <- paste0(
      "Binary outcome was modeled with a probit link function"
    )
  } else if (identical(outcome_model$family, "zinb")) {
    outcome_model_summary <- paste0(
      "Count outcome was modeled as zero-inflated negative binomial: a log-link count part (Polya-Gamma augmentation)",
      ifelse(isTRUE(x$model_params$has_offset), " with an offset", ""),
      " and a probit zero part (",
      x$model_params$num_trees_zero,
      " trees",
      ifelse(isTRUE(x$model_params$has_zero_rfx), " and random effects", ""),
      ")"
    )
  } else if (identical(outcome_model$outcome, "count")) {
    outcome_model_summary <- paste0(
      "Count outcome was modeled as negative binomial with a log link (Polya-Gamma augmentation)",
      ifelse(isTRUE(x$model_params$has_offset), " and an offset", "")
    )
  } else {
    outcome_model_summary <- paste0(
      "Continuous outcome was modeled as Gaussian"
    )
  }
  if (x$model_params$leaf_regression) {
    summary_message <- paste0(
      summary_message,
      "\n",
      outcome_model_summary,
      " with a leaf regression prior with ",
      x$model_params$leaf_dimension,
      " bases for the mean forest"
    )
  } else if (x$model_params$include_mean_forest) {
    summary_message <- paste0(
      summary_message,
      "\n",
      outcome_model_summary,
      " with a constant leaf prior for the mean forest"
    )
  } else {
    summary_message <- paste0(
      summary_message,
      "\n",
      outcome_model_summary,
    )
  }

  # Standardization
  if (x$model_params$standardize) {
    summary_message <- paste0(
      summary_message,
      "\n",
      "Outcome was standardized"
    )
  }

  # Random effects details
  if (x$model_params$has_rfx) {
    if (x$model_params$rfx_model_spec == "custom") {
      summary_message <- paste0(
        summary_message,
        "\n",
        "Random effects were fit with a user-supplied basis"
      )
    } else if (x$model_params$rfx_model_spec == "intercept_only") {
      summary_message <- paste0(
        summary_message,
        "\n",
        "Random effects were fit with an 'intercept-only' parameterization"
      )
    }
  }

  # Sampler details
  summary_message <- paste0(
    summary_message,
    "\n",
    "The sampler was run for ",
    x$model_params$num_gfr,
    " GFR iterations, with ",
    x$model_params$num_chains,
    ifelse(
      x$model_params$num_chains == 1,
      " chain of ",
      " chains of "
    ),
    x$model_params$num_burnin,
    " burn-in iterations and ",
    x$model_params$num_mcmc,
    " MCMC iterations, ",
    ifelse(
      x$model_params$keep_every == 1,
      "retaining every iteration (i.e. no thinning)",
      paste0(
        "retaining every ",
        x$model_params$keep_every,
        "th iteration (i.e. thinning)"
      )
    )
  )

  # Print the model details
  cat(summary_message, "\n")

  # Return bart_model invisibly
  invisible(x)
}

#' @title Summarize BART Model Fit and Parameters
#' @description Summarize a BART fit with a description of the model that was fit and numeric summaries of any sampled quantities.
#' @param object The BART model object
#' @param ... Additional arguments
#' @export
#' @return BART model object unchanged after summarizing
summary.bartmodel <- function(object, ...) {
  # First, print the BART model
  tmp <- print(object)

  # Summarize any sampled quantities

  # Global error scale
  if (object$model_params$sample_sigma2_global) {
    sigma2_samples <- object$sigma2_global_samples
    n_samples <- length(sigma2_samples)
    mean_sigma2 <- mean(sigma2_samples)
    sd_sigma2 <- sd(sigma2_samples)
    quantiles_sigma2 <- quantile(
      sigma2_samples,
      probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
    )
    cat(sprintf(
      "Summary of sigma^2 posterior: \n%d samples, mean = %.3f, standard deviation = %.3f, quantiles:\n",
      n_samples,
      mean_sigma2,
      sd_sigma2
    ))
    print(quantiles_sigma2)
  }

  # Leaf scale
  if (object$model_params$sample_sigma2_leaf) {
    sigma2_leaf_samples <- object$sigma2_leaf_samples
    n_samples <- length(sigma2_leaf_samples)
    mean_sigma2 <- mean(sigma2_leaf_samples)
    sd_sigma2 <- sd(sigma2_leaf_samples)
    quantiles_sigma2 <- quantile(
      sigma2_leaf_samples,
      probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
    )
    cat(sprintf(
      "Summary of leaf scale posterior: \n%d samples, mean = %.3f, standard deviation = %.3f, quantiles:\n",
      n_samples,
      mean_sigma2,
      sd_sigma2
    ))
    print(quantiles_sigma2)
  }

  # Negative-binomial dispersion
  if (!is.null(object$negbin_dispersion_samples)) {
    r_samples <- object$negbin_dispersion_samples
    cat(sprintf(
      "Summary of negative-binomial dispersion (size r) posterior: \n%d samples, mean = %.3f, standard deviation = %.3f, quantiles:\n",
      length(r_samples),
      mean(r_samples),
      sd(r_samples)
    ))
    print(quantile(r_samples, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
  }

  # Zero-inflated family: structural-zero probability
  if (!is.null(object$zero_hat_train)) {
    zero_probability_mean <- rowMeans(pnorm(as.matrix(object$zero_hat_train)))
    cat(sprintf(
      "Summary of in-sample structural-zero probability (pi) posterior means: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n",
      length(zero_probability_mean),
      mean(zero_probability_mean),
      sd(zero_probability_mean)
    ))
    print(quantile(
      zero_probability_mean,
      probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
    ))
  }

  # Determine whether outcome model is binary / ordinal
  outcome_model <- object$model_params$outcome_model
  is_probit <- (outcome_model$link == "probit" &&
    outcome_model$outcome == "binary")
  is_binary_cloglog <- (outcome_model$link == "cloglog" &&
    outcome_model$outcome == "binary")
  is_ordinal_cloglog <- (outcome_model$link == "cloglog" &&
    outcome_model$outcome == "ordinal")
  non_continuous_outcome <- (is_probit ||
    is_binary_cloglog ||
    is_ordinal_cloglog)

  # In-sample predictions
  if (!is.null(object$y_hat_train)) {
    y_hat_train_mean <- rowMeans(object$y_hat_train)
    n_y_hat_train <- length(y_hat_train_mean)
    mean_y_hat_train <- mean(y_hat_train_mean)
    sd_y_hat_train <- sd(y_hat_train_mean)
    quantiles_y_hat_train <- quantile(
      y_hat_train_mean,
      probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
    )
    if (identical(outcome_model$outcome, "count")) {
      summary_text <- "Summary of in-sample log-rate (linear predictor) posterior means: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n"
    } else if (non_continuous_outcome) {
      summary_text <- "Summary of in-sample inverse-link-scale posterior predictions: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n"
    } else {
      summary_text <- "Summary of in-sample posterior mean predictions: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n"
    }
    cat(sprintf(
      summary_text,
      n_y_hat_train,
      mean_y_hat_train,
      sd_y_hat_train
    ))
    print(quantiles_y_hat_train)
  }

  # Test-set predictions
  if (!is.null(object$y_hat_test)) {
    y_hat_test_mean <- rowMeans(object$y_hat_test)
    n_y_hat_test <- length(y_hat_test_mean)
    mean_y_hat_test <- mean(y_hat_test_mean)
    sd_y_hat_test <- sd(y_hat_test_mean)
    quantiles_y_hat_test <- quantile(
      y_hat_test_mean,
      probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
    )
    if (identical(outcome_model$outcome, "count")) {
      summary_text <- "Summary of test-set log-rate (linear predictor) posterior means: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n"
    } else if (non_continuous_outcome) {
      summary_text <- "Summary of test-set inverse-link-scale posterior predictions: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n"
    } else {
      summary_text <- "Summary of test-set posterior mean predictions: \n%d observations, mean = %.3f, standard deviation = %.3f, quantiles:\n"
    }
    cat(sprintf(
      summary_text,
      n_y_hat_test,
      mean_y_hat_test,
      sd_y_hat_test
    ))
    print(quantiles_y_hat_test)
  }

  # Random effects
  if (object$model_params$has_rfx) {
    rfx_samples <- getRandomEffectSamples(object)
    rfx_beta_samples <- rfx_samples$beta_samples
    if (length(dim(rfx_beta_samples)) > 2) {
      cat(
        "Random effects summary of variance components across groups and posterior draws:\n"
      )
      rfx_component_means <- apply(rfx_beta_samples, 1, mean)
      rfx_component_sds <- apply(rfx_beta_samples, 1, sd)
      cat("Variance component means: ", rfx_component_means, "\n")
      cat("Variance component standard deviations: ", rfx_component_sds, "\n")
      quantile_summary <- t(apply(
        rfx_beta_samples,
        1,
        quantile,
        probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
      ))
      cat("Variance component quantiles:\n")
      print(quantile_summary)
    } else {
      cat(
        "Random effects summary of variance components across groups and posterior draws:\n"
      )
      rfx_component_means <- mean(rfx_beta_samples)
      rfx_component_sds <- sd(rfx_beta_samples)
      cat("Random effects overall mean: ", rfx_component_means, "\n")
      cat(
        "Random effects overall standard deviation: ",
        rfx_component_sds,
        "\n"
      )
      cat("Random effects overall quantiles:\n")
      quantile_summary <- quantile(
        rfx_beta_samples,
        probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
      )
      cat("Random effects overall quantiles:\n")
      print(quantile_summary)
    }
  }

  # Return bart_model invisibly
  invisible(object)
}

#' @title Plot BART Model Fit
#' @description Plot the BART model fit and any relevant sampled quantities. This will default to a traceplot of the global error scale and the in-sample mean forest predictions for the first train set observation. Since `nashtree::bart()` is flexible and it's possible to sample a model with a fixed global error scale and no mean forest, this procedure is adaptive and will attempt to plot a trace of whichever model terms are included if these two default terms are omitted.
#' @param x The BART model object
#' @param ... Additional arguments
#' @export
#' @return BART model object unchanged after summarizing
plot.bartmodel <- function(x, ...) {
  # Check if model has global error scale samples
  has_sigma2_samples <- x$model_params$sample_sigma2_global
  has_mean_forest_preds <- !is.null(x$y_hat_train)

  # Check if model is ordinal / binary
  is_probit <- (x$model_params$outcome_model$link == "probit" &&
    x$model_params$outcome_model$outcome == "binary")
  is_binary_cloglog <- (x$model_params$outcome_model$link == "cloglog" &&
    x$model_params$outcome_model$outcome == "binary")
  is_ordinal_cloglog <- (x$model_params$outcome_model$link == "cloglog" &&
    x$model_params$outcome_model$outcome == "ordinal")
  non_continuous_outcome <- (is_probit ||
    is_binary_cloglog ||
    is_ordinal_cloglog)

  # First try combinations of sigma2 and mean forest predictions
  if (has_sigma2_samples || has_mean_forest_preds) {
    if (has_sigma2_samples) {
      plot(
        x$sigma2_global_samples,
        type = "l",
        ylab = "Sigma^2",
        main = "Global error scale traceplot"
      )
    } else if (has_mean_forest_preds) {
      if (non_continuous_outcome) {
        plot_text <- "In-sample inverse-link-scale prediction trace for the first train set observation"
      } else {
        plot_text <- "In-sample mean function trace for the first train set observation"
      }
      plot(
        x$y_hat_train[1, ],
        type = "l",
        ylab = "Predictions",
        main = plot_text
      )
    }
  } else {
    stop(
      "This model does not have enough model terms / parameter traces to produce stochtree's default plots. See `predict.bartmodel()` for examples of how to further investigate your model."
    )
  }

  # Return x invisibly
  invisible(x)
}

#' @title Extract Random Effects Samples from BART Model
#' @description
#' Extract raw sample values for each of the random effect parameter terms.
#'
#' @param object Object of type `bartmodel` containing draws of a BART model and associated sampling outputs.
#' @param ... Other parameters to be used in random effects extraction
#' @return List of arrays. The alpha array has dimension (`num_components`, `num_samples`) and is simply a vector if `num_components = 1`.
#' The xi and beta arrays have dimension (`num_components`, `num_groups`, `num_samples`) and is simply a matrix if `num_components = 1`.
#' The sigma array has dimension (`num_components`, `num_samples`) and is simply a vector if `num_components = 1`.
#' @export
#'
#' @examples
#' n <- 100
#' p <- 5
#' X <- matrix(runif(n*p), ncol = p)
#' f_XW <- (
#'     ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) +
#'     ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) +
#'     ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) +
#'     ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
#' )
#' snr <- 3
#' group_ids <- rep(c(1,2), n %/% 2)
#' rfx_coefs <- matrix(c(-1, -1, 1, 1), nrow=2, byrow=TRUE)
#' rfx_basis <- cbind(1, runif(n, -1, 1))
#' rfx_term <- rowSums(rfx_coefs[group_ids,] * rfx_basis)
#' E_y <- f_XW + rfx_term
#' y <- E_y + rnorm(n, 0, 1)*(sd(E_y)/snr)
#' test_set_pct <- 0.2
#' n_test <- round(test_set_pct*n)
#' n_train <- n - n_test
#' test_inds <- sort(sample(1:n, n_test, replace = FALSE))
#' train_inds <- (1:n)[!((1:n) %in% test_inds)]
#' X_test <- X[test_inds,]
#' X_train <- X[train_inds,]
#' y_test <- y[test_inds]
#' y_train <- y[train_inds]
#' rfx_group_ids_test <- group_ids[test_inds]
#' rfx_group_ids_train <- group_ids[train_inds]
#' rfx_basis_test <- rfx_basis[test_inds,]
#' rfx_basis_train <- rfx_basis[train_inds,]
#' rfx_term_test <- rfx_term[test_inds]
#' rfx_term_train <- rfx_term[train_inds]
#' bart_model <- bart(X_train = X_train, y_train = y_train, X_test = X_test,
#'                    rfx_group_ids_train = rfx_group_ids_train,
#'                    rfx_group_ids_test = rfx_group_ids_test,
#'                    rfx_basis_train = rfx_basis_train,
#'                    rfx_basis_test = rfx_basis_test,
#'                    num_gfr = 10, num_burnin = 0, num_mcmc = 10)
#' rfx_samples <- getRandomEffectSamples(bart_model)
getRandomEffectSamples.bartmodel <- function(object, ...) {
  result = list()

  if (!object$model_params$has_rfx) {
    warning("This model has no RFX terms, returning an empty list")
    return(result)
  }

  # Extract the samples
  result <- object$rfx_samples$extract_parameter_samples()

  # Scale by sd(y_train)
  result$beta_samples <- result$beta_samples *
    object$model_params$outcome_scale
  result$xi_samples <- result$xi_samples * object$model_params$outcome_scale
  result$alpha_samples <- result$alpha_samples *
    object$model_params$outcome_scale
  result$sigma_samples <- result$sigma_samples *
    (object$model_params$outcome_scale^2)

  return(result)
}

#' @title Extract BART Parameter Samples
#' @description Extract a vector, matrix or array of parameter samples from a BART model by name.
#' Random effects are handled by a separate `getRandomEffectSamples` function due to the complexity of the random effects parameters.
#' If the requested model term is not found, an error is thrown.
#' The following conventions are used for parameter names:
#' - Global error variance: `"sigma2"`, `"global_error_scale"`, `"sigma2_global"`
#' - Leaf scale: `"sigma2_leaf"`, `"leaf_scale"`
#' - In-sample mean function predictions: `"y_hat_train"`
#' - Test set mean function predictions: `"y_hat_test"`
#' - In-sample variance forest predictions: `"sigma2_x_train"`, `"var_x_train"`
#' - Test set variance forest predictions: `"sigma2_x_test"`, `"var_x_test"`
#' - Ordinal model cutpoints (valid only for ordinal cloglog models): `"cloglog_cutpoints"`, `"cutpoints"`
#' - Negative-binomial size (dispersion) parameter `r` (valid only for count models): `"negbin_dispersion"`, `"dispersion"`
#' - Zero-inflated count models: probit index of the zero probability, in sample and on the test set: `"zero_hat_train"`, `"zero_hat_test"`; the zero probability `pi` itself: `"zero_probability_train"`, `"zero_probability_test"`; zero forest leaf scale draws (if sampled): `"sigma2_leaf_zero"`
#'
#' @param object Object of type `bartmodel` containing draws of a BART model and associated sampling outputs.
#' @param term Name of the parameter to extract (e.g., `"sigma2"`, `"y_hat_train"`, etc.)
#' @return Array of parameter samples. If the underlying parameter is a scalar, this will be a vector of length `num_samples`.
#' If the underlying parameter is vector-valued, this will be (`parameter_dimension` x `num_samples`) matrix, and if the underlying
#' parameter is multidimensional, this will be an array of dimension (`parameter_dimension_1` x `parameter_dimension_2` x ... x `num_samples`).
#' @export
#'
#' @examples
#' n <- 100
#' p <- 5
#' X <- matrix(runif(n*p), ncol = p)
#' f_XW <- (
#'     ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) +
#'     ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) +
#'     ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) +
#'     ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
#' )
#' snr <- 3
#' group_ids <- rep(c(1,2), n %/% 2)
#' rfx_coefs <- matrix(c(-1, -1, 1, 1), nrow=2, byrow=TRUE)
#' rfx_basis <- cbind(1, runif(n, -1, 1))
#' rfx_term <- rowSums(rfx_coefs[group_ids,] * rfx_basis)
#' E_y <- f_XW + rfx_term
#' y <- E_y + rnorm(n, 0, 1)*(sd(E_y)/snr)
#' test_set_pct <- 0.2
#' n_test <- round(test_set_pct*n)
#' n_train <- n - n_test
#' test_inds <- sort(sample(1:n, n_test, replace = FALSE))
#' train_inds <- (1:n)[!((1:n) %in% test_inds)]
#' X_test <- X[test_inds,]
#' X_train <- X[train_inds,]
#' y_test <- y[test_inds]
#' y_train <- y[train_inds]
#' rfx_group_ids_test <- group_ids[test_inds]
#' rfx_group_ids_train <- group_ids[train_inds]
#' rfx_basis_test <- rfx_basis[test_inds,]
#' rfx_basis_train <- rfx_basis[train_inds,]
#' rfx_term_test <- rfx_term[test_inds]
#' rfx_term_train <- rfx_term[train_inds]
#' bart_model <- bart(X_train = X_train, y_train = y_train, X_test = X_test,
#'                    rfx_group_ids_train = rfx_group_ids_train,
#'                    rfx_group_ids_test = rfx_group_ids_test,
#'                    rfx_basis_train = rfx_basis_train,
#'                    rfx_basis_test = rfx_basis_test,
#'                    num_gfr = 10, num_burnin = 0, num_mcmc = 10)
#' sigma2_samples <- extractParameter(bart_model, "sigma2")
extractParameter.bartmodel <- function(object, term) {
  if (term %in% c("sigma2", "global_error_scale", "sigma2_global")) {
    if (!is.null(object$sigma2_global_samples)) {
      return(object$sigma2_global_samples)
    } else {
      stop("This model does not have global variance parameter samples")
    }
  }

  if (term %in% c("sigma2_leaf", "leaf_scale")) {
    if (!is.null(object$sigma2_leaf_samples)) {
      return(object$sigma2_leaf_samples)
    } else {
      stop("This model does not have leaf variance parameter samples")
    }
  }

  if (term %in% c("y_hat_train")) {
    if (!is.null(object$y_hat_train)) {
      return(object$y_hat_train)
    } else {
      stop(
        "This model does not have in-sample mean function prediction samples"
      )
    }
  }

  if (term %in% c("y_hat_test")) {
    if (!is.null(object$y_hat_test)) {
      return(object$y_hat_test)
    } else {
      stop("This model does not have test set mean function prediction samples")
    }
  }

  if (term %in% c("sigma2_x_train", "var_x_train")) {
    if (!is.null(object$sigma2_x_hat_train)) {
      return(object$sigma2_x_hat_train)
    } else {
      stop("This model does not have in-sample variance forest predictions")
    }
  }

  if (term %in% c("sigma2_x_test", "var_x_test")) {
    if (!is.null(object$sigma2_x_hat_test)) {
      return(object$sigma2_x_hat_test)
    } else {
      stop("This model does not have test set variance forest predictions")
    }
  }

  if (term %in% c("cloglog_cutpoints", "cutpoints")) {
    if (!is.null(object$cloglog_cutpoint_samples)) {
      return(object$cloglog_cutpoint_samples)
    } else {
      stop("This model does not have ordinal cutpoint samples")
    }
  }

  if (term %in% c("zero_hat_train", "zero_probability_train")) {
    if (!is.null(object$zero_hat_train)) {
      zero_hat <- object$zero_hat_train
      return(if (term == "zero_hat_train") zero_hat else pnorm(zero_hat))
    } else {
      stop("This model does not have in-sample zero forest predictions")
    }
  }

  if (term %in% c("zero_hat_test", "zero_probability_test")) {
    if (!is.null(object$zero_hat_test)) {
      zero_hat <- object$zero_hat_test
      return(if (term == "zero_hat_test") zero_hat else pnorm(zero_hat))
    } else {
      stop("This model does not have test set zero forest predictions")
    }
  }

  if (term == "sigma2_leaf_zero") {
    if (!is.null(object$sigma2_leaf_zero_samples)) {
      return(object$sigma2_leaf_zero_samples)
    } else {
      stop("This model does not have zero forest leaf scale samples")
    }
  }

  if (term %in% c("negbin_dispersion", "dispersion")) {
    if (!is.null(object$negbin_dispersion_samples)) {
      return(object$negbin_dispersion_samples)
    } else {
      stop("This model does not have negative-binomial dispersion samples")
    }
  }

  stop(paste0("term ", term, " is not a valid BART model term"))
}

#' @title Convert BART Model to JSON
#' @rdname BARTSerialization
#' @param object Object of type `bartmodel` containing draws of a BART model and associated sampling outputs.
#' @export
saveBARTModelToJson <- function(object) {
  jsonobj <- createCppJson()

  if (!inherits(object, "bartmodel")) {
    stop("`object` must be a BART model")
  }

  if (is.null(object$model_params)) {
    stop("This BCF model has not yet been sampled")
  }

  # Add the forests
  if (object$model_params$include_mean_forest) {
    jsonobj$add_forest(object$mean_forests)
  }
  if (object$model_params$include_variance_forest) {
    jsonobj$add_forest(object$variance_forests)
  }
  is_zinb_model <- identical(object$model_params$outcome_model$family, "zinb")
  if (is_zinb_model) {
    jsonobj$add_forest(object$zero_forests)
  }

  # Add metadata
  jsonobj$add_scalar(
    "num_numeric_vars",
    object$train_set_metadata$num_numeric_vars
  )
  jsonobj$add_scalar(
    "num_ordered_cat_vars",
    object$train_set_metadata$num_ordered_cat_vars
  )
  jsonobj$add_scalar(
    "num_unordered_cat_vars",
    object$train_set_metadata$num_unordered_cat_vars
  )
  if (object$train_set_metadata$num_numeric_vars > 0) {
    jsonobj$add_string_vector(
      "numeric_vars",
      object$train_set_metadata$numeric_vars
    )
  }
  if (object$train_set_metadata$num_ordered_cat_vars > 0) {
    jsonobj$add_string_vector(
      "ordered_cat_vars",
      object$train_set_metadata$ordered_cat_vars
    )
    jsonobj$add_string_list(
      "ordered_unique_levels",
      object$train_set_metadata$ordered_unique_levels
    )
  }
  if (object$train_set_metadata$num_unordered_cat_vars > 0) {
    jsonobj$add_string_vector(
      "unordered_cat_vars",
      object$train_set_metadata$unordered_cat_vars
    )
    jsonobj$add_string_list(
      "unordered_unique_levels",
      object$train_set_metadata$unordered_unique_levels
    )
  }

  # Add version stamp and global parameters
  jsonobj$add_string("stochtree_version", getStochtreeVersion())
  jsonobj$add_scalar("outcome_scale", object$model_params$outcome_scale)
  jsonobj$add_scalar("outcome_mean", object$model_params$outcome_mean)
  jsonobj$add_boolean("standardize", object$model_params$standardize)
  jsonobj$add_scalar("sigma2_init", object$model_params$sigma2_init)
  jsonobj$add_boolean(
    "sample_sigma2_global",
    object$model_params$sample_sigma2_global
  )
  jsonobj$add_boolean(
    "sample_sigma2_leaf",
    object$model_params$sample_sigma2_leaf
  )
  jsonobj$add_boolean(
    "include_mean_forest",
    object$model_params$include_mean_forest
  )
  jsonobj$add_boolean(
    "include_variance_forest",
    object$model_params$include_variance_forest
  )
  jsonobj$add_boolean("has_rfx", object$model_params$has_rfx)
  jsonobj$add_boolean("has_rfx_basis", object$model_params$has_rfx_basis)
  jsonobj$add_scalar("num_rfx_basis", object$model_params$num_rfx_basis)
  jsonobj$add_scalar("num_gfr", object$model_params$num_gfr)
  jsonobj$add_scalar("num_burnin", object$model_params$num_burnin)
  jsonobj$add_scalar("num_mcmc", object$model_params$num_mcmc)
  jsonobj$add_scalar("num_samples", object$model_params$num_samples)
  jsonobj$add_scalar("num_covariates", object$model_params$num_covariates)
  jsonobj$add_scalar("num_basis", object$model_params$num_basis)
  jsonobj$add_scalar("num_chains", object$model_params$num_chains)
  jsonobj$add_scalar("keep_every", object$model_params$keep_every)
  jsonobj$add_boolean("requires_basis", object$model_params$requires_basis)
  jsonobj$add_string(
    "outcome",
    object$model_params$outcome_model$outcome,
    "outcome_model"
  )
  jsonobj$add_string(
    "link",
    object$model_params$outcome_model$link,
    "outcome_model"
  )
  jsonobj$add_boolean(
    "probit_outcome_model",
    object$model_params$probit_outcome_model
  )
  jsonobj$add_string(
    "rfx_model_spec",
    object$model_params$rfx_model_spec
  )
  jsonobj$add_string(
    "rfx_unseen_groups",
    validateRfxUnseenGroups(object$model_params[["rfx_unseen_groups"]])
  )
  if (object$model_params$outcome_model$link == "cloglog") {
    jsonobj$add_scalar(
      "cloglog_num_categories",
      object$model_params$cloglog_num_categories
    )
    if (object$model_params$outcome_model$outcome == "ordinal") {
      for (i in 1:(object$model_params$cloglog_num_categories - 1)) {
        jsonobj$add_vector(
          paste0("cloglog_cutpoint_samples_", i),
          object$cloglog_cutpoint_samples[i, ],
          "parameters"
        )
      }
    }
  }
  if (identical(object$model_params$outcome_model$outcome, "count")) {
    jsonobj$add_string(
      "family",
      object$model_params$outcome_model$family,
      "outcome_model"
    )
    if (is_zinb_model) {
      jsonobj$add_scalar(
        "zero_outcome_mean",
        object$model_params$zero_outcome_mean
      )
      jsonobj$add_scalar(
        "zinb_zero_probability_init",
        object$model_params$zinb_zero_probability_init
      )
      jsonobj$add_scalar("num_trees_zero", object$model_params$num_trees_zero)
      jsonobj$add_boolean(
        "sample_sigma2_leaf_zero",
        isTRUE(object$model_params$sample_sigma2_leaf_zero)
      )
      jsonobj$add_scalar(
        "sigma2_leaf_zero_init",
        object$model_params$sigma2_leaf_zero_init
      )
      if (isTRUE(object$model_params$sample_sigma2_leaf_zero)) {
        jsonobj$add_vector(
          "sigma2_leaf_zero_samples",
          object$sigma2_leaf_zero_samples,
          "parameters"
        )
      }
      jsonobj$add_boolean(
        "has_zero_rfx",
        isTRUE(object$model_params$has_zero_rfx)
      )
    }
    jsonobj$add_boolean("has_offset", isTRUE(object$model_params$has_offset))
    jsonobj$add_scalar(
      "negbin_dispersion_prior_shape",
      object$model_params$negbin_dispersion_prior_shape
    )
    jsonobj$add_scalar(
      "negbin_dispersion_prior_rate",
      object$model_params$negbin_dispersion_prior_rate
    )
    jsonobj$add_scalar(
      "negbin_dispersion_init",
      object$model_params$negbin_dispersion_init
    )
    jsonobj$add_boolean(
      "sample_negbin_dispersion",
      isTRUE(object$model_params$sample_negbin_dispersion)
    )
    jsonobj$add_vector(
      "negbin_dispersion_samples",
      object$negbin_dispersion_samples,
      "parameters"
    )
  }
  if (object$model_params$sample_sigma2_global) {
    jsonobj$add_vector(
      "sigma2_global_samples",
      object$sigma2_global_samples,
      "parameters"
    )
  }
  if (object$model_params$sample_sigma2_leaf) {
    jsonobj$add_vector(
      "sigma2_leaf_samples",
      object$sigma2_leaf_samples,
      "parameters"
    )
  }

  # Add random effects (if present)
  if (object$model_params$has_rfx) {
    jsonobj$add_random_effects(object$rfx_samples)
    jsonobj$add_string_vector(
      "rfx_unique_group_ids",
      object$rfx_unique_group_ids
    )
  }
  # Zero-part random effects (zero-inflated family): the second random
  # effects container, random_effect_container_1
  if (isTRUE(object$model_params$has_zero_rfx)) {
    jsonobj$add_random_effects(object$zero_rfx_samples)
  }

  # Add covariate preprocessor metadata
  preprocessor_metadata_string <- savePreprocessorToJsonString(
    object$train_set_metadata
  )
  jsonobj$add_string("preprocessor_metadata", preprocessor_metadata_string)

  return(jsonobj)
}

#' @title Save BART Model to JSON File
#' @rdname BARTSerialization
#' @param object Object of type `bartmodel` containing draws of a BART model and associated sampling outputs.
#' @param filename String of filepath, must end in ".json"
#'
#' @export
saveBARTModelToJsonFile <- function(object, filename) {
  # Convert to Json
  jsonobj <- saveBARTModelToJson(object)

  # Save to file
  jsonobj$save_file(filename)
}

#' @title Convert BART Model to JSON String
#' @rdname BARTSerialization
#' @param object Object of type `bartmodel` containing draws of a BART model and associated sampling outputs.
#' @export
saveBARTModelToJsonString <- function(object) {
  # Convert to Json
  jsonobj <- saveBARTModelToJson(object)

  # Dump to string
  return(jsonobj$return_json_string())
}

# Reconstruct mean-forest leaf-model metadata that is not serialized directly.
# These fields are fully determined by num_basis (a basis/regression leaf model
# has num_basis > 0), so we derive them on load rather than store redundant
# copies. Restoring them keeps print()/summary()/predict() working on a
# reloaded model (they were previously NULL and errored).
.reconstructBartLeafModelFields <- function(model_params) {
  num_basis <- model_params[["num_basis"]]
  has_basis <- isTRUE(num_basis > 0)
  model_params[["has_basis"]] <- has_basis
  model_params[["leaf_regression"]] <- has_basis
  model_params[["is_leaf_constant"]] <- !has_basis
  model_params[["leaf_dimension"]] <- if (has_basis) num_basis else 1
  model_params
}

#' @title Convert JSON to BART Model
#' @rdname BARTSerialization
#' @param json_object Object of type `CppJson` containing Json representation of a BART model
#' @export
createBARTModelFromJson <- function(json_object) {
  # Initialize the BCF model
  output <- list()

  # Helpers for optional-field presence checks
  .ver <- inferStochtreeJsonVersion(json_object)
  has_field <- function(name) {
    json_contains_field_cpp(json_object$json_ptr, name)
  }
  has_subfolder_field <- function(subfolder, name) {
    json_contains_field_subfolder_cpp(json_object$json_ptr, subfolder, name)
  }

  # Unpack the forests
  include_mean_forest <- json_object$get_boolean("include_mean_forest")
  include_variance_forest <- json_object$get_boolean(
    "include_variance_forest"
  )
  if (include_mean_forest) {
    output[["mean_forests"]] <- loadForestContainerJson(
      json_object,
      "forest_0"
    )
    if (include_variance_forest) {
      output[["variance_forests"]] <- loadForestContainerJson(
        json_object,
        "forest_1"
      )
    }
  } else {
    output[["variance_forests"]] <- loadForestContainerJson(
      json_object,
      "forest_0"
    )
  }

  # Unpack metadata
  train_set_metadata = list()
  train_set_metadata[["num_numeric_vars"]] <- json_object$get_scalar(
    "num_numeric_vars"
  )
  train_set_metadata[["num_ordered_cat_vars"]] <- json_object$get_scalar(
    "num_ordered_cat_vars"
  )
  train_set_metadata[["num_unordered_cat_vars"]] <- json_object$get_scalar(
    "num_unordered_cat_vars"
  )
  if (train_set_metadata[["num_numeric_vars"]] > 0) {
    train_set_metadata[["numeric_vars"]] <- json_object$get_string_vector(
      "numeric_vars"
    )
  }
  if (train_set_metadata[["num_ordered_cat_vars"]] > 0) {
    train_set_metadata[[
      "ordered_cat_vars"
    ]] <- json_object$get_string_vector("ordered_cat_vars")
    train_set_metadata[[
      "ordered_unique_levels"
    ]] <- json_object$get_string_list(
      "ordered_unique_levels",
      train_set_metadata[["ordered_cat_vars"]]
    )
  }
  if (train_set_metadata[["num_unordered_cat_vars"]] > 0) {
    train_set_metadata[[
      "unordered_cat_vars"
    ]] <- json_object$get_string_vector("unordered_cat_vars")
    train_set_metadata[[
      "unordered_unique_levels"
    ]] <- json_object$get_string_list(
      "unordered_unique_levels",
      train_set_metadata[["unordered_cat_vars"]]
    )
  }
  output[["train_set_metadata"]] <- train_set_metadata

  # Unpack model params
  model_params = list()
  model_params[["outcome_scale"]] <- json_object$get_scalar("outcome_scale")
  model_params[["outcome_mean"]] <- json_object$get_scalar("outcome_mean")
  model_params[["standardize"]] <- json_object$get_boolean("standardize")
  model_params[["sigma2_init"]] <- json_object$get_scalar("sigma2_init")
  model_params[["sample_sigma2_global"]] <- json_object$get_boolean(
    "sample_sigma2_global"
  )
  model_params[["sample_sigma2_leaf"]] <- json_object$get_boolean(
    "sample_sigma2_leaf"
  )
  model_params[["include_mean_forest"]] <- include_mean_forest
  model_params[["include_variance_forest"]] <- include_variance_forest
  model_params[["has_rfx"]] <- json_object$get_boolean("has_rfx")

  if (has_field("has_rfx_basis")) {
    model_params[["has_rfx_basis"]] <- json_object$get_boolean("has_rfx_basis")
    model_params[["num_rfx_basis"]] <- json_object$get_scalar("num_rfx_basis")
  } else {
    model_params[["has_rfx_basis"]] <- FALSE
    model_params[["num_rfx_basis"]] <- 1
    warning(paste0(
      "Fields 'has_rfx_basis' and 'num_rfx_basis' not found in JSON (model appears to have been ",
      "serialized under stochtree ",
      .ver,
      "). Defaulting to FALSE / 1. ",
      "Re-save your model to suppress this warning."
    ))
  }

  model_params[["num_gfr"]] <- json_object$get_scalar("num_gfr")
  model_params[["num_burnin"]] <- json_object$get_scalar("num_burnin")
  model_params[["num_mcmc"]] <- json_object$get_scalar("num_mcmc")
  model_params[["num_samples"]] <- json_object$get_scalar("num_samples")
  model_params[["num_covariates"]] <- if (has_field("num_covariates")) {
    json_object$get_scalar("num_covariates")
  } else {
    NA_real_
  }
  model_params[["num_basis"]] <- json_object$get_scalar("num_basis")
  model_params <- .reconstructBartLeafModelFields(model_params)
  model_params[["requires_basis"]] <- json_object$get_boolean("requires_basis")

  if (has_field("num_chains")) {
    model_params[["num_chains"]] <- json_object$get_scalar("num_chains")
  } else {
    model_params[["num_chains"]] <- 1
    warning(paste0(
      "Field 'num_chains' not found in JSON (model appears to have been serialized under stochtree ",
      .ver,
      "). Defaulting to 1. Re-save your model to suppress this warning."
    ))
  }

  if (has_field("keep_every")) {
    model_params[["keep_every"]] <- json_object$get_scalar("keep_every")
  } else {
    model_params[["keep_every"]] <- 1
    warning(paste0(
      "Field 'keep_every' not found in JSON (model appears to have been serialized under stochtree ",
      .ver,
      "). Defaulting to 1. Re-save your model to suppress this warning."
    ))
  }

  model_params[["probit_outcome_model"]] <- if (
    has_field("probit_outcome_model")
  ) {
    json_object$get_boolean("probit_outcome_model")
  } else {
    FALSE
  }

  if (
    has_subfolder_field("outcome_model", "outcome") &&
      has_subfolder_field("outcome_model", "link")
  ) {
    outcome_model_outcome <- json_object$get_string("outcome", "outcome_model")
    outcome_model_link <- json_object$get_string("link", "outcome_model")
  } else {
    outcome_model_outcome <- "continuous"
    outcome_model_link <- "identity"
    warning(paste0(
      "Fields 'outcome' and 'link' not found under 'outcome_model' in JSON (model appears to have ",
      "been serialized under stochtree ",
      .ver,
      "). Defaulting to outcome='continuous', ",
      "link='identity'. Re-save your model to suppress this warning."
    ))
  }
  model_params[["outcome_model"]] <- OutcomeModel(
    outcome = outcome_model_outcome,
    link = outcome_model_link,
    family = .readOutcomeFamilyFromJson(json_object, outcome_model_outcome)
  )

  if (has_field("rfx_model_spec")) {
    model_params[["rfx_model_spec"]] <- json_object$get_string("rfx_model_spec")
  } else {
    model_params[["rfx_model_spec"]] <- ""
    if (model_params[["has_rfx"]]) {
      warning(paste0(
        "Field 'rfx_model_spec' not found in JSON (model appears to have been serialized under ",
        "stochtree ",
        .ver,
        "). Defaulting to ''. Re-save your model to suppress this warning."
      ))
    }
  }
  if (has_field("rfx_unseen_groups")) {
    model_params[["rfx_unseen_groups"]] <- json_object$get_string(
      "rfx_unseen_groups"
    )
  } else {
    # Models serialized before this option always refused unseen group labels.
    model_params[["rfx_unseen_groups"]] <- "error"
  }
  if (model_params[["outcome_model"]]$link == "cloglog") {
    cloglog_num_categories <- json_object$get_scalar("cloglog_num_categories")
    model_params[["cloglog_num_categories"]] <- cloglog_num_categories
  } else {
    model_params[["cloglog_num_categories"]] <- 0
  }

  if (identical(model_params[["outcome_model"]]$outcome, "count")) {
    negbin_fields <- .loadBartNegbinFieldsFromJson(
      list(json_object),
      model_params,
      output
    )
    model_params <- negbin_fields$model_params
    output <- negbin_fields$output
    if (identical(model_params[["outcome_model"]]$family, "zinb")) {
      zinb_fields <- .loadBartZinbFieldsFromJson(
        list(json_object),
        model_params,
        output
      )
      model_params <- zinb_fields$model_params
      output <- zinb_fields$output
      output[["zero_forests"]] <- loadForestContainerJson(
        json_object,
        .zeroForestJsonLabel(include_mean_forest, include_variance_forest)
      )
    }
  }
  output[["model_params"]] <- model_params

  # Unpack sampled parameters
  if (model_params[["sample_sigma2_global"]]) {
    output[["sigma2_global_samples"]] <- json_object$get_vector(
      "sigma2_global_samples",
      "parameters"
    )
  }
  if (model_params[["sample_sigma2_leaf"]]) {
    output[["sigma2_leaf_samples"]] <- json_object$get_vector(
      "sigma2_leaf_samples",
      "parameters"
    )
  }
  if (
    model_params[["outcome_model"]]$link == "cloglog" &&
      model_params[["outcome_model"]]$outcome == "ordinal"
  ) {
    cloglog_cutpoint_samples <- matrix(
      NA_real_,
      model_params[["cloglog_num_categories"]] - 1,
      model_params[["num_samples"]]
    )
    for (i in 1:(model_params[["cloglog_num_categories"]] - 1)) {
      cloglog_cutpoint_samples[i, ] <- json_object$get_vector(
        paste0("cloglog_cutpoint_samples_", i),
        "parameters"
      )
    }
    output[["cloglog_cutpoint_samples"]] <- cloglog_cutpoint_samples
  }

  # Unpack random effects
  if (model_params[["has_rfx"]]) {
    output[["rfx_unique_group_ids"]] <- json_object$get_string_vector(
      "rfx_unique_group_ids"
    )
    output[["rfx_samples"]] <- loadRandomEffectSamplesJson(json_object, 0)
  }

  # Unpack covariate preprocessor
  if (has_field("preprocessor_metadata")) {
    preprocessor_metadata_string <- json_object$get_string(
      "preprocessor_metadata"
    )
    output[["train_set_metadata"]] <- createPreprocessorFromJsonString(
      preprocessor_metadata_string
    )
  } else {
    output[["train_set_metadata"]] <- NULL
    warning(paste0(
      "Field 'preprocessor_metadata' not found in JSON (model appears to have been serialized ",
      "under stochtree ",
      .ver,
      "). DataFrame covariates will not be supported for prediction. ",
      "Re-save your model to suppress this warning."
    ))
  }

  class(output) <- "bartmodel"
  return(output)
}

#' @title Convert JSON File to BART Model
#' @rdname BARTSerialization
#' @param json_filename String of filepath, must end in ".json"
#' @export
createBARTModelFromJsonFile <- function(json_filename) {
  # Load a `CppJson` object from file
  bart_json <- createCppJsonFile(json_filename)

  # Create and return the BART object
  bart_object <- createBARTModelFromJson(bart_json)

  return(bart_object)
}

#' @title Convert JSON String to BART Model
#' @rdname BARTSerialization
#' @param json_string JSON string dump
#' @export
createBARTModelFromJsonString <- function(json_string) {
  # Load a `CppJson` object from string
  bart_json <- createCppJsonString(json_string)

  # Create and return the BART object
  bart_object <- createBARTModelFromJson(bart_json)

  return(bart_object)
}

#' @title Convert JSON List to Single BART Model
#' @rdname BARTSerialization
#' @param json_object_list List of objects of type `CppJson` containing Json representation of a BART model
#' @export
createBARTModelFromCombinedJson <- function(json_object_list) {
  # Initialize the BCF model
  output <- list()

  # For scalar / preprocessing details which aren't sample-dependent,
  # defer to the first json
  json_object_default <- json_object_list[[1]]

  # Helpers for optional-field presence checks
  .ver <- inferStochtreeJsonVersion(json_object_default)
  has_field <- function(name) {
    json_contains_field_cpp(json_object_default$json_ptr, name)
  }
  has_subfolder_field <- function(subfolder, name) {
    json_contains_field_subfolder_cpp(
      json_object_default$json_ptr,
      subfolder,
      name
    )
  }

  # Unpack the forests
  include_mean_forest <- json_object_default$get_boolean(
    "include_mean_forest"
  )
  include_variance_forest <- json_object_default$get_boolean(
    "include_variance_forest"
  )
  if (include_mean_forest) {
    output[["mean_forests"]] <- loadForestContainerCombinedJson(
      json_object_list,
      "forest_0"
    )
    if (include_variance_forest) {
      output[["variance_forests"]] <- loadForestContainerCombinedJson(
        json_object_list,
        "forest_1"
      )
    }
  } else {
    output[["variance_forests"]] <- loadForestContainerCombinedJson(
      json_object_list,
      "forest_0"
    )
  }

  # Unpack metadata
  train_set_metadata = list()
  train_set_metadata[["num_numeric_vars"]] <- json_object_default$get_scalar(
    "num_numeric_vars"
  )
  train_set_metadata[[
    "num_ordered_cat_vars"
  ]] <- json_object_default$get_scalar("num_ordered_cat_vars")
  train_set_metadata[[
    "num_unordered_cat_vars"
  ]] <- json_object_default$get_scalar("num_unordered_cat_vars")
  if (train_set_metadata[["num_numeric_vars"]] > 0) {
    train_set_metadata[[
      "numeric_vars"
    ]] <- json_object_default$get_string_vector("numeric_vars")
  }
  if (train_set_metadata[["num_ordered_cat_vars"]] > 0) {
    train_set_metadata[[
      "ordered_cat_vars"
    ]] <- json_object_default$get_string_vector("ordered_cat_vars")
    train_set_metadata[[
      "ordered_unique_levels"
    ]] <- json_object_default$get_string_list(
      "ordered_unique_levels",
      train_set_metadata[["ordered_cat_vars"]]
    )
  }
  if (train_set_metadata[["num_unordered_cat_vars"]] > 0) {
    train_set_metadata[[
      "unordered_cat_vars"
    ]] <- json_object_default$get_string_vector("unordered_cat_vars")
    train_set_metadata[[
      "unordered_unique_levels"
    ]] <- json_object_default$get_string_list(
      "unordered_unique_levels",
      train_set_metadata[["unordered_cat_vars"]]
    )
  }
  output[["train_set_metadata"]] <- train_set_metadata

  # Unpack model params
  model_params = list()
  model_params[["outcome_scale"]] <- json_object_default$get_scalar(
    "outcome_scale"
  )
  model_params[["outcome_mean"]] <- json_object_default$get_scalar(
    "outcome_mean"
  )
  model_params[["standardize"]] <- json_object_default$get_boolean(
    "standardize"
  )
  model_params[["sigma2_init"]] <- json_object_default$get_scalar(
    "sigma2_init"
  )
  model_params[["sample_sigma2_global"]] <- json_object_default$get_boolean(
    "sample_sigma2_global"
  )
  model_params[["sample_sigma2_leaf"]] <- json_object_default$get_boolean(
    "sample_sigma2_leaf"
  )
  model_params[["include_mean_forest"]] <- include_mean_forest
  model_params[["include_variance_forest"]] <- include_variance_forest
  model_params[["has_rfx"]] <- json_object_default$get_boolean("has_rfx")

  if (has_field("has_rfx_basis")) {
    model_params[["has_rfx_basis"]] <- json_object_default$get_boolean(
      "has_rfx_basis"
    )
    model_params[["num_rfx_basis"]] <- json_object_default$get_scalar(
      "num_rfx_basis"
    )
  } else {
    model_params[["has_rfx_basis"]] <- FALSE
    model_params[["num_rfx_basis"]] <- 1
    warning(paste0(
      "Fields 'has_rfx_basis' and 'num_rfx_basis' not found in JSON (model appears to have been ",
      "serialized under stochtree ",
      .ver,
      "). Defaulting to FALSE / 1. ",
      "Re-save your model to suppress this warning."
    ))
  }

  model_params[["num_covariates"]] <- if (has_field("num_covariates")) {
    json_object_default$get_scalar("num_covariates")
  } else {
    NA_real_
  }
  model_params[["num_basis"]] <- json_object_default$get_scalar("num_basis")
  model_params <- .reconstructBartLeafModelFields(model_params)
  model_params[["requires_basis"]] <- json_object_default$get_boolean(
    "requires_basis"
  )

  model_params[["probit_outcome_model"]] <- if (
    has_field("probit_outcome_model")
  ) {
    json_object_default$get_boolean("probit_outcome_model")
  } else {
    FALSE
  }

  if (
    has_subfolder_field("outcome_model", "outcome") &&
      has_subfolder_field("outcome_model", "link")
  ) {
    outcome_model_outcome <- json_object_default$get_string(
      "outcome",
      "outcome_model"
    )
    outcome_model_link <- json_object_default$get_string(
      "link",
      "outcome_model"
    )
  } else {
    outcome_model_outcome <- "continuous"
    outcome_model_link <- "identity"
    warning(paste0(
      "Fields 'outcome' and 'link' not found under 'outcome_model' in JSON (model appears to have ",
      "been serialized under stochtree ",
      .ver,
      "). Defaulting to outcome='continuous', ",
      "link='identity'. Re-save your model to suppress this warning."
    ))
  }
  model_params[["outcome_model"]] <- OutcomeModel(
    outcome = outcome_model_outcome,
    link = outcome_model_link,
    family = .readOutcomeFamilyFromJson(json_object_default, outcome_model_outcome)
  )

  if (has_field("rfx_model_spec")) {
    model_params[["rfx_model_spec"]] <- json_object_default$get_string(
      "rfx_model_spec"
    )
  } else {
    model_params[["rfx_model_spec"]] <- ""
    if (model_params[["has_rfx"]]) {
      warning(paste0(
        "Field 'rfx_model_spec' not found in JSON (model appears to have been serialized under ",
        "stochtree ",
        .ver,
        "). Defaulting to ''. Re-save your model to suppress this warning."
      ))
    }
  }
  if (has_field("rfx_unseen_groups")) {
    model_params[["rfx_unseen_groups"]] <- json_object_default$get_string(
      "rfx_unseen_groups"
    )
  } else {
    # Models serialized before this option always refused unseen group labels.
    model_params[["rfx_unseen_groups"]] <- "error"
  }

  if (has_field("num_chains")) {
    model_params[["num_chains"]] <- json_object_default$get_scalar("num_chains")
  } else {
    model_params[["num_chains"]] <- 1
    warning(paste0(
      "Field 'num_chains' not found in JSON (model appears to have been serialized under stochtree ",
      .ver,
      "). Defaulting to 1. Re-save your model to suppress this warning."
    ))
  }

  if (has_field("keep_every")) {
    model_params[["keep_every"]] <- json_object_default$get_scalar("keep_every")
  } else {
    model_params[["keep_every"]] <- 1
    warning(paste0(
      "Field 'keep_every' not found in JSON (model appears to have been serialized under stochtree ",
      .ver,
      "). Defaulting to 1. Re-save your model to suppress this warning."
    ))
  }
  if (model_params[["outcome_model"]]$link == "cloglog") {
    cloglog_num_categories <- json_object_default$get_scalar(
      "cloglog_num_categories"
    )
    model_params[["cloglog_num_categories"]] <- cloglog_num_categories
  } else {
    model_params[["cloglog_num_categories"]] <- 0
  }

  # Combine values that are sample-specific
  for (i in 1:length(json_object_list)) {
    json_object <- json_object_list[[i]]
    if (i == 1) {
      model_params[["num_gfr"]] <- json_object$get_scalar("num_gfr")
      model_params[["num_burnin"]] <- json_object$get_scalar("num_burnin")
      model_params[["num_mcmc"]] <- json_object$get_scalar("num_mcmc")
      model_params[["num_samples"]] <- json_object$get_scalar(
        "num_samples"
      )
    } else {
      prev_json <- json_object_list[[i - 1]]
      model_params[["num_gfr"]] <- model_params[["num_gfr"]] +
        json_object$get_scalar("num_gfr")
      model_params[["num_burnin"]] <- model_params[["num_burnin"]] +
        json_object$get_scalar("num_burnin")
      model_params[["num_mcmc"]] <- model_params[["num_mcmc"]] +
        json_object$get_scalar("num_mcmc")
      model_params[["num_samples"]] <- model_params[["num_samples"]] +
        json_object$get_scalar("num_samples")
    }
  }
  if (identical(model_params[["outcome_model"]]$outcome, "count")) {
    negbin_fields <- .loadBartNegbinFieldsFromJson(
      json_object_list,
      model_params,
      output
    )
    model_params <- negbin_fields$model_params
    output <- negbin_fields$output
    if (identical(model_params[["outcome_model"]]$family, "zinb")) {
      zinb_fields <- .loadBartZinbFieldsFromJson(
        json_object_list,
        model_params,
        output
      )
      model_params <- zinb_fields$model_params
      output <- zinb_fields$output
      output[["zero_forests"]] <- loadForestContainerCombinedJson(
        json_object_list,
        .zeroForestJsonLabel(include_mean_forest, include_variance_forest)
      )
    }
  }
  output[["model_params"]] <- model_params

  # Unpack sampled parameters
  if (model_params[["sample_sigma2_global"]]) {
    for (i in 1:length(json_object_list)) {
      json_object <- json_object_list[[i]]
      if (i == 1) {
        output[["sigma2_global_samples"]] <- json_object$get_vector(
          "sigma2_global_samples",
          "parameters"
        )
      } else {
        output[["sigma2_global_samples"]] <- c(
          output[["sigma2_global_samples"]],
          json_object$get_vector(
            "sigma2_global_samples",
            "parameters"
          )
        )
      }
    }
  }
  if (model_params[["sample_sigma2_leaf"]]) {
    for (i in 1:length(json_object_list)) {
      json_object <- json_object_list[[i]]
      if (i == 1) {
        output[["sigma2_leaf_samples"]] <- json_object$get_vector(
          "sigma2_leaf_samples",
          "parameters"
        )
      } else {
        output[["sigma2_leaf_samples"]] <- c(
          output[["sigma2_leaf_samples"]],
          json_object$get_vector("sigma2_leaf_samples", "parameters")
        )
      }
    }
  }
  if (
    model_params[["outcome_model"]]$link == "cloglog" &&
      model_params[["outcome_model"]]$outcome == "ordinal"
  ) {
    cloglog_cutpoint_samples <- matrix(
      NA_real_,
      model_params[["cloglog_num_categories"]] - 1,
      model_params[["num_samples"]]
    )
    index_start <- 1
    for (i in 1:length(json_object_list)) {
      json_object <- json_object_list[[i]]
      num_samples <- json_object$get_scalar("num_samples")
      subset_inds <- index_start:(index_start + num_samples - 1)
      for (j in 1:(model_params[["cloglog_num_categories"]] - 1)) {
        cloglog_cutpoint_samples[j, subset_inds] <- json_object$get_vector(
          paste0("cloglog_cutpoint_samples_", j),
          "parameters"
        )
      }
    }
    output[["cloglog_cutpoint_samples"]] <- cloglog_cutpoint_samples
  }

  # Unpack random effects
  if (model_params[["has_rfx"]]) {
    output[[
      "rfx_unique_group_ids"
    ]] <- json_object_default$get_string_vector("rfx_unique_group_ids")
    output[["rfx_samples"]] <- loadRandomEffectSamplesCombinedJson(
      json_object_list,
      0
    )
  }

  # Unpack covariate preprocessor
  if (has_field("preprocessor_metadata")) {
    preprocessor_metadata_string <- json_object_default$get_string(
      "preprocessor_metadata"
    )
    output[["train_set_metadata"]] <- createPreprocessorFromJsonString(
      preprocessor_metadata_string
    )
  } else {
    output[["train_set_metadata"]] <- NULL
    warning(paste0(
      "Field 'preprocessor_metadata' not found in JSON (model appears to have been serialized ",
      "under stochtree ",
      .ver,
      "). DataFrame covariates will not be supported for prediction. ",
      "Re-save your model to suppress this warning."
    ))
  }

  class(output) <- "bartmodel"
  return(output)
}

#' @title Convert JSON String List to Single BART Model
#' @rdname BARTSerialization
#' @param json_string_list List of JSON strings which can be parsed to objects of type `CppJson` containing Json representation of a BART model
#' @export
createBARTModelFromCombinedJsonString <- function(json_string_list) {
  # Initialize the BCF model
  output <- list()

  # Convert JSON strings
  json_object_list <- list()
  for (i in 1:length(json_string_list)) {
    json_string <- json_string_list[[i]]
    json_object_list[[i]] <- createCppJsonString(json_string)
  }

  # For scalar / preprocessing details which aren't sample-dependent,
  # defer to the first json
  json_object_default <- json_object_list[[1]]

  # Helpers for optional-field presence checks
  .ver <- inferStochtreeJsonVersion(json_object_default)
  has_field <- function(name) {
    json_contains_field_cpp(json_object_default$json_ptr, name)
  }
  has_subfolder_field <- function(subfolder, name) {
    json_contains_field_subfolder_cpp(
      json_object_default$json_ptr,
      subfolder,
      name
    )
  }

  # Unpack the forests
  include_mean_forest <- json_object_default$get_boolean(
    "include_mean_forest"
  )
  include_variance_forest <- json_object_default$get_boolean(
    "include_variance_forest"
  )
  if (include_mean_forest) {
    output[["mean_forests"]] <- loadForestContainerCombinedJson(
      json_object_list,
      "forest_0"
    )
    if (include_variance_forest) {
      output[["variance_forests"]] <- loadForestContainerCombinedJson(
        json_object_list,
        "forest_1"
      )
    }
  } else {
    output[["variance_forests"]] <- loadForestContainerCombinedJson(
      json_object_list,
      "forest_0"
    )
  }

  # Unpack metadata
  train_set_metadata = list()
  train_set_metadata[["num_numeric_vars"]] <- json_object_default$get_scalar(
    "num_numeric_vars"
  )
  train_set_metadata[[
    "num_ordered_cat_vars"
  ]] <- json_object_default$get_scalar("num_ordered_cat_vars")
  train_set_metadata[[
    "num_unordered_cat_vars"
  ]] <- json_object_default$get_scalar("num_unordered_cat_vars")
  if (train_set_metadata[["num_numeric_vars"]] > 0) {
    train_set_metadata[[
      "numeric_vars"
    ]] <- json_object_default$get_string_vector("numeric_vars")
  }
  if (train_set_metadata[["num_ordered_cat_vars"]] > 0) {
    train_set_metadata[[
      "ordered_cat_vars"
    ]] <- json_object_default$get_string_vector("ordered_cat_vars")
    train_set_metadata[[
      "ordered_unique_levels"
    ]] <- json_object_default$get_string_list(
      "ordered_unique_levels",
      train_set_metadata[["ordered_cat_vars"]]
    )
  }
  if (train_set_metadata[["num_unordered_cat_vars"]] > 0) {
    train_set_metadata[[
      "unordered_cat_vars"
    ]] <- json_object_default$get_string_vector("unordered_cat_vars")
    train_set_metadata[[
      "unordered_unique_levels"
    ]] <- json_object_default$get_string_list(
      "unordered_unique_levels",
      train_set_metadata[["unordered_cat_vars"]]
    )
  }
  output[["train_set_metadata"]] <- train_set_metadata

  # Unpack model params
  model_params = list()
  model_params[["outcome_scale"]] <- json_object_default$get_scalar(
    "outcome_scale"
  )
  model_params[["outcome_mean"]] <- json_object_default$get_scalar(
    "outcome_mean"
  )
  model_params[["standardize"]] <- json_object_default$get_boolean(
    "standardize"
  )
  model_params[["sigma2_init"]] <- json_object_default$get_scalar(
    "sigma2_init"
  )
  model_params[["sample_sigma2_global"]] <- json_object_default$get_boolean(
    "sample_sigma2_global"
  )
  model_params[["sample_sigma2_leaf"]] <- json_object_default$get_boolean(
    "sample_sigma2_leaf"
  )
  model_params[["include_mean_forest"]] <- include_mean_forest
  model_params[["include_variance_forest"]] <- include_variance_forest
  model_params[["has_rfx"]] <- json_object_default$get_boolean("has_rfx")

  if (has_field("has_rfx_basis")) {
    model_params[["has_rfx_basis"]] <- json_object_default$get_boolean(
      "has_rfx_basis"
    )
    model_params[["num_rfx_basis"]] <- json_object_default$get_scalar(
      "num_rfx_basis"
    )
  } else {
    model_params[["has_rfx_basis"]] <- FALSE
    model_params[["num_rfx_basis"]] <- 1
    warning(paste0(
      "Fields 'has_rfx_basis' and 'num_rfx_basis' not found in JSON (model appears to have been ",
      "serialized under stochtree ",
      .ver,
      "). Defaulting to FALSE / 1. ",
      "Re-save your model to suppress this warning."
    ))
  }

  model_params[["num_covariates"]] <- if (has_field("num_covariates")) {
    json_object_default$get_scalar("num_covariates")
  } else {
    NA_real_
  }
  model_params[["num_basis"]] <- json_object_default$get_scalar("num_basis")
  model_params <- .reconstructBartLeafModelFields(model_params)
  model_params[["requires_basis"]] <- json_object_default$get_boolean(
    "requires_basis"
  )

  model_params[["probit_outcome_model"]] <- if (
    has_field("probit_outcome_model")
  ) {
    json_object_default$get_boolean("probit_outcome_model")
  } else {
    FALSE
  }

  if (
    has_subfolder_field("outcome_model", "outcome") &&
      has_subfolder_field("outcome_model", "link")
  ) {
    outcome_model_outcome <- json_object_default$get_string(
      "outcome",
      "outcome_model"
    )
    outcome_model_link <- json_object_default$get_string(
      "link",
      "outcome_model"
    )
  } else {
    outcome_model_outcome <- "continuous"
    outcome_model_link <- "identity"
    warning(paste0(
      "Fields 'outcome' and 'link' not found under 'outcome_model' in JSON (model appears to have ",
      "been serialized under stochtree ",
      .ver,
      "). Defaulting to outcome='continuous', ",
      "link='identity'. Re-save your model to suppress this warning."
    ))
  }
  model_params[["outcome_model"]] <- OutcomeModel(
    outcome = outcome_model_outcome,
    link = outcome_model_link,
    family = .readOutcomeFamilyFromJson(json_object_default, outcome_model_outcome)
  )

  if (has_field("rfx_model_spec")) {
    model_params[["rfx_model_spec"]] <- json_object_default$get_string(
      "rfx_model_spec"
    )
  } else {
    model_params[["rfx_model_spec"]] <- ""
    if (model_params[["has_rfx"]]) {
      warning(paste0(
        "Field 'rfx_model_spec' not found in JSON (model appears to have been serialized under ",
        "stochtree ",
        .ver,
        "). Defaulting to ''. Re-save your model to suppress this warning."
      ))
    }
  }
  if (has_field("rfx_unseen_groups")) {
    model_params[["rfx_unseen_groups"]] <- json_object_default$get_string(
      "rfx_unseen_groups"
    )
  } else {
    # Models serialized before this option always refused unseen group labels.
    model_params[["rfx_unseen_groups"]] <- "error"
  }

  if (has_field("num_chains")) {
    model_params[["num_chains"]] <- json_object_default$get_scalar("num_chains")
  } else {
    model_params[["num_chains"]] <- 1
    warning(paste0(
      "Field 'num_chains' not found in JSON (model appears to have been serialized under stochtree ",
      .ver,
      "). Defaulting to 1. Re-save your model to suppress this warning."
    ))
  }

  if (has_field("keep_every")) {
    model_params[["keep_every"]] <- json_object_default$get_scalar("keep_every")
  } else {
    model_params[["keep_every"]] <- 1
    warning(paste0(
      "Field 'keep_every' not found in JSON (model appears to have been serialized under stochtree ",
      .ver,
      "). Defaulting to 1. Re-save your model to suppress this warning."
    ))
  }

  if (model_params[["outcome_model"]]$link == "cloglog") {
    cloglog_num_categories <- json_object_default$get_scalar(
      "cloglog_num_categories"
    )
    model_params[["cloglog_num_categories"]] <- cloglog_num_categories
  } else {
    model_params[["cloglog_num_categories"]] <- 0
  }

  # Combine values that are sample-specific
  for (i in 1:length(json_object_list)) {
    json_object <- json_object_list[[i]]
    if (i == 1) {
      model_params[["num_gfr"]] <- json_object$get_scalar("num_gfr")
      model_params[["num_burnin"]] <- json_object$get_scalar("num_burnin")
      model_params[["num_mcmc"]] <- json_object$get_scalar("num_mcmc")
      model_params[["num_samples"]] <- json_object$get_scalar(
        "num_samples"
      )
    } else {
      prev_json <- json_object_list[[i - 1]]
      model_params[["num_gfr"]] <- model_params[["num_gfr"]] +
        json_object$get_scalar("num_gfr")
      model_params[["num_burnin"]] <- model_params[["num_burnin"]] +
        json_object$get_scalar("num_burnin")
      model_params[["num_mcmc"]] <- model_params[["num_mcmc"]] +
        json_object$get_scalar("num_mcmc")
      model_params[["num_samples"]] <- model_params[["num_samples"]] +
        json_object$get_scalar("num_samples")
    }
  }
  if (identical(model_params[["outcome_model"]]$outcome, "count")) {
    negbin_fields <- .loadBartNegbinFieldsFromJson(
      json_object_list,
      model_params,
      output
    )
    model_params <- negbin_fields$model_params
    output <- negbin_fields$output
    if (identical(model_params[["outcome_model"]]$family, "zinb")) {
      zinb_fields <- .loadBartZinbFieldsFromJson(
        json_object_list,
        model_params,
        output
      )
      model_params <- zinb_fields$model_params
      output <- zinb_fields$output
      output[["zero_forests"]] <- loadForestContainerCombinedJson(
        json_object_list,
        .zeroForestJsonLabel(include_mean_forest, include_variance_forest)
      )
    }
  }
  output[["model_params"]] <- model_params

  # Unpack sampled parameters
  if (model_params[["sample_sigma2_global"]]) {
    for (i in 1:length(json_object_list)) {
      json_object <- json_object_list[[i]]
      if (i == 1) {
        output[["sigma2_global_samples"]] <- json_object$get_vector(
          "sigma2_global_samples",
          "parameters"
        )
      } else {
        output[["sigma2_global_samples"]] <- c(
          output[["sigma2_global_samples"]],
          json_object$get_vector(
            "sigma2_global_samples",
            "parameters"
          )
        )
      }
    }
  }
  if (model_params[["sample_sigma2_leaf"]]) {
    for (i in 1:length(json_object_list)) {
      json_object <- json_object_list[[i]]
      if (i == 1) {
        output[["sigma2_leaf_samples"]] <- json_object$get_vector(
          "sigma2_leaf_samples",
          "parameters"
        )
      } else {
        output[["sigma2_leaf_samples"]] <- c(
          output[["sigma2_leaf_samples"]],
          json_object$get_vector("sigma2_leaf_samples", "parameters")
        )
      }
    }
  }
  if (
    model_params[["outcome_model"]]$link == "cloglog" &&
      model_params[["outcome_model"]]$outcome == "ordinal"
  ) {
    cloglog_cutpoint_samples <- matrix(
      NA_real_,
      model_params[["cloglog_num_categories"]] - 1,
      model_params[["num_samples"]]
    )
    index_start <- 1
    for (i in 1:length(json_object_list)) {
      json_object <- json_object_list[[i]]
      num_samples <- json_object$get_scalar("num_samples")
      subset_inds <- index_start:(index_start + num_samples - 1)
      for (j in 1:(model_params[["cloglog_num_categories"]] - 1)) {
        cloglog_cutpoint_samples[j, subset_inds] <- json_object$get_vector(
          paste0("cloglog_cutpoint_samples_", j),
          "parameters"
        )
      }
    }
    output[["cloglog_cutpoint_samples"]] <- cloglog_cutpoint_samples
  }

  # Unpack random effects
  if (model_params[["has_rfx"]]) {
    output[[
      "rfx_unique_group_ids"
    ]] <- json_object_default$get_string_vector("rfx_unique_group_ids")
    output[["rfx_samples"]] <- loadRandomEffectSamplesCombinedJson(
      json_object_list,
      0
    )
  }

  # Unpack covariate preprocessor
  if (has_field("preprocessor_metadata")) {
    preprocessor_metadata_string <- json_object_default$get_string(
      "preprocessor_metadata"
    )
    output[["train_set_metadata"]] <- createPreprocessorFromJsonString(
      preprocessor_metadata_string
    )
  } else {
    output[["train_set_metadata"]] <- NULL
    warning(paste0(
      "Field 'preprocessor_metadata' not found in JSON (model appears to have been serialized ",
      "under stochtree ",
      .ver,
      "). DataFrame covariates will not be supported for prediction. ",
      "Re-save your model to suppress this warning."
    ))
  }

  class(output) <- "bartmodel"
  return(output)
}
