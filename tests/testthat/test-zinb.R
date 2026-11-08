# Zero-inflated negative-binomial (count, log link, family = "zinb") outcome
# of bart(). See R/negbin.R and the "Zero-inflated count outcomes" section of
# ?bart.

zinb_model <- function() OutcomeModel(outcome = "count", link = "log", family = "zinb")

# Count part as in test-negbin.R; zero part pi = Phi(g) with g smooth in X3 and
# X4, mean pi about 0.38.
simulate_zinb <- function(seed, r, with_offset, with_rfx, n = 2000, p = 5) {
  set.seed(seed)
  X <- matrix(runif(n * p), ncol = p)
  f <- 0.8 + sin(2 * pi * X[, 1]) + 0.8 * X[, 2]^2
  g <- -0.35 + 0.9 * cos(pi * X[, 3]) + 0.6 * (X[, 4] - 0.5)
  offset <- if (with_offset) log(runif(n, 0.5, 3)) else rep(0, n)
  groups <- sample(1:20, n, replace = TRUE)
  u <- if (with_rfx) rnorm(20, 0, 0.5) else rep(0, 20)
  log_mu <- offset + f + u[groups]
  pi <- pnorm(g)
  s <- rbinom(n, 1, pi)
  y <- ifelse(s == 1, 0, rnbinom(n, mu = exp(log_mu), size = r))
  list(X = X, y = y, offset = offset, groups = groups, u = u, log_mu = log_mu, pi = pi)
}

fit_zinb <- function(d, with_offset, with_rfx, seed, ...) {
  args <- list(
    X_train = d$X,
    y_train = d$y,
    num_gfr = 10,
    num_burnin = 150,
    num_mcmc = 250,
    general_params = list(outcome_model = zinb_model(), random_seed = seed),
    mean_forest_params = list(num_trees = 100),
    ...
  )
  if (with_offset) args$offset_train <- d$offset
  if (with_rfx) {
    args$rfx_group_ids_train <- d$groups
    args$random_effects_params <- list(model_spec = "intercept_only")
  }
  do.call(bart, args)
}

small_zinb_data <- function(seed = 3, n = 300, n_test = 50) {
  set.seed(seed)
  X <- matrix(runif((n + n_test) * 4), ncol = 4)
  offset <- log(runif(n + n_test, 0.5, 2))
  groups <- sample(letters[1:6], n + n_test, replace = TRUE)
  log_mu <- offset + 1 + X[, 1] - X[, 2] +
    c(a = -0.4, b = 0, c = 0.3, d = 0.1, e = -0.2, f = 0.5)[groups]
  pi <- pnorm(-0.4 + X[, 3])
  y <- ifelse(runif(n + n_test) < pi, 0, rnbinom(n + n_test, mu = exp(log_mu), size = 4))
  train <- seq_len(n)
  test <- n + seq_len(n_test)
  list(
    X = X[train, ], y = y[train], offset = offset[train], groups = groups[train],
    X_test = X[test, ], offset_test = offset[test], groups_test = groups[test]
  )
}

fit_small_zinb <- function(d, seed, num_chains = 1, num_gfr = 5, ...) {
  bart(
    d$X, d$y, X_test = d$X_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = num_gfr, num_burnin = 10, num_mcmc = 20,
    general_params = list(outcome_model = zinb_model(), random_seed = seed,
      num_chains = num_chains),
    zero_forest_params = list(num_trees = 20, ...)
  )
}

# ---------------------------------------------------------------------------
# Parameter recovery
#
# Thresholds come from a calibration run of 8 independent data sets per
# scenario (simulation seeds 1001-1008, fit seeds 1-8, the settings of
# fit_zinb(); n = 2,000, mean pi 0.38-0.40, 49-54% zeros at r = 2 and 43-51% at
# r = 10). Over those 32 fits:
#   - correlation of posterior mean rates mu with the truth: min 0.952
#   - RMSE of the log rate: max 0.219
#   - correlation of posterior mean zero probabilities pi with the truth:
#     min 0.934; RMSE of pi: max 0.083
#   - RMSE of the posterior mean response (1 - pi) mu, divided by the sd of the
#     true response: max 0.365
#   - random-intercept correlation with the truth: min 0.954
#   - true r inside the 99% interval of its draws: 31 of 32 (95%: 27 of 32).
# The miss is r = 10 without offset, with random effects (seed 1004): interval
# 11.2 to 47.0. At r = 10 the posterior of r is wide and shifted up (median of
# medians 11.3 and 12.5 in the two r = 10 scenarios), as for the plain NB model
# (test-negbin.R), where the forest absorbs part of the gamma frailty.
# The checks use a rate correlation floor of 0.94, a log-rate RMSE ceiling of
# 0.25, a pi correlation floor of 0.92, a pi RMSE ceiling of 0.1, a relative
# response RMSE ceiling of 0.4, a random-intercept correlation floor of 0.9 and
# the 99% interval for r. Each test uses the first calibration seed of its
# scenario.
#
# What the checks catch: a plain NB fit (family = "negbin") on the same data,
# i.e. a sampler that ignores the zero part, gives rate correlations of 0.82 to
# 0.85, log-rate RMSEs of 0.70 to 0.76 and a median r of 0.40 to 0.61 on the
# four first-seed data sets. Its response-scale error (0.20 to 0.25 relative)
# is no worse than the ZINB fit's, because an NB mean forest fits E[y] =
# (1 - pi) mu directly; the response check is a sanity bound, the rate and pi
# checks are the ones that separate the models.
# ---------------------------------------------------------------------------

zinb_recovery_checks <- function(m, d, r) {
  rate <- rowMeans(exp(m$y_hat_train))
  expect_gt(cor(rate, exp(d$log_mu)), 0.94)
  expect_lt(sqrt(mean((log(rate) - d$log_mu)^2)), 0.25)
  pi_hat <- rowMeans(pnorm(m$zero_hat_train))
  expect_gt(cor(pi_hat, d$pi), 0.92)
  expect_lt(sqrt(mean((pi_hat - d$pi)^2)), 0.1)
  response_true <- (1 - d$pi) * exp(d$log_mu)
  response <- rowMeans(exp(m$y_hat_train) * (1 - pnorm(m$zero_hat_train)))
  expect_lt(sqrt(mean((response - response_true)^2)) / sd(response_true), 0.4)
  q <- quantile(m$negbin_dispersion_samples, c(0.005, 0.995))
  expect_true(q[[1]] <= r && r <= q[[2]], label = sprintf(
    "true r = %g inside 99%% interval [%.2f, %.2f]", r, q[[1]], q[[2]]
  ))
}

test_that("ZINB BART recovers rates, zero probabilities and dispersion: r = 2, no offset, no random effects", {
  skip_on_cran()
  d <- simulate_zinb(1001, r = 2, with_offset = FALSE, with_rfx = FALSE)
  m <- fit_zinb(d, FALSE, FALSE, seed = 1)
  zinb_recovery_checks(m, d, 2)
  expect_equal(dim(m$zero_hat_train), c(2000, 250))
  expect_length(m$negbin_dispersion_samples, 250)
})

test_that("ZINB BART recovers rates, zero probabilities and dispersion: r = 10 with an offset", {
  skip_on_cran()
  d <- simulate_zinb(1001, r = 10, with_offset = TRUE, with_rfx = FALSE)
  m <- fit_zinb(d, TRUE, FALSE, seed = 1)
  zinb_recovery_checks(m, d, 10)
  pred <- predict(m, d$X, offset = d$offset, terms = c("y_hat", "zero_forest"))
  expect_equal(pred$y_hat, m$y_hat_train)
  expect_equal(pred$zero_forest_predictions, m$zero_hat_train)
})

test_that("ZINB BART recovers a random intercept: r = 2 with an offset", {
  skip_on_cran()
  d <- simulate_zinb(1001, r = 2, with_offset = TRUE, with_rfx = TRUE)
  m <- fit_zinb(d, TRUE, TRUE, seed = 1)
  zinb_recovery_checks(m, d, 2)
  u_hat <- rowMeans(m$rfx_samples$extract_parameter_samples()$beta_samples)
  expect_gt(cor(u_hat, d$u), 0.9)
})

test_that("ZINB BART recovers a random intercept: r = 10, no offset", {
  skip_on_cran()
  d <- simulate_zinb(1001, r = 10, with_offset = FALSE, with_rfx = TRUE)
  m <- fit_zinb(d, FALSE, TRUE, seed = 1)
  zinb_recovery_checks(m, d, 10)
  u_hat <- rowMeans(m$rfx_samples$extract_parameter_samples()$beta_samples)
  expect_gt(cor(u_hat, d$u), 0.9)
})

# ---------------------------------------------------------------------------
# Structural zeros drop out of the count part exactly
#
# The count forest gives a structural zero an infinite variance weight. The
# Gaussian leaf sufficient statistics add 1 / w = 0 and y / w = 0, so such a row
# changes neither a split's marginal likelihood nor a leaf draw, and a leaf
# holding only such rows has sum_w = 0: its marginal likelihood is that of the
# prior (finite) and its value is a draw from the leaf prior N(0, tau). The
# checks below drive the forest and random-effects samplers directly.
# ---------------------------------------------------------------------------

run_weighted_forest <- function(X, y, w, gfr, num_iter = 30, seed = 1, num_trees = 10, tau = 0.1) {
  dataset <- createForestDataset(X, variance_weights = w)
  residual <- createOutcome(y)
  global_config <- createGlobalModelConfig(global_error_variance = 1)
  forest_config <- createForestModelConfig(
    feature_types = rep(0L, ncol(X)), num_trees = num_trees,
    num_features = ncol(X), num_observations = nrow(X),
    variable_weights = rep(1 / ncol(X), ncol(X)), min_samples_leaf = 1,
    max_depth = 6, leaf_model_type = 0, leaf_model_scale = matrix(tau)
  )
  forest_model <- createForestModel(dataset, forest_config, global_config)
  samples <- createForestSamples(num_trees, 1, TRUE)
  active <- createForest(num_trees, 1, TRUE)
  active$prepare_for_sampler(dataset, residual, forest_model, 0, 0)
  rng <- createCppRNG(seed)
  for (i in seq_len(num_iter)) {
    forest_model$sample_one_iteration(dataset, residual, samples, active, rng,
      forest_config, global_config, keep_forest = TRUE, gfr = gfr)
  }
  samples$predict(dataset)
}

test_that("Rows with infinite variance weight leave forest draws unchanged", {
  skip_on_cran()
  set.seed(31)
  n <- 200
  X <- matrix(runif(n * 2), ncol = 2)
  excluded <- X[, 1] > 0.6
  w <- ifelse(excluded, Inf, runif(n, 0.5, 2))
  y <- 2 * X[, 2] + rnorm(n, 0, 0.3)
  for (gfr in c(TRUE, FALSE)) {
    pred <- run_weighted_forest(X, replace(y, excluded, 0), w, gfr)
    expect_true(all(is.finite(pred)))
    # The outcomes of excluded rows never enter: changing them to huge values
    # gives the same draws bitwise.
    pred_huge <- run_weighted_forest(X, replace(y, excluded, 1e6), w, gfr)
    expect_identical(pred_huge, pred)
  }
})

test_that("A forest whose rows are all excluded draws its leaves from the prior", {
  skip_on_cran()
  set.seed(32)
  n <- 100
  X <- matrix(runif(n * 2), ncol = 2)
  tau <- 0.1
  num_trees <- 10
  for (gfr in c(TRUE, FALSE)) {
    pred <- run_weighted_forest(X, rep(1e6, n), rep(Inf, n), gfr,
      num_iter = 400, num_trees = num_trees, tau = tau)
    expect_true(all(is.finite(pred)))
    # Each prediction is a sum of num_trees independent N(0, tau) leaves
    draws <- pred[1, -(1:50)]
    expect_lt(abs(mean(draws)), 5 * sqrt(num_trees * tau / length(draws)))
    expect_gt(sd(draws), 0.6 * sqrt(num_trees * tau))
    expect_lt(sd(draws), 1.4 * sqrt(num_trees * tau))
  }
})

test_that("Rows with infinite variance weight drop out of the random-effects conditional", {
  skip_on_cran()
  set.seed(33)
  groups <- rep(1:4, each = 50)
  n <- length(groups)
  # Group 4 is entirely excluded, groups 1-3 partly
  excluded <- groups == 4 | (runif(n) < 0.3)
  w <- ifelse(excluded, Inf, 1)
  y <- c(-1, 0, 1, 0)[groups] + rnorm(n, 0, 0.2)
  run <- function(y_values) {
    dataset <- createRandomEffectsDataset(groups, matrix(1, n, 1), w)
    outcome <- createOutcome(y_values)
    tracker <- createRandomEffectsTracker(groups)
    model <- createRandomEffectsModel(1, 4)
    model$set_working_parameter(1)
    model$set_group_parameters(matrix(0, 1, 4))
    model$set_working_parameter_cov(matrix(1))
    model$set_group_parameter_cov(matrix(1))
    model$set_variance_prior_shape(1)
    model$set_variance_prior_scale(1)
    samples <- createRandomEffectSamples(1, 4, tracker)
    rng <- createCppRNG(7)
    for (i in 1:300) {
      model$sample_random_effect(dataset, outcome, tracker, samples, TRUE, 1.0, rng)
    }
    samples$extract_parameter_samples()$beta_samples
  }
  beta <- run(replace(y, excluded, 0))
  expect_true(all(is.finite(beta)))
  expect_identical(run(replace(y, excluded, 1e6)), beta)
  # The fully excluded group has no data, so its draws are prior draws (centred
  # on 0) while the others track their group means
  post_mean <- rowMeans(beta[, -(1:50)])
  expect_lt(max(abs(post_mean[1:3] - c(-1, 0, 1))), 0.15)
})

# ---------------------------------------------------------------------------
# family = "negbin" (the default) is bitwise unchanged
# ---------------------------------------------------------------------------

test_that("Negative-binomial draws are bitwise identical to cma-nb (da675755)", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "negbin-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-negbin-reference.R"), local = TRUE)
  current <- negbin_reference_fits()
  for (name in setdiff(names(ref), c("platform", "generated_from"))) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
  expect_identical(OutcomeModel("count", "log")$family, "negbin")
  expect_identical(
    OutcomeModel("count", "log"),
    OutcomeModel("count", "log", family = "negbin")
  )
})

# ---------------------------------------------------------------------------
# Prediction, posterior predictive, serialization, multiple chains
# ---------------------------------------------------------------------------

test_that("ZINB predict() scales and terms are consistent", {
  skip_on_cran()
  d <- small_zinb_data()
  m <- fit_small_zinb(d, seed = 9)
  expect_equal(dim(m$zero_hat_train), c(300, 20))
  expect_equal(dim(m$zero_hat_test), c(50, 20))
  expect_null(m$sigma2_leaf_zero_samples)
  pred <- function(...) {
    predict(m, d$X_test, rfx_group_ids = d$groups_test, offset = d$offset_test, ...)
  }
  linear <- pred()
  expect_named(linear, c("y_hat", "mean_forest_predictions", "rfx_predictions",
    "zero_forest_predictions"))
  expect_equal(linear$y_hat, m$y_hat_test)
  expect_equal(linear$zero_forest_predictions, m$zero_hat_test)
  expect_equal(
    linear$y_hat,
    linear$mean_forest_predictions + linear$rfx_predictions + d$offset_test
  )
  pi <- pnorm(linear$zero_forest_predictions)

  rate <- pred(scale = "rate")
  expect_equal(rate$y_hat, exp(linear$y_hat))
  expect_equal(rate$mean_forest_predictions, exp(linear$mean_forest_predictions))
  expect_equal(rate$rfx_predictions, exp(linear$rfx_predictions))
  expect_equal(rate$zero_forest_predictions, pi)

  response <- pred(scale = "response")
  expect_equal(response$y_hat, (1 - pi) * exp(linear$y_hat))
  expect_equal(response$mean_forest_predictions, rate$mean_forest_predictions)
  expect_equal(response$zero_forest_predictions, pi)

  zero_probability <- pred(scale = "zero_probability")
  expect_named(zero_probability, c("y_hat", "zero_forest_predictions"))
  expect_equal(zero_probability$y_hat, pi)
  expect_equal(zero_probability$zero_forest_predictions, pi)
  # No offset is needed (and no warning raised) on the zero-probability scale
  expect_silent(zp <- predict(m, d$X_test, scale = "zero_probability", terms = "y_hat"))
  expect_equal(zp, pi)

  expect_equal(pred(terms = "zero_forest"), linear$zero_forest_predictions)
  expect_equal(
    pred(terms = "y_hat", scale = "response", type = "mean"),
    rowMeans((1 - pi) * exp(linear$y_hat))
  )
  expect_equal(
    pred(terms = "zero_forest", scale = "zero_probability", type = "mean"),
    rowMeans(pi)
  )
  expect_equal(extractParameter(m, "zero_hat_train"), m$zero_hat_train)
  expect_equal(extractParameter(m, "zero_probability_test"), pnorm(m$zero_hat_test))
  expect_error(extractParameter(m, "sigma2_leaf_zero"), "leaf scale")
  expect_output(print(m), "zero-inflated negative binomial")
  expect_output(print(m), "zero forest")
  expect_output(summary(m), "structural-zero probability")
})

test_that("ZINB posterior predictive draws include the zero mass", {
  skip_on_cran()
  d <- small_zinb_data()
  m <- bart(
    d$X, d$y, offset_train = d$offset,
    num_gfr = 5, num_burnin = 10, num_mcmc = 40,
    general_params = list(outcome_model = zinb_model(), random_seed = 10),
    zero_forest_params = list(num_trees = 20)
  )
  ppd <- sampleBARTPosteriorPredictive(m, X = d$X_test, offset = d$offset_test,
    num_draws_per_sample = 50)
  expect_equal(dim(ppd), c(50, 40, 50))
  expect_true(all(ppd >= 0 & ppd == round(ppd)))
  mu <- predict(m, d$X_test, offset = d$offset_test, terms = "y_hat", scale = "rate")
  pi <- predict(m, d$X_test, terms = "y_hat", scale = "zero_probability")
  r <- matrix(m$negbin_dispersion_samples, nrow(mu), ncol(mu), byrow = TRUE)
  # 2,000 draws per row; compare averages over rows to avoid per-row noise
  expect_equal(mean(ppd), mean((1 - pi) * mu), tolerance = 0.05)
  expect_equal(
    mean(ppd == 0),
    mean(pi + (1 - pi) * (r / (r + mu))^r),
    tolerance = 0.03
  )
})

test_that("ZINB models survive the JSON round trips", {
  skip_on_cran()
  d <- small_zinb_data()
  m <- fit_small_zinb(d, seed = 11, sample_sigma2_leaf = TRUE)
  expect_length(m$sigma2_leaf_zero_samples, 20)
  pred <- function(model) {
    predict(model, d$X_test, rfx_group_ids = d$groups_test,
      offset = d$offset_test, scale = "response")
  }
  reference <- pred(m)
  json_file <- tempfile(fileext = ".json")
  on.exit(unlink(json_file))
  saveBARTModelToJsonFile(m, json_file)
  roundtrips <- list(
    object = createBARTModelFromJson(saveBARTModelToJson(m)),
    string = createBARTModelFromJsonString(saveBARTModelToJsonString(m)),
    file = createBARTModelFromJsonFile(json_file),
    combined = createBARTModelFromCombinedJson(list(saveBARTModelToJson(m))),
    combined_string = createBARTModelFromCombinedJsonString(list(saveBARTModelToJsonString(m)))
  )
  for (name in names(roundtrips)) {
    rt <- roundtrips[[name]]
    expect_equal(pred(rt), reference, label = name)
    expect_identical(rt$model_params$outcome_model, zinb_model(), label = name)
    expect_equal(rt$negbin_dispersion_samples, m$negbin_dispersion_samples, label = name)
    expect_equal(rt$sigma2_leaf_zero_samples, m$sigma2_leaf_zero_samples, label = name)
    for (field in c("zero_outcome_mean", "zinb_zero_probability_init",
      "num_trees_zero", "sample_sigma2_leaf_zero", "sigma2_leaf_zero_init",
      "outcome_mean", "has_offset")) {
      expect_equal(rt$model_params[[field]], m$model_params[[field]],
        label = paste(name, field))
    }
  }
  # A plain NB model records its family and loads back as "negbin"
  nb <- bart(d$X, d$y, num_gfr = 2, num_mcmc = 5,
    general_params = list(outcome_model = OutcomeModel("count", "log"), random_seed = 1))
  nb_rt <- createBARTModelFromJsonString(saveBARTModelToJsonString(nb))
  expect_identical(nb_rt$model_params$outcome_model$family, "negbin")
  expect_null(nb_rt$zero_forests)
})

test_that("ZINB multi-chain sampling, combined JSON and warm starts have the right shapes", {
  skip_on_cran()
  d <- small_zinb_data()
  m3 <- fit_small_zinb(d, seed = 12, num_chains = 3, num_gfr = 6, sample_sigma2_leaf = TRUE)
  expect_equal(dim(m3$y_hat_train), c(300, 60))
  expect_equal(dim(m3$zero_hat_train), c(300, 60))
  expect_equal(dim(m3$zero_hat_test), c(50, 60))
  expect_length(m3$negbin_dispersion_samples, 60)
  expect_length(m3$sigma2_leaf_zero_samples, 60)
  expect_equal(m3$zero_forests$num_samples(), 60)
  expect_true(all(is.finite(m3$zero_hat_train)))

  m_root <- fit_small_zinb(d, seed = 13, num_chains = 2, num_gfr = 0,
    sample_sigma2_leaf = TRUE)
  expect_equal(dim(m_root$zero_hat_train), c(300, 40))
  expect_equal(m_root$zero_forests$num_samples(), 40)

  combined <- createBARTModelFromCombinedJsonString(list(
    saveBARTModelToJsonString(m3), saveBARTModelToJsonString(m_root)
  ))
  expect_equal(combined$model_params$num_samples, 100)
  pred <- predict(combined, d$X_test, terms = "zero_forest")
  expect_equal(pred, cbind(m3$zero_hat_test, m_root$zero_hat_test))
  expect_equal(
    combined$sigma2_leaf_zero_samples,
    c(m3$sigma2_leaf_zero_samples, m_root$sigma2_leaf_zero_samples)
  )
  m_fixed <- fit_small_zinb(d, seed = 15, num_gfr = 2)
  expect_error(
    createBARTModelFromCombinedJsonString(list(
      saveBARTModelToJsonString(m3), saveBARTModelToJsonString(m_fixed)
    )),
    "differ in zero_forest_params\\$sample_sigma2_leaf"
  )

  warm <- bart(d$X, d$y, offset_train = d$offset, num_gfr = 0, num_mcmc = 5,
    rfx_group_ids_train = d$groups,
    random_effects_params = list(model_spec = "intercept_only"),
    previous_model_json = saveBARTModelToJsonString(m3),
    previous_model_warmstart_sample_num = 60,
    general_params = list(outcome_model = zinb_model(), random_seed = 14),
    zero_forest_params = list(num_trees = 20, sample_sigma2_leaf = TRUE))
  expect_equal(dim(warm$zero_hat_train), c(300, 5))
  expect_error(
    bart(d$X, d$y, num_gfr = 0, num_mcmc = 5,
      previous_model_json = saveBARTModelToJsonString(m3),
      previous_model_warmstart_sample_num = 60,
      general_params = list(outcome_model = zinb_model(), random_seed = 14)),
    "must match the zero forest"
  )
  expect_error(
    bart(d$X, d$y, num_gfr = 0, num_mcmc = 5,
      previous_model_json = saveBARTModelToJsonString(m3),
      previous_model_warmstart_sample_num = 60,
      general_params = list(outcome_model = OutcomeModel("count", "log"), random_seed = 14)),
    "family = 'zinb'"
  )
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("Zero-inflated count inputs are validated", {
  skip_on_cran()
  d <- small_zinb_data(n = 60, n_test = 10)
  zinb <- function(..., y = d$y, general_params = list()) {
    bart(d$X, y, num_gfr = 0, num_burnin = 0, num_mcmc = 2,
      general_params = c(list(outcome_model = zinb_model()), general_params), ...)
  }
  expect_error(OutcomeModel("count", "log", family = "poisson"), "'negbin' or 'zinb'")
  expect_error(OutcomeModel("continuous", family = "zinb"), "only available for count")
  expect_error(OutcomeModel("binary", "probit", family = "negbin"), "only available for count")
  expect_null(OutcomeModel("binary", "probit")$family)

  expect_error(zinb(y = d$y + 1), "at least one zero")
  expect_error(zinb(observation_weights = rep(1, 60)), "observation_weights")
  expect_error(zinb(zero_forest_params = list(num_trees = 0)), "zero forest")
  expect_error(zinb(zero_forest_params = list(sigma2_leaf_init = -1)), "sigma2_leaf_init")
  expect_error(
    zinb(zero_forest_params = list(keep_vars = 9)),
    "zero_forest_params\\$keep_vars"
  )
  expect_error(zinb(variance_forest_params = list(num_trees = 5)), "variance forest")
  expect_error(
    bart(d$X, d$y, num_gfr = 0, num_mcmc = 2,
      general_params = list(outcome_model = OutcomeModel("count", "log")),
      zero_forest_params = list(num_trees = 5)),
    "only used by zero-inflated"
  )
  expect_error(
    bart(d$X, as.numeric(d$y), num_gfr = 0, num_mcmc = 2,
      zero_forest_params = list(num_trees = 5)),
    "only used by zero-inflated"
  )

  m <- zinb(zero_forest_params = list(num_trees = 5, drop_vars = 1))
  expect_error(predict(m, d$X_test, scale = "zero_probability", terms = "mean_forest"),
    "only defined for")
  expect_error(predict(m, d$X_test, scale = "probability"), "count model")
  nb <- bart(d$X, d$y, num_gfr = 0, num_mcmc = 2,
    general_params = list(outcome_model = OutcomeModel("count", "log")))
  expect_error(predict(nb, d$X_test, scale = "zero_probability"), "zero-inflated")
  expect_warning(
    nb_zero <- predict(nb, d$X_test, terms = "zero_forest"),
    "None of the requested model terms"
  )
  expect_null(nb_zero)
  expect_error(
    bcf(d$X, Z_train = rbinom(60, 1, 0.5), y_train = d$y, propensity_train = rep(0.5, 60),
      num_gfr = 0, num_mcmc = 2, general_params = list(outcome_model = zinb_model()),
      variance_forest_params = list(num_trees = 5)),
    "variance forest"
  )
})
