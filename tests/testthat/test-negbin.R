# Negative-binomial (count, log link) outcome of bart(), fitted by Pólya-Gamma
# augmentation. See R/negbin.R and the "Count outcomes" section of ?bart.

count_model <- function() OutcomeModel(outcome = "count", link = "log")

simulate_negbin <- function(seed, r, with_offset, with_rfx, n = 2000, p = 5) {
  set.seed(seed)
  X <- matrix(runif(n * p), ncol = p)
  f <- 0.8 + sin(2 * pi * X[, 1]) + 0.8 * X[, 2]^2
  offset <- if (with_offset) log(runif(n, 0.5, 3)) else rep(0, n)
  groups <- sample(1:20, n, replace = TRUE)
  u <- if (with_rfx) rnorm(20, 0, 0.5) else rep(0, 20)
  log_mu <- offset + f + u[groups]
  y <- rnbinom(n, mu = exp(log_mu), size = r)
  list(X = X, y = y, offset = offset, groups = groups, u = u, log_mu = log_mu)
}

fit_negbin <- function(d, with_offset, with_rfx, seed, ...) {
  args <- list(
    X_train = d$X,
    y_train = d$y,
    num_gfr = 10,
    num_burnin = 150,
    num_mcmc = 250,
    general_params = list(outcome_model = count_model(), random_seed = seed),
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

small_count_data <- function(seed = 3, n = 300, n_test = 50) {
  set.seed(seed)
  X <- matrix(runif((n + n_test) * 4), ncol = 4)
  offset <- log(runif(n + n_test, 0.5, 2))
  groups <- sample(letters[1:6], n + n_test, replace = TRUE)
  log_mu <- offset + 1 + X[, 1] - X[, 2] + c(a = -0.4, b = 0, c = 0.3, d = 0.1, e = -0.2, f = 0.5)[groups]
  y <- rnbinom(n + n_test, mu = exp(log_mu), size = 4)
  train <- seq_len(n)
  test <- n + seq_len(n_test)
  list(
    X = X[train, ], y = y[train], offset = offset[train], groups = groups[train],
    X_test = X[test, ], offset_test = offset[test], groups_test = groups[test]
  )
}

# ---------------------------------------------------------------------------
# Parameter recovery
#
# Thresholds come from a calibration run of 8 independent data sets per
# scenario (simulation seeds 1001-1008, fit seeds 1-8, the settings below;
# script and table in the PR description). Over those 32 fits:
#   - correlation of posterior mean rates with true rates: min 0.972
#   - RMSE of the log rate: max 0.174
#   - random-intercept correlation with the truth: min 0.951
#   - true r inside the 99% interval of its draws: 28 of 32 (95%: 27 of 32).
# The four misses all lie just below the interval. With mu fixed at the truth
# the r step reproduces the exact conditional posterior (mean 9.449 vs 9.453,
# sd 0.755 vs 0.746 on one data set), and the posterior of r rises with forest
# flexibility (median 8.3 with 1 tree, 9.9 with 200, oracle 9.45), so the
# upward shift is the forest absorbing part of the gamma frailty (roughly
# RMSE^2 of the log rate, 0.01 to 0.03 on the 1 / r scale), not a sampler
# error. The checks use a correlation floor of 0.95, an RMSE ceiling of 0.2,
# a random-intercept correlation floor of 0.9 and the 99% interval for r.
# Each test uses the first calibration seed of its scenario. Mutations fail
# the checks: ignoring the offset gives correlation 0.82 and RMSE 0.49 (r = 10
# scenario); replacing the Polya-Gamma weights by ones overflows the log rates
# at r = 2 and gives RMSE 0.24 with a median r of 5.9 at r = 10.
# ---------------------------------------------------------------------------

recovery_checks <- function(m, d, r) {
  rate <- rowMeans(exp(m$y_hat_train))
  expect_gt(cor(rate, exp(d$log_mu)), 0.95)
  expect_lt(sqrt(mean((log(rate) - d$log_mu)^2)), 0.2)
  q <- quantile(m$negbin_dispersion_samples, c(0.005, 0.995))
  expect_true(q[[1]] <= r && r <= q[[2]], label = sprintf(
    "true r = %g inside 99%% interval [%.2f, %.2f]", r, q[[1]], q[[2]]
  ))
}

test_that("NB BART recovers rates and dispersion: r = 2, no offset, no random effects", {
  skip_on_cran()
  d <- simulate_negbin(1001, r = 2, with_offset = FALSE, with_rfx = FALSE)
  m <- fit_negbin(d, FALSE, FALSE, seed = 1)
  recovery_checks(m, d, 2)
  expect_length(m$negbin_dispersion_samples, 250)
  expect_equal(dim(m$y_hat_train), c(2000, 250))
})

test_that("NB BART recovers rates and dispersion: r = 10 with an offset", {
  skip_on_cran()
  d <- simulate_negbin(1001, r = 10, with_offset = TRUE, with_rfx = FALSE)
  m <- fit_negbin(d, TRUE, FALSE, seed = 1)
  recovery_checks(m, d, 10)
  # The centering constant is the constant-rate MLE with exposure
  expect_equal(m$model_params$outcome_mean, log(sum(d$y) / sum(exp(d$offset))))
  # predict() with the training offset reproduces y_hat_train
  pred <- predict(m, d$X, offset = d$offset, terms = "y_hat")
  expect_equal(pred, m$y_hat_train)
})

test_that("NB BART recovers rates, dispersion and a random intercept: r = 2 with an offset", {
  skip_on_cran()
  d <- simulate_negbin(1001, r = 2, with_offset = TRUE, with_rfx = TRUE)
  m <- fit_negbin(d, TRUE, TRUE, seed = 1)
  recovery_checks(m, d, 2)
  u_hat <- rowMeans(m$rfx_samples$extract_parameter_samples()$beta_samples)
  expect_gt(cor(u_hat, d$u), 0.9)
})

test_that("NB BART recovers rates, dispersion and a random intercept: r = 10, no offset", {
  skip_on_cran()
  d <- simulate_negbin(1001, r = 10, with_offset = FALSE, with_rfx = TRUE)
  m <- fit_negbin(d, FALSE, TRUE, seed = 1)
  recovery_checks(m, d, 10)
  u_hat <- rowMeans(m$rfx_samples$extract_parameter_samples()$beta_samples)
  expect_gt(cor(u_hat, d$u), 0.9)
})

# ---------------------------------------------------------------------------
# Existing outcome paths are bitwise unchanged
# ---------------------------------------------------------------------------

test_that("Gaussian, probit and cloglog draws are bitwise identical to 0ac0e7d0", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "outcome-paths-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-outcome-reference.R"), local = TRUE)
  current <- suppressWarnings(outcome_reference_fits())
  # Weighted fits pass 1 / w where the reference passed w (observation_weights
  # became precision weights in 0.4.6.9013). 1 / (1 / w) differs from w in the
  # last bit, so those fits match to rounding, not bitwise.
  weighted <- "gaussian_weighted"
  for (name in setdiff(names(ref), c("platform", "generated_from", weighted))) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
  expect_equal(current[[weighted]], ref[[weighted]], tolerance = 1e-10,
    label = weighted)
})

test_that("Unweighted random-effects draws are bitwise identical to 0ac0e7d0", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "rfx-unweighted-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-rfx-reference.R"), local = TRUE)
  current <- rfx_reference_fits()
  for (name in c("bart_custom", "bart_intercept", "bcf_intercept")) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
})

# ---------------------------------------------------------------------------
# Prediction, posterior predictive, serialization, multiple chains
# ---------------------------------------------------------------------------

test_that("NB predict() scales, terms and offsets are consistent", {
  skip_on_cran()
  d <- small_count_data()
  m <- bart(
    d$X, d$y, X_test = d$X_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 5, num_burnin = 10, num_mcmc = 20,
    general_params = list(outcome_model = count_model(), random_seed = 9)
  )
  linear <- predict(m, d$X_test, rfx_group_ids = d$groups_test, offset = d$offset_test)
  expect_equal(linear$y_hat, m$y_hat_test)
  expect_equal(
    linear$y_hat,
    linear$mean_forest_predictions + linear$rfx_predictions + d$offset_test
  )
  rate <- predict(m, d$X_test, rfx_group_ids = d$groups_test,
    offset = d$offset_test, scale = "rate")
  expect_equal(rate$y_hat, exp(linear$y_hat))
  expect_equal(rate$mean_forest_predictions, exp(linear$mean_forest_predictions))
  expect_equal(rate$rfx_predictions, exp(linear$rfx_predictions))
  response <- predict(m, d$X_test, rfx_group_ids = d$groups_test,
    offset = d$offset_test, scale = "response")
  expect_identical(response, rate)
  rate_mean <- predict(m, d$X_test, rfx_group_ids = d$groups_test,
    offset = d$offset_test, scale = "rate", type = "mean", terms = "y_hat")
  expect_equal(rate_mean, rowMeans(exp(linear$y_hat)))

  # Missing offset: unit exposure, with a warning
  expect_warning(
    no_offset <- predict(m, d$X_test, rfx_group_ids = d$groups_test, terms = "y_hat"),
    "offset"
  )
  expect_equal(no_offset, linear$y_hat - d$offset_test)

  expect_equal(extractParameter(m, "dispersion"), m$negbin_dispersion_samples)
  expect_equal(extractParameter(m, "negbin_dispersion"), m$negbin_dispersion_samples)
  expect_output(print(m), "negative binomial")
  expect_output(summary(m), "dispersion")
})

test_that("NB posterior predictive draws are counts with the right mean", {
  skip_on_cran()
  d <- small_count_data()
  m <- bart(
    d$X, d$y, offset_train = d$offset,
    num_gfr = 5, num_burnin = 10, num_mcmc = 40,
    general_params = list(outcome_model = count_model(), random_seed = 10)
  )
  ppd <- sampleBARTPosteriorPredictive(m, X = d$X_test, offset = d$offset_test,
    num_draws_per_sample = 50)
  expect_equal(dim(ppd), c(50, 40, 50))
  expect_true(all(ppd >= 0 & ppd == round(ppd)))
  rate <- predict(m, d$X_test, offset = d$offset_test, terms = "y_hat",
    scale = "rate", type = "mean")
  # 2,000 draws per row; compare the average over rows to avoid per-row noise
  expect_equal(mean(ppd), mean(rate), tolerance = 0.05)
  ppd_default <- sampleBARTPosteriorPredictive(m, X = d$X_test, offset = d$offset_test)
  expect_equal(dim(ppd_default), c(50, 40, 25))
})

test_that("NB models survive the JSON round trips", {
  skip_on_cran()
  d <- small_count_data()
  m <- bart(
    d$X, d$y, X_test = d$X_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 5, num_burnin = 10, num_mcmc = 20,
    general_params = list(
      outcome_model = count_model(), random_seed = 11,
      negbin_dispersion_prior_shape = 2, negbin_dispersion_prior_rate = 0.05
    )
  )
  pred <- function(model) {
    predict(model, d$X_test, rfx_group_ids = d$groups_test,
      offset = d$offset_test, scale = "rate")
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
    expect_equal(rt$negbin_dispersion_samples, m$negbin_dispersion_samples, label = name)
    expect_identical(rt$model_params$outcome_model$outcome, "count")
    expect_identical(rt$model_params$outcome_model$link, "log")
    expect_true(rt$model_params$has_offset)
    expect_equal(rt$model_params$negbin_dispersion_prior_shape, 2)
    expect_equal(rt$model_params$negbin_dispersion_prior_rate, 0.05)
    expect_equal(rt$model_params$outcome_mean, m$model_params$outcome_mean)
  }
})

test_that("NB multi-chain sampling and combined JSON have the right shapes", {
  skip_on_cran()
  d <- small_count_data()
  fit <- function(seed, num_gfr, num_chains) {
    bart(
      d$X, d$y, X_test = d$X_test, offset_train = d$offset, offset_test = d$offset_test,
      num_gfr = num_gfr, num_burnin = 5, num_mcmc = 10,
      general_params = list(outcome_model = count_model(), random_seed = seed,
        num_chains = num_chains)
    )
  }
  m3 <- fit(12, num_gfr = 6, num_chains = 3)
  expect_equal(dim(m3$y_hat_train), c(300, 30))
  expect_equal(dim(m3$y_hat_test), c(50, 30))
  expect_length(m3$negbin_dispersion_samples, 30)
  expect_true(all(is.finite(m3$negbin_dispersion_samples) & m3$negbin_dispersion_samples > 0))

  m_root <- fit(13, num_gfr = 0, num_chains = 2)
  expect_equal(dim(m_root$y_hat_train), c(300, 20))
  expect_length(m_root$negbin_dispersion_samples, 20)

  combined <- createBARTModelFromCombinedJsonString(list(
    saveBARTModelToJsonString(m3), saveBARTModelToJsonString(m_root)
  ))
  expect_equal(combined$model_params$num_samples, 50)
  expect_equal(
    combined$negbin_dispersion_samples,
    c(m3$negbin_dispersion_samples, m_root$negbin_dispersion_samples)
  )
  pred <- predict(combined, d$X_test, offset = d$offset_test, terms = "y_hat")
  expect_equal(pred, cbind(m3$y_hat_test, m_root$y_hat_test))
})

test_that("NB sampler warm-starts from a previous model and honours fixed dispersion", {
  skip_on_cran()
  d <- small_count_data()
  m <- bart(d$X, d$y, num_gfr = 5, num_mcmc = 10,
    general_params = list(outcome_model = count_model(), random_seed = 14))
  m2 <- bart(d$X, d$y, num_gfr = 0, num_mcmc = 10,
    previous_model_json = saveBARTModelToJsonString(m),
    previous_model_warmstart_sample_num = 10,
    general_params = list(outcome_model = count_model(), random_seed = 15))
  expect_length(m2$negbin_dispersion_samples, 10)
  expect_error(
    bart(d$X, d$y, num_gfr = 0, num_mcmc = 5,
      previous_model_json = saveBARTModelToJsonString(m),
      previous_model_warmstart_sample_num = 10,
      general_params = list(random_seed = 15)),
    "count outcome"
  )

  fixed <- bart(d$X, d$y, num_gfr = 2, num_mcmc = 10,
    general_params = list(outcome_model = count_model(), random_seed = 16,
      sample_negbin_dispersion = FALSE, negbin_dispersion_init = 3.5))
  expect_equal(fixed$negbin_dispersion_samples, rep(3.5, 10))
  expect_false(fixed$model_params$sample_negbin_dispersion)
})

test_that("NB observation weights multiply the Polya-Gamma precision", {
  skip_on_cran()
  d <- small_count_data()
  fit <- function(w) {
    bart(d$X, d$y, observation_weights = w, num_gfr = 2, num_mcmc = 10,
      general_params = list(outcome_model = count_model(), random_seed = 17))
  }
  unweighted <- fit(NULL)
  unit <- fit(rep(1, length(d$y)))
  expect_identical(unit$y_hat_train, unweighted$y_hat_train)
  expect_identical(unit$negbin_dispersion_samples, unweighted$negbin_dispersion_samples)
  weighted <- fit(runif(length(d$y), 0.5, 2))
  expect_false(isTRUE(all.equal(weighted$y_hat_train, unweighted$y_hat_train)))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("Count outcome inputs are validated", {
  skip_on_cran()
  d <- small_count_data(n = 60, n_test = 10)
  nb <- function(..., y = d$y, general_params = list()) {
    bart(d$X, y, num_gfr = 0, num_burnin = 0, num_mcmc = 2,
      general_params = c(list(outcome_model = count_model()), general_params), ...)
  }
  expect_error(OutcomeModel("count", "identity"), "'log' for count")
  expect_error(OutcomeModel("continuous", "log"), "identity")
  expect_identical(OutcomeModel("count")$link, "log")

  expect_error(nb(y = replace(d$y, 1, -1)), "non-negative integers")
  expect_error(nb(y = replace(d$y, 1, 1.5)), "non-negative integers")
  expect_error(nb(y = replace(d$y, 1, NA)), "missing")
  expect_error(nb(y = rep(0, length(d$y))), "positive count")
  expect_error(nb(variance_forest_params = list(num_trees = 10)), "variance forest")
  expect_error(nb(general_params = list(sample_sigma2_global = TRUE)), "sample_sigma2_global")
  expect_error(nb(mean_forest_params = list(num_trees = 0)), "mean forest")
  expect_error(nb(offset_train = d$offset[-1]), "length\\(offset_train\\)")
  expect_error(nb(offset_train = replace(d$offset, 2, NA)), "finite")
  expect_error(nb(offset_train = as.character(d$offset)), "numeric")
  expect_error(nb(offset_test = d$offset_test), "without X_test")
  expect_error(nb(X_test = d$X_test, offset_test = d$offset_test[-1]), "length\\(offset_test\\)")
  expect_warning(nb(offset_train = d$offset, X_test = d$X_test), "offset_test")
  expect_error(nb(general_params = list(negbin_dispersion_prior_shape = 0)), "prior_shape")
  expect_error(nb(general_params = list(negbin_dispersion_init = -1)), "negbin_dispersion_init")
  expect_error(
    nb(general_params = list(sample_negbin_dispersion = FALSE)),
    "negbin_dispersion_init must be supplied"
  )

  # Offsets are refused for the other outcome models
  expect_error(
    bart(d$X, as.numeric(d$y), offset_train = d$offset, num_gfr = 0, num_mcmc = 2),
    "only supported for count"
  )

  m <- nb()
  expect_error(predict(m, d$X_test, scale = "probability"), "count model")
  expect_error(predict(m, d$X_test, scale = "class"), "count model")
  expect_error(predict(m, d$X_test, offset = 1:3, terms = "y_hat"), "length\\(offset\\)")
  g <- bart(d$X, as.numeric(d$y), num_gfr = 0, num_mcmc = 2)
  expect_error(predict(g, d$X_test, scale = "rate"), "count")
  expect_error(predict(g, d$X_test, offset = d$offset_test), "only supported for count")
  expect_error(sampleBARTPosteriorPredictive(g, X = d$X_test, offset = d$offset_test), "only supported for count")

  expect_error(
    bcf(d$X, Z_train = rbinom(60, 1, 0.5), y_train = d$y, propensity_train = rep(0.5, 60),
      num_gfr = 0, num_mcmc = 2, general_params = list(outcome_model = count_model()),
      variance_forest_params = list(num_trees = 5)),
    "variance forest"
  )
})
