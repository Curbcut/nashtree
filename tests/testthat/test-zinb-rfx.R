# Random intercept on the zero part of the zero-inflated negative-binomial
# family (general_params$zero_random_effects, on by default whenever
# rfx_group_ids_train is given), in bart() and bcf(). See the "Zero-inflated
# count outcomes" sections of ?bart and ?bcf.

zinb_model <- function() OutcomeModel(outcome = "count", link = "log", family = "zinb")

# Panel of G groups (tracts) with m rows each. Count part:
# log mu = 0.9 + sin(2 pi x1) + 0.6 x2^2 + u_c[g], u_c ~ N(0, 0.4^2); zero
# part: pi = Phi(-0.6 + 0.8 cos(pi x3) + 0.4 (x4 - 0.5) + u_z[g]),
# u_z ~ N(0, 0.7^2); r = 5. Mean pi about 0.3.
#
# With treatment = TRUE a binary treatment is added whose assignment depends
# on the zero-part tract intercept, P(z = 1) = Phi(a(x) + 0.9 u_z[g]) with
# a(x) = 0.6 (x1 - 0.5), and has known effects: a log rate ratio of 0.3 on the
# count part and -0.5 on the probit scale of the zero part. The propensity
# handed to bcf() is the one an analyst could estimate from x alone,
# Phi(a(x) / sqrt(1 + 0.9^2 0.7^2)), so the forests cannot recover u_z from it.
simulate_zinb_panel <- function(seed, G = 60, m = 30, treatment = FALSE,
                                sd_count = 0.4, sd_zero = 0.7, r = 5) {
  set.seed(seed)
  n <- G * m
  groups <- rep(seq_len(G), each = m)
  X <- matrix(runif(n * 4), ncol = 4)
  u_count <- rnorm(G, 0, sd_count)
  u_zero <- rnorm(G, 0, sd_zero)
  tau_count <- 0.3
  tau_zero <- -0.5
  a <- 0.6 * (X[, 1] - 0.5)
  if (treatment) {
    Z <- rbinom(n, 1, pnorm(a + 0.9 * u_zero[groups]))
    propensity <- pnorm(a / sqrt(1 + 0.9^2 * sd_zero^2))
  } else {
    Z <- rep(0, n)
    propensity <- NULL
  }
  f <- 0.9 + sin(2 * pi * X[, 1]) + 0.6 * X[, 2]^2
  g <- -0.6 + 0.8 * cos(pi * X[, 3]) + 0.4 * (X[, 4] - 0.5)
  log_mu0 <- f + u_count[groups]
  zero_index0 <- g + u_zero[groups]
  log_mu <- log_mu0 + tau_count * Z
  pi <- pnorm(zero_index0 + tau_zero * Z)
  y <- ifelse(runif(n) < pi, 0, rnbinom(n, mu = exp(log_mu), size = r))
  list(
    X = X, y = y, Z = Z, propensity = propensity, groups = groups,
    u_count = u_count, u_zero = u_zero, log_mu = log_mu, pi = pi,
    response = (1 - pi) * exp(log_mu), tau_count = tau_count, tau_zero = tau_zero,
    response0 = (1 - pnorm(zero_index0)) * exp(log_mu0),
    response1 = (1 - pnorm(zero_index0 + tau_zero)) * exp(log_mu0 + tau_count)
  )
}

fit_zinb_panel_bart <- function(d, zero_rfx, seed = 1, num_gfr = 10,
                                num_burnin = 150, num_mcmc = 250,
                                general = list(), ...) {
  bart(d$X, d$y,
    rfx_group_ids_train = d$groups,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc,
    general_params = c(list(outcome_model = zinb_model(), random_seed = seed,
      zero_random_effects = zero_rfx), general),
    mean_forest_params = list(num_trees = 100),
    ...
  )
}

fit_zinb_panel_bcf <- function(d, zero_rfx, seed = 1, num_gfr = 10,
                               num_burnin = 150, num_mcmc = 250,
                               general = list(), ...) {
  bcf(d$X, d$Z, d$y, propensity_train = d$propensity,
    rfx_group_ids_train = d$groups,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc,
    general_params = c(list(outcome_model = zinb_model(), random_seed = seed,
      zero_random_effects = zero_rfx), general),
    prognostic_forest_params = list(num_trees = 100),
    treatment_effect_forest_params = list(num_trees = 40),
    ...
  )
}

# Recovery summaries of a bart() fit on simulate_zinb_panel() data
zinb_panel_bart_summary <- function(m, d) {
  pi_hat <- rowMeans(pnorm(m$zero_hat_train))
  response_hat <- rowMeans(predict(m, d$X, rfx_group_ids = d$groups,
    terms = "y_hat", scale = "response"))
  out <- list(
    pi_cor = cor(pi_hat, d$pi),
    pi_rmse = sqrt(mean((pi_hat - d$pi)^2)),
    response_rmse = sqrt(mean((response_hat - d$response)^2))
  )
  if (isTRUE(m$model_params$has_zero_rfx)) {
    u_hat <- rowMeans(m$zero_rfx_samples$extract_parameter_samples()$beta_samples)
    out$u_zero_cor <- cor(u_hat, d$u_zero)
  }
  out
}

# Error of the posterior mean of the average zero-part effect (probit scale)
# of a bcf() fit, and of the average log rate ratio
zinb_panel_bcf_summary <- function(m, d) {
  out <- list(
    zero_effect_error = abs(mean(m$zero_tau_hat_train) - d$tau_zero),
    count_effect_error = abs(mean(m$tau_hat_train) - d$tau_count)
  )
  if (isTRUE(m$model_params$has_zero_rfx)) {
    u_hat <- rowMeans(m$zero_rfx_samples$extract_parameter_samples()$beta_samples)
    out$u_zero_cor <- cor(u_hat, d$u_zero)
  }
  out
}

# ---------------------------------------------------------------------------
# Recovery with and without the zero-part random intercept
#
# Thresholds come from a calibration run of 8 data sets (simulation seeds
# 3001-3008, fit seeds 1-8, the settings of fit_zinb_panel_*(); 60 groups x 30
# rows, 41-48% zeros). Over those 8 pairs of fits:
#   bart(): correlation of the posterior mean pi with the truth 0.882-0.927
#     with the zero-part intercept and 0.535-0.677 without (gain >= 0.236);
#     RMSE of the posterior mean response (1 - pi) mu with / without
#     <= 0.834; correlation of the zero-part group intercepts with the truth
#     0.870-0.924.
#   bcf() (treatment assignment correlated 0.35-0.49 with u_z): error of the
#     posterior mean average zero-part effect 0.076-0.296 with and
#     0.556-0.824 without (difference >= 0.340); intercept correlation
#     0.833-0.900; error of the average log rate ratio <= 0.11 either way.
# The tests below use the first data set.
# ---------------------------------------------------------------------------

test_that("The zero-part random intercept improves pi and the response in bart()", {
  skip_on_cran()
  d <- simulate_zinb_panel(3001)
  with <- zinb_panel_bart_summary(fit_zinb_panel_bart(d, TRUE, seed = 1), d)
  without <- zinb_panel_bart_summary(fit_zinb_panel_bart(d, FALSE, seed = 1), d)
  expect_gt(with$pi_cor, 0.85)
  expect_gt(with$pi_cor - without$pi_cor, 0.15)
  expect_lt(with$pi_rmse, without$pi_rmse)
  expect_lt(with$response_rmse, 0.92 * without$response_rmse)
  expect_gt(with$u_zero_cor, 0.8)
})

test_that("The zero-part random intercept removes most of the confounding of the zero-part effect in bcf()", {
  skip_on_cran()
  d <- simulate_zinb_panel(3001, treatment = TRUE)
  # Treatment assignment is correlated with the zero-part tract intercept
  expect_gt(cor(d$Z, d$u_zero[d$groups]), 0.3)
  with <- zinb_panel_bcf_summary(fit_zinb_panel_bcf(d, TRUE, seed = 1), d)
  without <- zinb_panel_bcf_summary(fit_zinb_panel_bcf(d, FALSE, seed = 1), d)
  expect_lt(with$zero_effect_error, 0.4)
  expect_lt(with$zero_effect_error, without$zero_effect_error - 0.2)
  expect_gt(with$u_zero_cor, 0.8)
  expect_lt(with$count_effect_error, 0.15)
  expect_lt(without$count_effect_error, 0.15)
})

# Small data set for the mechanics tests: 8 labelled groups, 25 rows each,
# a test set of 40 rows, an offset and a binary treatment
small_zinb_rfx_data <- function(seed = 5, G = 8, m = 25, n_test = 40) {
  set.seed(seed)
  n <- G * m + n_test
  groups <- c(rep(letters[seq_len(G)], each = m),
    sample(letters[seq_len(G)], n_test, replace = TRUE))
  X <- matrix(runif(n * 4), ncol = 4)
  offset <- log(runif(n, 0.5, 2))
  u_count <- setNames(rnorm(G, 0, 0.3), letters[seq_len(G)])
  u_zero <- setNames(rnorm(G, 0, 0.7), letters[seq_len(G)])
  propensity <- 0.3 + 0.4 * X[, 1]
  Z <- rbinom(n, 1, propensity)
  log_mu <- offset + 1 + X[, 1] - X[, 2] + 0.3 * Z + u_count[groups]
  pi <- pnorm(-0.4 + X[, 3] - 0.4 * Z + u_zero[groups])
  y <- ifelse(runif(n) < pi, 0, rnbinom(n, mu = exp(log_mu), size = 4))
  tr <- seq_len(G * m)
  te <- G * m + seq_len(n_test)
  list(
    X = X[tr, ], y = y[tr], Z = Z[tr], propensity = propensity[tr],
    offset = offset[tr], groups = groups[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = propensity[te],
    offset_test = offset[te], groups_test = groups[te]
  )
}

fit_small_bart <- function(d, seed = 1, num_chains = 1, num_gfr = 4,
                           num_mcmc = 15, general = list(), rfx = list(), ...) {
  bart(d$X, d$y, X_test = d$X_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = c(list(model_spec = "intercept_only"), rfx),
    num_gfr = num_gfr, num_burnin = 5, num_mcmc = num_mcmc,
    general_params = c(list(outcome_model = zinb_model(), random_seed = seed,
      num_chains = num_chains), general),
    mean_forest_params = list(num_trees = 30),
    zero_forest_params = list(num_trees = 15),
    ...
  )
}

fit_small_bcf <- function(d, seed = 1, num_chains = 1, num_gfr = 4,
                          num_mcmc = 15, general = list(), rfx = list(), ...) {
  bcf(d$X, d$Z, d$y, propensity_train = d$propensity,
    X_test = d$X_test, Z_test = d$Z_test, propensity_test = d$propensity_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = c(list(model_spec = "intercept_only"), rfx),
    num_gfr = num_gfr, num_burnin = 5, num_mcmc = num_mcmc,
    general_params = c(list(outcome_model = zinb_model(), random_seed = seed,
      num_chains = num_chains), general),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15),
    zero_prognostic_forest_params = list(num_trees = 15),
    zero_treatment_effect_forest_params = list(num_trees = 10),
    ...
  )
}

zero_rfx_beta <- function(m) {
  beta <- m$zero_rfx_samples$extract_parameter_samples()$beta_samples
  if (length(dim(beta)) == 3) beta[1, , ] else beta
}

# ---------------------------------------------------------------------------
# Paths without a zero-part random intercept are bitwise unchanged
# ---------------------------------------------------------------------------

test_that("ZINB without zero-part random effects and NB draws are bitwise identical to cma-bcf-count (d8a0cb7b)", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "zinb-rfx-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-zinb-rfx-reference.R"), local = TRUE)
  current <- suppressWarnings(zinb_rfx_reference_fits())
  for (name in setdiff(names(ref), c("platform", "generated_from"))) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
})

test_that("zero_random_effects is validated and only used with group ids and family = 'zinb'", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  for (bad in list("yes", NA, c(TRUE, FALSE), 1)) {
    expect_error(
      fit_small_bart(d, num_gfr = 0, num_mcmc = 2,
        general = list(zero_random_effects = bad)),
      "zero_random_effects must be TRUE or FALSE"
    )
    expect_error(
      fit_small_bcf(d, num_gfr = 0, num_mcmc = 2,
        general = list(zero_random_effects = bad)),
      "zero_random_effects must be TRUE or FALSE"
    )
  }
  on <- fit_small_bart(d, num_gfr = 0, num_mcmc = 2)
  expect_true(on$model_params$has_zero_rfx)
  off <- fit_small_bart(d, num_gfr = 0, num_mcmc = 2,
    general = list(zero_random_effects = FALSE))
  expect_false(off$model_params$has_zero_rfx)
  expect_null(off$zero_rfx_samples)
  no_groups <- bart(d$X, d$y, num_gfr = 0, num_burnin = 0, num_mcmc = 2,
    general_params = list(outcome_model = zinb_model(), random_seed = 1),
    mean_forest_params = list(num_trees = 10))
  expect_false(no_groups$model_params$has_zero_rfx)
  nb <- bart(d$X, d$y, rfx_group_ids_train = d$groups,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 0, num_burnin = 0, num_mcmc = 2,
    general_params = list(outcome_model = OutcomeModel("count", "log"),
      random_seed = 1),
    mean_forest_params = list(num_trees = 10))
  expect_null(nb$model_params$has_zero_rfx)
  expect_null(nb$zero_rfx_samples)
  bcf_on <- fit_small_bcf(d, num_gfr = 0, num_mcmc = 2)
  expect_true(bcf_on$model_params$has_zero_rfx)
})

# ---------------------------------------------------------------------------
# Prediction, posterior predictive, unseen groups
# ---------------------------------------------------------------------------

test_that("ZINB bart() predictions put the zero-part random intercept in pi", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  m <- fit_small_bart(d, seed = 2)
  expect_equal(dim(m$zero_rfx_preds_train), c(200, 15))
  expect_equal(dim(m$zero_rfx_preds_test), c(40, 15))
  # Intercept-only: every row of a group carries that group's draw
  beta <- zero_rfx_beta(m)
  expect_equal(m$zero_rfx_preds_test,
    beta[match(d$groups_test, m$rfx_unique_group_ids), ])

  pred <- function(...) {
    predict(m, d$X_test, rfx_group_ids = d$groups_test, offset = d$offset_test, ...)
  }
  linear <- pred()
  expect_named(linear, c("y_hat", "mean_forest_predictions", "rfx_predictions",
    "zero_forest_predictions", "zero_rfx_predictions"))
  expect_equal(linear$zero_rfx_predictions, m$zero_rfx_preds_test)
  # zero_hat_test is the full probit index c0 + g(x) + u_z
  expect_equal(linear$zero_forest_predictions + linear$zero_rfx_predictions,
    m$zero_hat_test)
  pi <- pnorm(m$zero_hat_test)

  zero_probability <- pred(scale = "zero_probability")
  expect_named(zero_probability, c("y_hat", "zero_forest_predictions",
    "zero_rfx_predictions"))
  expect_equal(zero_probability$y_hat, pi)
  expect_equal(zero_probability$zero_forest_predictions,
    pnorm(linear$zero_forest_predictions))
  # The zero-part term stays on the probit scale on every scale
  expect_equal(zero_probability$zero_rfx_predictions, linear$zero_rfx_predictions)

  response <- pred(scale = "response")
  expect_equal(response$y_hat, (1 - pi) * exp(linear$y_hat))
  expect_equal(pred(terms = "zero_rfx", type = "mean"),
    rowMeans(linear$zero_rfx_predictions))
  expect_equal(extractParameter(m, "zero_probability_train"), pnorm(m$zero_hat_train))

  # pi needs the group labels once the zero part has random effects
  expect_error(
    predict(m, d$X_test, terms = "y_hat", scale = "zero_probability"),
    "zero part has random effects"
  )
  # ... but the zero forest alone does not
  expect_equal(
    predict(m, d$X_test, terms = "zero_forest"),
    linear$zero_forest_predictions
  )
  expect_output(print(m), "zero-part random effects")
  expect_output(print(m), "trees and random effects")
})

test_that("ZINB bart() posterior predictive draws use pi with the zero-part intercept", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  m <- fit_small_bart(d, seed = 3, num_mcmc = 40)
  ppd <- sampleBARTPosteriorPredictive(m, X = d$X_test,
    rfx_group_ids = d$groups_test, offset = d$offset_test,
    num_draws_per_sample = 50)
  expect_equal(dim(ppd), c(40, 40, 50))
  pi <- predict(m, d$X_test, rfx_group_ids = d$groups_test, terms = "y_hat",
    scale = "zero_probability")
  mu <- predict(m, d$X_test, rfx_group_ids = d$groups_test,
    offset = d$offset_test, terms = "y_hat", scale = "rate")
  r <- matrix(m$negbin_dispersion_samples, nrow(mu), ncol(mu), byrow = TRUE)
  expect_equal(mean(ppd), mean((1 - pi) * mu), tolerance = 0.05)
  expect_equal(
    mean(ppd == 0),
    mean(pi + (1 - pi) * (r / (r + mu))^r),
    tolerance = 0.03
  )
  # Rows in groups with a high zero-part intercept get more zeros
  u <- rowMeans(m$zero_rfx_preds_test)
  hi <- u > median(u)
  p0 <- pi + (1 - pi) * (r / (r + mu))^r
  ppd_gap <- mean(ppd[hi, , ] == 0) - mean(ppd[!hi, , ] == 0)
  expect_gt(ppd_gap, 0)
  expect_lt(abs(ppd_gap - (mean(p0[hi, ]) - mean(p0[!hi, ]))), 0.03)
})

test_that("Unseen groups get the mean random effect on both parts in bart()", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  d$groups_test[1:5] <- "new"
  m <- fit_small_bart(d, seed = 4, rfx = list(unseen_groups = "mean"))
  count_beta <- m$rfx_samples$extract_parameter_samples()$beta_samples
  zero_beta <- zero_rfx_beta(m)
  for (i in 1:5) {
    expect_equal(m$rfx_preds_test[i, ], colMeans(count_beta) * m$model_params$outcome_scale)
    expect_equal(m$zero_rfx_preds_test[i, ], colMeans(zero_beta))
  }
  p <- predict(m, d$X_test, rfx_group_ids = d$groups_test, offset = d$offset_test)
  expect_equal(p$zero_rfx_predictions, m$zero_rfx_preds_test)
  expect_equal(p$zero_forest_predictions + p$zero_rfx_predictions, m$zero_hat_test)
  zp <- predict(m, d$X_test, rfx_group_ids = d$groups_test, terms = "y_hat",
    scale = "zero_probability")
  expect_equal(zp[1:5, ], pnorm(p$zero_forest_predictions[1:5, ] +
    matrix(colMeans(zero_beta), 5, ncol(zero_beta), byrow = TRUE)))

  # With the default policy, unseen labels are refused for both parts
  m_err <- fit_small_bart(small_zinb_rfx_data(), seed = 4)
  expect_error(
    predict(m_err, d$X_test, rfx_group_ids = d$groups_test, terms = "zero_rfx"),
    "must have been present"
  )
})

test_that("ZINB bcf() predictions and contrasts put the zero-part random intercept in pi", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  m <- fit_small_bcf(d, seed = 6)
  expect_equal(dim(m$zero_rfx_preds_train), c(200, 15))
  expect_equal(dim(m$zero_rfx_preds_test), c(40, 15))
  expect_equal(m$zero_rfx_preds_test,
    zero_rfx_beta(m)[match(d$groups_test, m$rfx_unique_group_ids), ])
  expect_equal(
    m$zero_hat_test,
    m$zero_mu_hat_test + m$zero_tau_hat_test * d$Z_test + m$zero_rfx_preds_test
  )
  pred <- function(...) {
    predict(m, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, ...)
  }
  zp <- pred(scale = "zero_probability")
  expect_named(zp, c("y_hat", "zero_prognostic", "zero_tau", "zero_forest",
    "zero_rfx_predictions"))
  pi <- pnorm(m$zero_hat_test)
  expect_equal(zp$y_hat, pi)
  expect_equal(zp$zero_rfx_predictions, m$zero_rfx_preds_test)
  rate <- pred(terms = "y_hat", scale = "rate")
  expect_equal(pred(terms = "y_hat", scale = "response"), (1 - pi) * rate)
  expect_equal(pred(terms = "zero_rfx", scale = "rate", type = "mean"),
    rowMeans(m$zero_rfx_preds_test))

  # Response-scale contrast: the group intercept enters pi in both arms
  n_test <- nrow(d$X_test)
  contrast <- computeContrastBCFModel(m, d$X_test, d$X_test, rep(0, n_test),
    rep(1, n_test), d$propensity_test, d$propensity_test, d$groups_test,
    d$groups_test, matrix(1, n_test, 1), matrix(1, n_test, 1),
    scale = "response", offset_0 = d$offset_test, offset_1 = d$offset_test)
  arm <- function(z) {
    predict(m, d$X_test, rep(z, n_test), d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, terms = "y_hat",
      scale = "response")
  }
  expect_equal(contrast, arm(1) - arm(0))
  pi0 <- pnorm(m$zero_mu_hat_test + m$zero_rfx_preds_test)
  expect_equal(
    arm(0),
    (1 - pi0) * predict(m, d$X_test, rep(0, n_test), d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, terms = "y_hat",
      scale = "rate")
  )

  ppd <- sampleBCFPosteriorPredictive(m, d$X_test, d$Z_test, d$propensity_test,
    rfx_group_ids = d$groups_test, offset = d$offset_test,
    num_draws_per_sample = 50)
  expect_equal(dim(ppd), c(40, 15, 50))
  expect_equal(mean(ppd == 0),
    mean(pi + (1 - pi) * (m$negbin_dispersion_samples[col(pi)] /
      (m$negbin_dispersion_samples[col(pi)] + rate))^m$negbin_dispersion_samples[col(pi)]),
    tolerance = 0.05)

  # bcf() has no unseen-group policy: unseen labels are refused, as for the
  # count part
  expect_error(
    predict(m, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = replace(d$groups_test, 1, "new"), terms = "zero_rfx"),
    "must have been present"
  )
  expect_output(print(m), "zero-part random effects")
})

# ---------------------------------------------------------------------------
# Serialization, multiple chains, warm starts
# ---------------------------------------------------------------------------

test_that("ZINB bart() with zero-part random effects survives JSON round trips", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  m <- fit_small_bart(d, seed = 7, num_chains = 2)
  all_terms <- function(model) {
    predict(model, d$X_test, rfx_group_ids = d$groups_test, offset = d$offset_test)
  }
  pi_of <- function(model) {
    predict(model, d$X_test, rfx_group_ids = d$groups_test, terms = "y_hat",
      scale = "zero_probability")
  }
  json_file <- tempfile(fileext = ".json")
  on.exit(unlink(json_file))
  saveBARTModelToJsonFile(m, json_file)
  round_trips <- list(
    string = createBARTModelFromJsonString(saveBARTModelToJsonString(m)),
    file = createBARTModelFromJsonFile(json_file),
    object = createBARTModelFromJson(saveBARTModelToJson(m)),
    combined = createBARTModelFromCombinedJson(list(saveBARTModelToJson(m))),
    combined_string = createBARTModelFromCombinedJsonString(
      list(saveBARTModelToJsonString(m))
    )
  )
  for (name in names(round_trips)) {
    rt <- round_trips[[name]]
    expect_true(rt$model_params$has_zero_rfx, label = name)
    expect_identical(
      rt$zero_rfx_samples$extract_parameter_samples(),
      m$zero_rfx_samples$extract_parameter_samples(),
      label = name
    )
    expect_identical(
      rt$rfx_samples$extract_parameter_samples(),
      m$rfx_samples$extract_parameter_samples(),
      label = name
    )
    expect_identical(all_terms(rt), all_terms(m), label = name)
    expect_identical(pi_of(rt), pi_of(m), label = name)
  }
  # The zero part's draws are not the count part's
  expect_false(isTRUE(all.equal(
    m$zero_rfx_samples$extract_parameter_samples()$beta_samples,
    m$rfx_samples$extract_parameter_samples()$beta_samples
  )))
  # Models saved without the field load as having no zero-part random effects
  off <- fit_small_bart(d, seed = 7, general = list(zero_random_effects = FALSE))
  js <- saveBARTModelToJson(off)
  expect_false(createBARTModelFromJson(js)$model_params$has_zero_rfx)
})

test_that("ZINB bcf() with zero-part random effects survives JSON round trips", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  m <- fit_small_bcf(d, seed = 8, num_chains = 2)
  all_terms <- function(model, scale) {
    predict(model, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, scale = scale)
  }
  json_file <- tempfile(fileext = ".json")
  on.exit(unlink(json_file))
  saveBCFModelToJsonFile(m, json_file)
  round_trips <- list(
    string = createBCFModelFromJsonString(saveBCFModelToJsonString(m)),
    file = createBCFModelFromJsonFile(json_file),
    object = createBCFModelFromJson(saveBCFModelToJson(m)),
    combined = createBCFModelFromCombinedJson(list(saveBCFModelToJson(m))),
    combined_string = createBCFModelFromCombinedJsonString(
      list(saveBCFModelToJsonString(m))
    )
  )
  for (name in names(round_trips)) {
    rt <- round_trips[[name]]
    expect_true(rt$model_params$has_zero_rfx, label = name)
    expect_identical(
      rt$zero_rfx_samples$extract_parameter_samples(),
      m$zero_rfx_samples$extract_parameter_samples(),
      label = name
    )
    expect_identical(all_terms(rt, "linear"), all_terms(m, "linear"), label = name)
    expect_identical(all_terms(rt, "zero_probability"),
      all_terms(m, "zero_probability"), label = name)
  }
})

test_that("Zero-part random effects have the right shapes across chains, combined models and warm starts", {
  skip_on_cran()
  d <- small_zinb_rfx_data()
  m3 <- fit_small_bart(d, seed = 9, num_chains = 3, num_gfr = 6, num_mcmc = 10)
  expect_equal(m3$model_params$num_samples, 30)
  expect_equal(dim(zero_rfx_beta(m3)), c(8, 30))
  expect_equal(dim(m3$zero_rfx_preds_train), c(200, 30))
  expect_equal(dim(m3$zero_hat_train), c(200, 30))
  m_root <- fit_small_bart(d, seed = 10, num_chains = 2, num_gfr = 0, num_mcmc = 10)
  expect_equal(dim(zero_rfx_beta(m_root)), c(8, 20))

  combined <- createBARTModelFromCombinedJsonString(list(
    saveBARTModelToJsonString(m3), saveBARTModelToJsonString(m_root)
  ))
  expect_equal(combined$model_params$num_samples, 50)
  expect_equal(zero_rfx_beta(combined), cbind(zero_rfx_beta(m3), zero_rfx_beta(m_root)))
  expect_equal(
    predict(combined, d$X_test, rfx_group_ids = d$groups_test, terms = "zero_rfx"),
    cbind(m3$zero_rfx_preds_test, m_root$zero_rfx_preds_test)
  )
  m_off <- fit_small_bart(d, seed = 11, general = list(zero_random_effects = FALSE))
  expect_error(
    createBARTModelFromCombinedJsonString(list(
      saveBARTModelToJsonString(m3), saveBARTModelToJsonString(m_off)
    )),
    "with and without zero-part random effects"
  )

  warm <- fit_small_bart(d, seed = 12, num_gfr = 0, num_mcmc = 5,
    previous_model_json = saveBARTModelToJsonString(m3),
    previous_model_warmstart_sample_num = 30)
  expect_equal(dim(zero_rfx_beta(warm)), c(8, 5))
  expect_warning(
    fit_small_bart(d, seed = 12, num_gfr = 0, num_mcmc = 5,
      previous_model_json = saveBARTModelToJsonString(m_off),
      previous_model_warmstart_sample_num = 10),
    "zero-part random effects samples"
  )

  b3 <- fit_small_bcf(d, seed = 13, num_chains = 3, num_gfr = 6, num_mcmc = 10)
  expect_equal(dim(zero_rfx_beta(b3)), c(8, 30))
  expect_equal(dim(b3$zero_rfx_preds_train), c(200, 30))
  b_root <- fit_small_bcf(d, seed = 14, num_chains = 2, num_gfr = 0, num_mcmc = 10)
  b_combined <- createBCFModelFromCombinedJsonString(list(
    saveBCFModelToJsonString(b3), saveBCFModelToJsonString(b_root)
  ))
  expect_equal(zero_rfx_beta(b_combined), cbind(zero_rfx_beta(b3), zero_rfx_beta(b_root)))
  b_warm <- fit_small_bcf(d, seed = 15, num_gfr = 0, num_mcmc = 5,
    previous_model_json = saveBCFModelToJsonString(b3),
    previous_model_warmstart_sample_num = 30)
  expect_equal(dim(zero_rfx_beta(b_warm)), c(8, 5))
})
