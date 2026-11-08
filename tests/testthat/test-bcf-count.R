# Count outcomes (negative binomial and zero-inflated negative binomial) in
# bcf(). See R/bcf_count.R and the "Count outcomes" section of ?bcf.

count_model <- function(family = "negbin") {
  OutcomeModel(outcome = "count", link = "log", family = family)
}

# Confounded binary treatment with a known propensity. Count part:
# log mu = offset + mu_c(x) + tau_c(x) z with tau_c(x) = 0.3 + 0.4 x2, a log
# rate ratio with mean 0.5; x1 drives both the propensity and mu_c. Zero part
# (zinb): pi = Phi(g0(x) + tau_z z) with tau_z = -0.5, g0 also rising in x1.
simulate_bcf_count <- function(seed, r = 3, zinb = FALSE, n = 3000, with_offset = TRUE,
                               with_rfx = FALSE) {
  set.seed(seed)
  X <- matrix(runif(n * 5), ncol = 5)
  propensity <- 0.15 + 0.7 * plogis(3 * (X[, 1] - 0.5))
  Z <- rbinom(n, 1, propensity)
  mu_c <- 0.6 + 1.2 * X[, 1] + 0.5 * sin(2 * pi * X[, 3])
  tau_c <- 0.3 + 0.4 * X[, 2]
  offset <- if (with_offset) log(runif(n, 0.5, 2)) else rep(0, n)
  groups <- sample(1:20, n, replace = TRUE)
  u <- if (with_rfx) rnorm(20, 0, 0.4) else rep(0, 20)
  log_mu <- offset + mu_c + tau_c * Z + u[groups]
  g0 <- -0.6 + 0.8 * X[, 1] - 0.6 * X[, 4]
  tau_z <- -0.5
  pi <- if (zinb) pnorm(g0 + tau_z * Z) else rep(0, n)
  s <- rbinom(n, 1, pi)
  y <- ifelse(s == 1, 0, rnbinom(n, mu = exp(log_mu), size = r))
  list(X = X, Z = Z, y = y, propensity = propensity, offset = offset, groups = groups,
    u = u, tau_c = tau_c, log_mu = log_mu, tau_z = tau_z, pi = pi,
    pi0 = pnorm(g0), pi1 = pnorm(g0 + tau_z))
}

fit_bcf_count <- function(d, family, seed) {
  bcf(d$X, d$Z, d$y, propensity_train = d$propensity, offset_train = d$offset,
    num_gfr = 10, num_burnin = 150, num_mcmc = 250,
    general_params = list(outcome_model = count_model(family), random_seed = seed))
}

small_bcf_count <- function(seed = 4, n = 300, family = "zinb") {
  d <- simulate_bcf_count(seed, zinb = family == "zinb", n = n, with_rfx = TRUE)
  te <- 1:40
  d$X_test <- d$X[te, ]
  d$Z_test <- d$Z[te]
  d$propensity_test <- d$propensity[te]
  d$offset_test <- d$offset[te]
  d$groups_test <- d$groups[te]
  d
}

fit_small_bcf_count <- function(d, family = "zinb", seed = 1, num_chains = 1,
                                num_gfr = 4, general = list(), ...) {
  bcf(d$X, d$Z, d$y, propensity_train = d$propensity, offset_train = d$offset,
    X_test = d$X_test, Z_test = d$Z_test, propensity_test = d$propensity_test,
    offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = num_gfr, num_burnin = 5, num_mcmc = 10,
    general_params = c(list(outcome_model = count_model(family), random_seed = seed,
      num_chains = num_chains), general),
    prognostic_forest_params = list(num_trees = 40),
    treatment_effect_forest_params = list(num_trees = 20),
    ...)
}

zinb_small_args <- list(
  zero_prognostic_forest_params = list(num_trees = 15),
  zero_treatment_effect_forest_params = list(num_trees = 10)
)

# ---------------------------------------------------------------------------
# Recovery of the log rate ratio and of the zero-part effect
#
# Thresholds come from a calibration run of 8 independent data sets per family
# (simulation seeds 2001-2008, fit seeds 1-8, the settings of fit_bcf_count();
# n = 3,000, r = 3, an offset, 8-10% zeros for "negbin", 29-32% for "zinb").
# The true sample-average log rate ratio is 0.495-0.503. Over those 16 fits:
#   - error of the posterior mean average log rate ratio: max 0.055 (negbin),
#     0.060 (zinb); truth inside the 95% interval of the average-CATE draws in
#     7 of 8 and 8 of 8 (the miss: interval 0.504 to 0.618 for truth 0.503)
#   - CATE RMSE: max 0.098 and 0.104; CATE correlation with tau_c(x) =
#     0.3 + 0.4 x2 (sd 0.115): min 0.74 and 0.46
#   - correlation of posterior mean rates with the truth: min 0.985 and 0.980
#   - true r = 3 inside the 99% interval of its draws: 16 of 16
#   - zinb: posterior mean of the average zero-part effect tau_z (truth -0.5):
#     max error 0.063, truth inside the 95% interval 8 of 8; average effect on
#     the zero probability, mean of pi(x, 1) - pi(x, 0) (truth -0.147): max
#     error 0.020, covered 8 of 8; pi correlation min 0.816, pi RMSE max 0.066.
# The checks use an ATE error ceiling of 0.1 plus the 95% interval, a CATE RMSE
# ceiling of 0.13, a CATE correlation floor of 0.4, a rate correlation floor of
# 0.97, the 99% interval for r and, for zinb, a tau_z error ceiling of 0.15 plus
# its 95% interval, a zero-probability-effect error ceiling of 0.04, a pi
# correlation floor of 0.78 and a pi RMSE ceiling of 0.08. Each test uses the
# first calibration seed.
#
# What the checks catch: the naive log rate ratio between the arms (exposure-
# adjusted) is 0.68-0.77 for negbin and 0.77-0.88 for zinb on these data sets,
# because assignment is confounded with x1; an estimator that ignored the
# prognostic adjustment would fail the ATE checks.
# ---------------------------------------------------------------------------

bcf_count_recovery_checks <- function(m, d) {
  ate_draws <- colMeans(m$tau_hat_train)
  ate_true <- mean(d$tau_c)
  expect_lt(abs(mean(ate_draws) - ate_true), 0.1)
  q <- quantile(ate_draws, c(0.025, 0.975))
  expect_true(q[[1]] <= ate_true && ate_true <= q[[2]], label = sprintf(
    "average log rate ratio %.3f inside 95%% interval [%.3f, %.3f]", ate_true, q[[1]], q[[2]]
  ))
  tau_hat <- rowMeans(m$tau_hat_train)
  expect_lt(sqrt(mean((tau_hat - d$tau_c)^2)), 0.13)
  expect_gt(cor(tau_hat, d$tau_c), 0.4)
  expect_gt(cor(rowMeans(exp(m$y_hat_train)), exp(d$log_mu)), 0.97)
  rq <- quantile(m$negbin_dispersion_samples, c(0.005, 0.995))
  expect_true(rq[[1]] <= 3 && 3 <= rq[[2]], label = sprintf(
    "true r = 3 inside 99%% interval [%.2f, %.2f]", rq[[1]], rq[[2]]
  ))
}

test_that("NB BCF recovers a log rate ratio under confounded assignment", {
  skip_on_cran()
  d <- simulate_bcf_count(2001, zinb = FALSE)
  m <- fit_bcf_count(d, "negbin", seed = 1)
  bcf_count_recovery_checks(m, d)
  expect_equal(dim(m$tau_hat_train), c(3000, 250))
})

test_that("ZINB BCF recovers the log rate ratio and the zero-part effect", {
  skip_on_cran()
  d <- simulate_bcf_count(2001, zinb = TRUE)
  m <- fit_bcf_count(d, "zinb", seed = 1)
  bcf_count_recovery_checks(m, d)
  tau_z_draws <- colMeans(m$zero_tau_hat_train)
  expect_lt(abs(mean(tau_z_draws) - d$tau_z), 0.15)
  tq <- quantile(tau_z_draws, c(0.025, 0.975))
  expect_true(tq[[1]] <= d$tau_z && d$tau_z <= tq[[2]], label = sprintf(
    "tau_z = %.2f inside 95%% interval [%.3f, %.3f]", d$tau_z, tq[[1]], tq[[2]]
  ))
  dpi_draws <- colMeans(
    pnorm(m$zero_mu_hat_train + m$zero_tau_hat_train) - pnorm(m$zero_mu_hat_train)
  )
  expect_lt(abs(mean(dpi_draws) - mean(d$pi1 - d$pi0)), 0.04)
  pi_hat <- rowMeans(pnorm(m$zero_hat_train))
  expect_gt(cor(pi_hat, d$pi), 0.78)
  expect_lt(sqrt(mean((pi_hat - d$pi)^2)), 0.08)
})

# ---------------------------------------------------------------------------
# Gaussian and probit bcf() are bitwise unchanged
# ---------------------------------------------------------------------------

test_that("Gaussian and probit BCF draws are bitwise identical to cma-zinb (242aae46)", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "bcf-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-bcf-reference.R"), local = TRUE)
  current <- suppressWarnings(bcf_reference_fits())
  for (name in setdiff(names(ref), c("platform", "generated_from"))) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
})

# ---------------------------------------------------------------------------
# Prediction, contrasts, posterior predictive
# ---------------------------------------------------------------------------

test_that("Count BCF predict() scales and terms are consistent", {
  skip_on_cran()
  d <- small_bcf_count()
  m <- do.call(fit_small_bcf_count, c(list(d), zinb_small_args))
  expect_equal(dim(m$zero_hat_train), c(300, 10))
  expect_equal(dim(m$zero_tau_hat_test), c(40, 10))
  expect_equal(m$zero_hat_train, m$zero_mu_hat_train + m$zero_tau_hat_train * d$Z)
  pred <- function(...) {
    predict(m, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, ...)
  }
  linear <- pred()
  expect_named(linear, c("y_hat", "mu_hat", "tau_hat", "prognostic_function",
    "cate", "rfx_predictions", "zero_prognostic", "zero_tau", "zero_forest"))
  expect_equal(linear$y_hat, m$y_hat_test)
  expect_equal(linear$mu_hat, m$mu_hat_test)
  expect_equal(linear$tau_hat, m$tau_hat_test)
  expect_equal(linear$zero_forest, m$zero_hat_test)
  expect_equal(linear$zero_prognostic, m$zero_mu_hat_test)
  expect_equal(linear$zero_tau, m$zero_tau_hat_test)
  expect_equal(
    linear$y_hat,
    linear$mu_hat + linear$tau_hat * d$Z_test + linear$rfx_predictions + d$offset_test
  )
  pi <- pnorm(linear$zero_forest)

  rate <- pred(scale = "rate")
  expect_equal(rate$y_hat, exp(linear$y_hat))
  expect_equal(rate$tau_hat, exp(linear$tau_hat))
  expect_equal(rate$cate, exp(linear$cate))
  expect_equal(rate$mu_hat, exp(linear$mu_hat))
  expect_equal(rate$rfx_predictions, exp(linear$rfx_predictions))
  expect_equal(rate$zero_forest, pi)
  expect_equal(rate$zero_tau, linear$zero_tau)

  response <- pred(scale = "response")
  expect_equal(response$y_hat, (1 - pi) * exp(linear$y_hat))
  expect_equal(response$tau_hat, rate$tau_hat)

  zero_probability <- pred(scale = "zero_probability")
  expect_named(zero_probability, c("y_hat", "zero_prognostic", "zero_tau", "zero_forest"))
  expect_equal(zero_probability$y_hat, pi)
  expect_equal(zero_probability$zero_prognostic, pnorm(linear$zero_prognostic))
  expect_equal(
    pred(terms = "y_hat", scale = "response", type = "mean"),
    rowMeans((1 - pi) * exp(linear$y_hat))
  )
  expect_warning(
    no_offset <- predict(m, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, terms = "y_hat"),
    "offset"
  )
  expect_equal(no_offset, linear$y_hat - d$offset_test)

  expect_equal(extractParameter(m, "dispersion"), m$negbin_dispersion_samples)
  expect_equal(extractParameter(m, "zero_tau_hat_train"), m$zero_tau_hat_train)
  expect_output(print(m), "zero-inflated negative binomial")
  expect_output(summary(m), "structural-zero probability")
})

test_that("Count BCF contrasts are differences of predictions on each scale", {
  skip_on_cran()
  d <- small_bcf_count()
  m <- do.call(fit_small_bcf_count, c(list(d), zinb_small_args))
  n_test <- nrow(d$X_test)
  ones <- matrix(1, n_test, 1)
  contrast <- function(scale) {
    computeContrastBCFModel(m, d$X_test, d$X_test, rep(0, n_test), rep(1, n_test),
      d$propensity_test, d$propensity_test, d$groups_test, d$groups_test, ones, ones,
      scale = scale, offset_0 = d$offset_test, offset_1 = d$offset_test)
  }
  arm <- function(z, scale) {
    predict(m, d$X_test, rep(z, n_test), d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test,
      terms = "y_hat", scale = scale)
  }
  # On the linear scale the contrast is the log rate ratio tau_0 + tau(x)
  cate <- predict(m, d$X_test, d$Z_test, d$propensity_test,
    rfx_group_ids = d$groups_test, terms = "cate")
  expect_equal(contrast("linear"), cate)
  for (scale in c("rate", "response", "zero_probability")) {
    expect_equal(contrast(scale), arm(1, scale) - arm(0, scale), label = scale)
  }
})

test_that("Count BCF posterior predictive draws are counts with the right mean and zero mass", {
  skip_on_cran()
  d <- small_bcf_count()
  m <- do.call(fit_small_bcf_count, c(list(d), zinb_small_args))
  ppd <- sampleBCFPosteriorPredictive(m, d$X_test, d$Z_test, d$propensity_test,
    rfx_group_ids = d$groups_test, offset = d$offset_test, num_draws_per_sample = 200)
  expect_equal(dim(ppd), c(40, 10, 200))
  expect_true(all(ppd >= 0 & ppd == round(ppd)))
  pred <- function(...) {
    predict(m, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, terms = "y_hat", ...)
  }
  mu <- pred(scale = "rate")
  pi <- pred(scale = "zero_probability")
  r <- matrix(m$negbin_dispersion_samples, nrow(mu), ncol(mu), byrow = TRUE)
  expect_equal(mean(ppd), mean((1 - pi) * mu), tolerance = 0.05)
  expect_equal(mean(ppd == 0), mean(pi + (1 - pi) * (r / (r + mu))^r), tolerance = 0.03)
})

# ---------------------------------------------------------------------------
# Serialization, multiple chains, warm starts, continuous dose
# ---------------------------------------------------------------------------

test_that("Count BCF models survive the JSON round trips", {
  skip_on_cran()
  d <- small_bcf_count()
  m <- do.call(fit_small_bcf_count, c(list(d, general = list(adaptive_coding = TRUE)),
    zinb_small_args))
  pred <- function(model) {
    predict(model, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, offset = d$offset_test, scale = "response")
  }
  reference <- pred(m)
  json_file <- tempfile(fileext = ".json")
  on.exit(unlink(json_file))
  saveBCFModelToJsonFile(m, json_file)
  roundtrips <- list(
    object = createBCFModelFromJson(saveBCFModelToJson(m)),
    string = createBCFModelFromJsonString(saveBCFModelToJsonString(m)),
    file = createBCFModelFromJsonFile(json_file),
    combined = createBCFModelFromCombinedJson(list(saveBCFModelToJson(m))),
    combined_string = createBCFModelFromCombinedJsonString(list(saveBCFModelToJsonString(m)))
  )
  for (name in names(roundtrips)) {
    rt <- roundtrips[[name]]
    expect_equal(pred(rt), reference, label = name)
    expect_identical(rt$model_params$outcome_model, count_model("zinb"), label = name)
    expect_equal(rt$negbin_dispersion_samples, m$negbin_dispersion_samples, label = name)
    for (field in c("zero_outcome_mean", "zinb_zero_probability_init",
      "num_trees_zero_mu", "num_trees_zero_tau", "sigma2_leaf_zero_mu",
      "sigma2_leaf_zero_tau", "outcome_mean", "has_offset",
      "log_rate_ratio_prior_sd", "negbin_dispersion_prior_rate")) {
      expect_equal(rt$model_params[[field]], m$model_params[[field]],
        label = paste(name, field))
    }
  }
  nb <- fit_small_bcf_count(small_bcf_count(family = "negbin"), family = "negbin")
  nb_rt <- createBCFModelFromJsonString(saveBCFModelToJsonString(nb))
  expect_identical(nb_rt$model_params$outcome_model, count_model("negbin"))
  expect_null(nb_rt$forests_zero_mu)
  expect_equal(nb_rt$negbin_dispersion_samples, nb$negbin_dispersion_samples)
})

test_that("Count BCF multi-chain sampling, combined JSON and warm starts have the right shapes", {
  skip_on_cran()
  d <- small_bcf_count()
  m3 <- do.call(fit_small_bcf_count, c(list(d, seed = 2, num_chains = 3, num_gfr = 6),
    zinb_small_args))
  expect_equal(dim(m3$y_hat_train), c(300, 30))
  expect_equal(dim(m3$zero_hat_test), c(40, 30))
  expect_length(m3$negbin_dispersion_samples, 30)
  expect_equal(m3$forests_zero_mu$num_samples(), 30)
  expect_equal(m3$forests_zero_tau$num_samples(), 30)
  expect_true(all(is.finite(m3$y_hat_train) & is.finite(m3$zero_hat_train)))

  m_root <- do.call(fit_small_bcf_count, c(list(d, seed = 3, num_chains = 2, num_gfr = 0),
    zinb_small_args))
  expect_equal(dim(m_root$tau_hat_train), c(300, 20))

  combined <- createBCFModelFromCombinedJsonString(list(
    saveBCFModelToJsonString(m3), saveBCFModelToJsonString(m_root)
  ))
  expect_equal(combined$model_params$num_samples, 50)
  expect_equal(
    combined$negbin_dispersion_samples,
    c(m3$negbin_dispersion_samples, m_root$negbin_dispersion_samples)
  )
  zero_pred <- predict(combined, d$X_test, d$Z_test, d$propensity_test,
    rfx_group_ids = d$groups_test, terms = "zero_forest")
  expect_equal(zero_pred, cbind(m3$zero_hat_test, m_root$zero_hat_test))

  warm <- do.call(fit_small_bcf_count, c(list(d, seed = 4, num_gfr = 0,
    previous_model_json = saveBCFModelToJsonString(m3),
    previous_model_warmstart_sample_num = 30), zinb_small_args))
  expect_equal(dim(warm$zero_hat_train), c(300, 10))
  expect_error(
    fit_small_bcf_count(d, seed = 4, num_gfr = 0,
      previous_model_json = saveBCFModelToJsonString(m3),
      previous_model_warmstart_sample_num = 30),
    "same number of trees"
  )
  expect_error(
    fit_small_bcf_count(d, family = "negbin", seed = 4, num_gfr = 0,
      previous_model_json = saveBCFModelToJsonString(m3),
      previous_model_warmstart_sample_num = 30),
    "same count family"
  )
})

test_that("NB BCF runs with a continuous dose, adaptive coding and observation weights", {
  skip_on_cran()
  set.seed(5)
  n <- 300
  X <- matrix(runif(n * 3), ncol = 3)
  dose <- runif(n)
  y <- rnbinom(n, mu = exp(0.5 + X[, 1] + 0.6 * dose), size = 4)
  m <- bcf(X, dose, y, propensity_train = rep(0.5, n), num_gfr = 3, num_mcmc = 10,
    general_params = list(outcome_model = count_model(), random_seed = 6))
  expect_equal(dim(m$tau_hat_train), c(n, 10))
  expect_true(all(is.finite(m$tau_hat_train)))
  expect_false(m$model_params$binary_treatment)

  Z <- rbinom(n, 1, 0.5)
  fit <- function(w) {
    bcf(X, Z, y, propensity_train = rep(0.5, n), observation_weights = w,
      num_gfr = 3, num_mcmc = 10,
      general_params = list(outcome_model = count_model(), random_seed = 7,
        adaptive_coding = TRUE))
  }
  unweighted <- fit(NULL)
  expect_length(unweighted$b_0_samples, 10)
  expect_true(all(is.finite(unweighted$b_1_samples)))
  expect_identical(fit(rep(1, n))$tau_hat_train, unweighted$tau_hat_train)
  expect_false(isTRUE(all.equal(fit(runif(n, 0.5, 2))$tau_hat_train, unweighted$tau_hat_train)))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("Count BCF inputs are validated", {
  skip_on_cran()
  d <- small_bcf_count(n = 80)
  cb <- function(..., y = d$y, Z = d$Z, family = "negbin", general_params = list()) {
    bcf(d$X, Z, y, propensity_train = d$propensity, num_gfr = 0, num_burnin = 0,
      num_mcmc = 2,
      general_params = c(list(outcome_model = count_model(family)), general_params), ...)
  }
  expect_error(cb(variance_forest_params = list(num_trees = 5)), "variance forest")
  expect_error(cb(general_params = list(sample_sigma2_global = TRUE)), "sample_sigma2_global")
  expect_error(cb(y = replace(d$y, 1, 2.5)), "non-negative integers")
  expect_error(cb(offset_train = d$offset[-1]), "length\\(offset_train\\)")
  expect_error(cb(offset_test = d$offset_test), "without X_test")
  expect_error(
    cb(treatment_effect_forest_params = list(log_rate_ratio_prior_sd = -1)),
    "log_rate_ratio_prior_sd"
  )
  expect_error(cb(zero_prognostic_forest_params = list(num_trees = 5)), "only used by zero-inflated")
  expect_error(cb(family = "zinb", y = d$y + 1), "at least one zero")
  expect_error(cb(family = "zinb", observation_weights = rep(1, 80)), "observation_weights")
  expect_error(
    cb(family = "zinb", Z = cbind(d$Z, runif(80))),
    "univariate"
  )
  expect_error(
    cb(family = "zinb", zero_treatment_effect_forest_params = list(num_trees = 0)),
    "num_trees >= 1"
  )

  g <- bcf(d$X, d$Z, as.numeric(d$y), propensity_train = d$propensity,
    num_gfr = 0, num_mcmc = 2)
  expect_error(
    bcf(d$X, d$Z, as.numeric(d$y), propensity_train = d$propensity,
      offset_train = d$offset, num_gfr = 0, num_mcmc = 2),
    "only supported for count"
  )
  expect_error(predict(g, d$X, d$Z, d$propensity, scale = "rate"), "count model")
  expect_error(predict(g, d$X, d$Z, d$propensity, offset = d$offset), "only supported for count")
  expect_error(
    sampleBCFPosteriorPredictive(g, d$X, d$Z, d$propensity, offset = d$offset),
    "only supported for count"
  )
  expect_error(
    computeContrastBCFModel(g, d$X, d$X, rep(0, 80), rep(1, 80), d$propensity,
      d$propensity, offset_0 = d$offset),
    "only supported for count"
  )
  nb <- cb()
  expect_error(predict(nb, d$X, d$Z, d$propensity, scale = "probability"), "count model")
  expect_error(predict(nb, d$X, d$Z, d$propensity, scale = "zero_probability"), "zero-inflated")
  zi <- cb(family = "zinb")
  expect_error(
    predict(zi, d$X, d$Z, d$propensity, scale = "zero_probability", terms = "tau"),
    "only defined for"
  )
})
