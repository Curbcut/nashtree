# A prior of its own for the zero part's random effects in the zero-inflated
# negative-binomial family (zero_random_effects_params in bart() and bcf()).
# The default, an empty list, gives the zero part the count part's prior as
# before; any key set there replaces the count part's value for the zero part
# only.

zinb_prior_model <- function() OutcomeModel(outcome = "count", link = "log", family = "zinb")

# Six groups of very different zero propensity, so the zero-part intercepts are
# far from zero unless the prior holds them there
zero_prior_data <- function(seed = 4242, n = 240, n_test = 30) {
  set.seed(seed)
  N <- n + n_test
  X <- matrix(runif(N * 4), ncol = 4)
  groups <- sample(1:6, N, replace = TRUE)
  u_count <- c(-0.5, 0.4, 0, 0.3, -0.3, 0.6)
  u_zero <- c(-1.2, 1, 0, 0.8, -0.8, 1.3)
  propensity <- 0.3 + 0.4 * X[, 1]
  Z <- rbinom(N, 1, propensity)
  offset <- log(runif(N, 0.5, 2))
  log_mu <- offset + 1 + X[, 1] - X[, 2] + 0.3 * Z + u_count[groups]
  pi <- pnorm(-0.4 + X[, 3] - 0.3 * Z + u_zero[groups])
  y <- ifelse(runif(N) < pi, 0, rnbinom(N, mu = exp(log_mu), size = 4))
  tr <- seq_len(n)
  te <- n + seq_len(n_test)
  list(
    X = X[tr, ], y = y[tr], Z = Z[tr], propensity = propensity[tr],
    offset = offset[tr], groups = groups[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = propensity[te],
    offset_test = offset[te], groups_test = groups[te]
  )
}

# general_params of the helpers below; entries of `general` replace the
# defaults whole (c() would keep the first of two outcome_model entries)
with_general <- function(seed, num_chains, general) {
  out <- list(outcome_model = zinb_prior_model(), random_seed = seed,
    num_chains = num_chains)
  out[names(general)] <- general
  out
}

prior_fit_bart <- function(d, seed = 1, num_chains = 1, num_gfr = 4,
                           num_mcmc = 15, rfx = list(), zero_rfx = list(),
                           general = list(), ...) {
  bart(d$X, d$y, X_test = d$X_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = c(list(model_spec = "intercept_only"), rfx),
    zero_random_effects_params = zero_rfx,
    num_gfr = num_gfr, num_burnin = 5, num_mcmc = num_mcmc,
    general_params = with_general(seed, num_chains, general),
    mean_forest_params = list(num_trees = 30),
    zero_forest_params = list(num_trees = 15),
    ...
  )
}

prior_fit_bcf <- function(d, seed = 1, num_chains = 1, num_gfr = 4,
                          num_mcmc = 15, rfx = list(), zero_rfx = list(),
                          general = list(), ...) {
  bcf(d$X, d$Z, d$y, propensity_train = d$propensity,
    X_test = d$X_test, Z_test = d$Z_test, propensity_test = d$propensity_test,
    offset_train = d$offset, offset_test = d$offset_test,
    rfx_group_ids_train = d$groups, rfx_group_ids_test = d$groups_test,
    random_effects_params = c(list(model_spec = "intercept_only"), rfx),
    zero_random_effects_params = zero_rfx,
    num_gfr = num_gfr, num_burnin = 5, num_mcmc = num_mcmc,
    general_params = with_general(seed, num_chains, general),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15),
    zero_prognostic_forest_params = list(num_trees = 15),
    zero_treatment_effect_forest_params = list(num_trees = 10),
    ...
  )
}

# Everything a fit draws, for bitwise comparisons
prior_fit_draws <- function(m) {
  list(
    y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
    zero_hat_train = m$zero_hat_train, zero_hat_test = m$zero_hat_test,
    tau_hat_train = m$tau_hat_train, zero_tau_hat_train = m$zero_tau_hat_train,
    dispersion = m$negbin_dispersion_samples,
    rfx = m$rfx_samples$extract_parameter_samples(),
    zero_rfx = m$zero_rfx_samples$extract_parameter_samples()
  )
}

# Largest absolute group effect over all draws
max_abs_rfx <- function(samples) {
  max(abs(samples$extract_parameter_samples()$beta_samples))
}

# A prior that holds the effects at about zero: the working parameter's prior
# SD is 1e-5, so |alpha * xi| stays far below 0.01 for any plausible xi
tight_prior <- list(working_parameter_prior_cov = 1e-10)
# The count part's default prior written out
default_prior <- list(
  working_parameter_prior_cov = 1, group_parameter_prior_cov = 1,
  variance_prior_shape = 1, variance_prior_scale = 1
)

# ---------------------------------------------------------------------------
# Paths that do not set a zero-part prior are bitwise unchanged
# ---------------------------------------------------------------------------

test_that("Gaussian, probit, cloglog, NB and default-prior ZINB draws are bitwise identical to cma-rfx (1273f128)", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "zinb-zero-rfx-prior-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-zinb-zero-rfx-prior-reference.R"), local = TRUE)
  current <- suppressWarnings(zinb_zero_rfx_prior_reference_fits())
  for (name in setdiff(names(ref), c("platform", "generated_from"))) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
})

test_that("A zero-part prior equal to the count part's gives the default draws bitwise", {
  skip_on_cran()
  d <- zero_prior_data()
  count_prior <- list(working_parameter_prior_cov = 0.5,
    group_parameter_prior_cov = 2, variance_prior_shape = 2,
    variance_prior_scale = 0.3)
  for (fit in list(prior_fit_bart, prior_fit_bcf)) {
    # Default prior written out in full, and key by key
    base <- prior_fit_draws(fit(d, seed = 3, num_chains = 2))
    expect_identical(
      prior_fit_draws(fit(d, seed = 3, num_chains = 2, zero_rfx = default_prior)),
      base
    )
    expect_identical(
      prior_fit_draws(fit(d, seed = 3, num_chains = 2,
        zero_rfx = list(variance_prior_shape = 1))),
      base
    )
    # A non-default count prior is inherited by the zero part
    inherited <- prior_fit_draws(fit(d, seed = 4, rfx = count_prior))
    expect_identical(
      prior_fit_draws(fit(d, seed = 4, rfx = count_prior, zero_rfx = count_prior)),
      inherited
    )
    expect_false(identical(prior_fit_draws(fit(d, seed = 4)), inherited))
    # NULL entries keep the count part's value
    expect_identical(
      prior_fit_draws(fit(d, seed = 4, rfx = count_prior,
        zero_rfx = list(variance_prior_scale = NULL))),
      inherited
    )
  }
})

# ---------------------------------------------------------------------------
# A distinct zero-part prior reaches the zero part, and only the zero part
# ---------------------------------------------------------------------------

test_that("A tight zero-part prior holds the zero-part effects at zero in every chain start in bart()", {
  skip_on_cran()
  d <- zero_prior_data()
  base <- prior_fit_bart(d, seed = 5, num_chains = 2)
  expect_gt(max_abs_rfx(base$zero_rfx_samples), 0.3)

  gfr <- prior_fit_bart(d, seed = 5, num_chains = 2, zero_rfx = tight_prior)
  root <- prior_fit_bart(d, seed = 6, num_chains = 2, num_gfr = 0,
    zero_rfx = tight_prior)
  for (m in list(gfr, root)) {
    expect_lt(max_abs_rfx(m$zero_rfx_samples), 0.01)
    expect_lt(max(abs(m$zero_rfx_preds_test)), 0.01)
    # The count part keeps its own (default) prior
    expect_gt(max_abs_rfx(m$rfx_samples), 0.2)
  }

  # Warm start from a model fitted under the default prior: the chain starts
  # from that model's effects but samples under the tight prior from then on
  warm <- prior_fit_bart(d, seed = 7, num_gfr = 0, num_mcmc = 10,
    zero_rfx = tight_prior,
    previous_model_json = saveBARTModelToJsonString(base),
    previous_model_warmstart_sample_num = 30)
  beta <- warm$zero_rfx_samples$extract_parameter_samples()$beta_samples
  expect_lt(max(abs(beta[, 10])), 0.01)

  # The reverse: a tight count prior with the default restored on the zero part
  rev <- prior_fit_bart(d, seed = 5, num_chains = 2, rfx = tight_prior,
    zero_rfx = default_prior)
  expect_lt(max_abs_rfx(rev$rfx_samples), 0.01)
  expect_gt(max_abs_rfx(rev$zero_rfx_samples), 0.3)
})

test_that("A tight zero-part prior holds the zero-part effects at zero in every chain start in bcf()", {
  skip_on_cran()
  d <- zero_prior_data()
  base <- prior_fit_bcf(d, seed = 8, num_chains = 2)
  expect_gt(max_abs_rfx(base$zero_rfx_samples), 0.3)

  gfr <- prior_fit_bcf(d, seed = 8, num_chains = 2, zero_rfx = tight_prior)
  root <- prior_fit_bcf(d, seed = 9, num_chains = 2, num_gfr = 0,
    zero_rfx = tight_prior)
  for (m in list(gfr, root)) {
    expect_lt(max_abs_rfx(m$zero_rfx_samples), 0.01)
    expect_lt(max(abs(m$zero_rfx_preds_test)), 0.01)
    expect_gt(max_abs_rfx(m$rfx_samples), 0.2)
  }

  warm <- prior_fit_bcf(d, seed = 10, num_gfr = 0, num_mcmc = 10,
    zero_rfx = tight_prior,
    previous_model_json = saveBCFModelToJsonString(base),
    previous_model_warmstart_sample_num = 30)
  beta <- warm$zero_rfx_samples$extract_parameter_samples()$beta_samples
  expect_lt(max(abs(beta[, 10])), 0.01)

  rev <- prior_fit_bcf(d, seed = 8, num_chains = 2, rfx = tight_prior,
    zero_rfx = default_prior)
  expect_lt(max_abs_rfx(rev$rfx_samples), 0.01)
  expect_gt(max_abs_rfx(rev$zero_rfx_samples), 0.3)
})

test_that("A half-Cauchy zero-part prior holds weak group effects and lets strong ones through", {
  skip_on_cran()
  d <- zero_prior_data()
  # Half-Cauchy(0, A) on the zero-part intercept SD: shape = scale = 1/2 on
  # the group variance, working-parameter variance A^2. The six zero-part
  # intercepts here have SD about 1; a heavy-tailed prior lets that through
  # unless its scale is tiny. Calibration (seeds 11 and 12, 60 draws): mean
  # realized SD 0.0007-0.0013 at A = 0.001, 0.86-1.10 at A = 0.02, 1.06-1.19 at
  # A = 0.1 and 1.05-1.20 at A = 1, with the group-variance draws rising to
  # about 30 at A = 0.02 to make up for the small working parameter.
  half_cauchy <- function(A) {
    list(working_parameter_prior_cov = A^2, variance_prior_shape = 0.5,
      variance_prior_scale = 0.5)
  }
  realized_sd <- function(m) {
    mean(apply(m$zero_rfx_samples$extract_parameter_samples()$beta_samples, 2, sd))
  }
  sds <- vapply(c(0.001, 0.1), function(A) {
    realized_sd(prior_fit_bcf(d, seed = 11, num_gfr = 4, num_mcmc = 60,
      zero_rfx = half_cauchy(A)))
  }, numeric(1))
  expect_lt(sds[1], 0.01)
  expect_gt(sds[2], 0.6)
})

test_that("The zero-part prior can be read back and survives JSON round trips", {
  skip_on_cran()
  d <- zero_prior_data()
  m <- prior_fit_bart(d, seed = 12, zero_rfx = list(variance_prior_shape = 0.5,
    variance_prior_scale = 0.5, working_parameter_prior_cov = 0.01))
  m2 <- createBARTModelFromJsonString(saveBARTModelToJsonString(m))
  expect_equal(
    predict(m2, d$X_test, rfx_group_ids = d$groups_test, terms = "zero_rfx"),
    m$zero_rfx_preds_test
  )
  b <- prior_fit_bcf(d, seed = 13, zero_rfx = list(working_parameter_prior_cov = 0.01))
  b2 <- createBCFModelFromJsonString(saveBCFModelToJsonString(b))
  expect_equal(
    predict(b2, d$X_test, d$Z_test, d$propensity_test,
      rfx_group_ids = d$groups_test, terms = "zero_rfx"),
    b$zero_rfx_preds_test
  )
})

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

test_that("zero_random_effects_params is validated", {
  skip_on_cran()
  d <- zero_prior_data(n = 80, n_test = 10)
  for (fit in list(prior_fit_bart, prior_fit_bcf)) {
    expect_error(fit(d, num_gfr = 0, num_mcmc = 2, zero_rfx = "tight"),
      "must be a list")
    expect_error(fit(d, num_gfr = 0, num_mcmc = 2, zero_rfx = list(1)),
      "named list")
    for (key in c("model_spec", "unseen_groups", "variance_prior_rate")) {
      bad <- setNames(list("intercept_only"), key)
      expect_error(fit(d, num_gfr = 0, num_mcmc = 2, zero_rfx = bad),
        "accepts only")
    }
    for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
      expect_error(
        fit(d, num_gfr = 0, num_mcmc = 2,
          zero_rfx = list(variance_prior_shape = bad)),
        "variance_prior_shape must be a single positive number"
      )
      expect_error(
        fit(d, num_gfr = 0, num_mcmc = 2,
          zero_rfx = list(variance_prior_scale = bad)),
        "variance_prior_scale must be a single positive number"
      )
    }
    # Only with family = "zinb", group ids and zero-part random effects on
    expect_error(
      fit(d, num_gfr = 0, num_mcmc = 2, zero_rfx = tight_prior,
        general = list(zero_random_effects = FALSE)),
      "applies only to family = 'zinb'"
    )
  }
  negbin <- OutcomeModel("count", "log")
  expect_error(
    bart(d$X, d$y, rfx_group_ids_train = d$groups,
      random_effects_params = list(model_spec = "intercept_only"),
      num_gfr = 0, num_mcmc = 2, zero_random_effects_params = tight_prior,
      general_params = list(outcome_model = negbin)),
    "applies only to family = 'zinb'"
  )
  expect_error(
    bcf(d$X, d$Z, d$y, propensity_train = d$propensity,
      rfx_group_ids_train = d$groups,
      random_effects_params = list(model_spec = "intercept_only"),
      num_gfr = 0, num_mcmc = 2, zero_random_effects_params = tight_prior,
      general_params = list(outcome_model = negbin)),
    "applies only to family = 'zinb'"
  )
  expect_error(
    bart(d$X, d$y, num_gfr = 0, num_mcmc = 2,
      general_params = list(outcome_model = zinb_prior_model()),
      zero_random_effects_params = tight_prior),
    "applies only to family = 'zinb'"
  )
  expect_error(
    bart(d$X, d$y, rfx_group_ids_train = d$groups,
      random_effects_params = list(model_spec = "intercept_only"),
      num_gfr = 0, num_mcmc = 2, zero_random_effects_params = tight_prior),
    "applies only to family = 'zinb'"
  )
  # An empty list or NULL is the default everywhere
  expect_no_error(
    bart(d$X, d$y, num_gfr = 0, num_mcmc = 2, zero_random_effects_params = NULL)
  )
})
