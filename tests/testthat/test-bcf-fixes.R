# Fixes to bcf() checked against seeded draws captured before them
# (fixtures/make-bcf-fixes-reference.R, cma-rfx 84f5ccbb).

# Fits whose settings no fix touches: must be bitwise identical to the reference.
bcf_fixes_unchanged <- c(
  "gaussian_intercept", "gaussian_adaptive", "gaussian_intercept_adaptive",
  "gaussian_intercept_adaptive_rfx_chains", "gaussian_variance_forest",
  "gaussian_tau_leaf_root", "probit_weighted_intercept_adaptive"
)
# Fits in the settings a fix targets: must differ from the reference.
bcf_fixes_changed <- c(
  "gaussian_tau_leaf_gfr", "gaussian_weighted_intercept",
  "gaussian_weighted_adaptive", "gaussian_weighted_intercept_adaptive"
)

test_that("bcf() draws match the pre-fix reference except where a fix applies", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "bcf-fixes-reference.rds"))
  skip_if_not(
    identical(ref$platform, R.version$platform),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-bcf-fixes-reference.R"), local = TRUE)
  current <- suppressWarnings(bcf_fixes_reference_fits())
  expect_setequal(
    c(bcf_fixes_unchanged, bcf_fixes_changed),
    setdiff(names(ref), c("platform", "generated_from"))
  )
  for (name in bcf_fixes_unchanged) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
  for (name in bcf_fixes_changed) {
    expect_false(identical(current[[name]], ref[[name]]), label = name)
  }
})

# ---------------------------------------------------------------------------
# GFR updates the treatment forest's leaf scale
# ---------------------------------------------------------------------------

# Runs `expr` with every ForestModelConfig created by bcf() instrumented to
# record each update_leaf_model_scale() call: the config's leaf model type,
# the value passed, and the scale the config holds afterwards.
record_leaf_scale_updates <- function(expr) {
  rec <- new.env()
  rec$updates <- list()
  ns <- asNamespace("nashtree")
  suppressMessages(trace(
    "createForestModelConfig",
    where = ns,
    print = FALSE,
    exit = bquote({
      cfg <- returnValue()
      rec <- .(rec)
      update <- cfg$update_leaf_model_scale
      unlockBinding("update_leaf_model_scale", cfg)
      cfg$update_leaf_model_scale <- function(leaf_model_scale) {
        update(leaf_model_scale)
        rec$updates[[length(rec$updates) + 1L]] <- list(
          type = cfg$leaf_model_type,
          value = as.numeric(leaf_model_scale),
          held = as.numeric(cfg$get_leaf_model_scale())
        )
      }
      lockBinding("update_leaf_model_scale", cfg)
    })
  ))
  on.exit(suppressMessages(untrace("createForestModelConfig", where = ns)))
  result <- expr
  list(result = result, updates = rec$updates)
}

test_that("GFR iterations update the treatment forest's leaf scale", {
  skip_on_cran()
  set.seed(8181)
  n <- 150
  X <- matrix(runif(n * 4), ncol = 4)
  pi_x <- 0.25 + 0.5 * X[, 1]
  Z <- rbinom(n, 1, pi_x)
  y <- 2 * sin(pi * X[, 2]) + (1 + X[, 4]) * Z + rnorm(n, 0, 0.5)
  num_gfr <- 4
  rec <- record_leaf_scale_updates(bcf(X, Z, y,
    propensity_train = pi_x, num_gfr = num_gfr, num_burnin = 0, num_mcmc = 3,
    general_params = list(random_seed = 9L, keep_gfr = TRUE),
    prognostic_forest_params = list(num_trees = 20),
    treatment_effect_forest_params = list(num_trees = 10,
      sample_sigma2_leaf = TRUE)))
  m <- rec$result
  by_type <- function(type) {
    u <- Filter(function(x) x$type == type, rec$updates)
    list(
      value = vapply(u, `[[`, numeric(1), "value"),
      held = vapply(u, `[[`, numeric(1), "held")
    )
  }
  # Expected sequence of updates for each forest: one per GFR draw, the reset
  # to the last GFR draw at the start of the chain, then one per MCMC draw.
  expected <- function(samples) {
    c(samples[seq_len(num_gfr)], samples[num_gfr], samples[-seq_len(num_gfr)])
  }
  tau <- by_type(1L)
  expect_length(m$sigma2_leaf_tau_samples, num_gfr + 3)
  expect_identical(tau$value, expected(m$sigma2_leaf_tau_samples))
  expect_identical(tau$held, tau$value)
  mu <- by_type(0L)
  expect_identical(mu$value, expected(m$sigma2_leaf_mu_samples))
  expect_identical(mu$held, mu$value)
})

# ---------------------------------------------------------------------------
# Gaussian tau_0 and adaptive coding use observation_weights
# ---------------------------------------------------------------------------

test_that("bcf_tau_0_posterior() is the weighted least-squares posterior", {
  set.seed(9191)
  n <- 50
  r <- rnorm(n)
  v <- runif(n, 0.2, 5)
  sigma2 <- 0.7
  prior_var <- 3
  # One basis column: closed form
  z <- rbinom(n, 1, 0.5)
  post <- nashtree:::bcf_tau_0_posterior(as.matrix(z), r, 1 / v, sigma2,
    prior_var)
  prec <- sum(z^2 / v) / sigma2 + 1 / prior_var
  expect_equal(post$mean, sum(z * r / v) / sigma2 / prec)
  expect_equal(as.numeric(post$cov), 1 / prec)
  # Two basis columns: Bayesian WLS as OLS on rows scaled by sqrt(1 / v),
  # with the prior as p pseudo-rows
  B <- cbind(z, runif(n))
  post <- nashtree:::bcf_tau_0_posterior(B, r, 1 / v, sigma2, prior_var)
  Xa <- rbind(B / sqrt(v * sigma2), diag(2) / sqrt(prior_var))
  ya <- c(r / sqrt(v * sigma2), 0, 0)
  expect_equal(post$mean, as.numeric(qr.solve(Xa, ya)))
  expect_equal(post$cov, solve(crossprod(Xa)), ignore_attr = TRUE)
  # Unit weights reproduce the unweighted path bitwise
  expect_identical(
    nashtree:::bcf_tau_0_posterior(B, r, rep(1, n), sigma2, prior_var),
    nashtree:::bcf_tau_0_posterior(B, r, NULL, sigma2, prior_var)
  )
})

test_that("bcf_coding_suff_stats() weights each row by its precision", {
  set.seed(9292)
  n <- 40
  tau_x <- rnorm(n, 1)
  r <- rnorm(n)
  Z <- rbinom(n, 1, 0.5)
  v <- runif(n, 0.2, 5)
  s <- nashtree:::bcf_coding_suff_stats(tau_x, r, Z, 1 / v)
  c0 <- Z == 0
  expect_equal(s$s_tt0, sum(tau_x[c0]^2 / v[c0]))
  expect_equal(s$s_tt1, sum(tau_x[!c0]^2 / v[!c0]))
  expect_equal(s$s_ty0, sum(tau_x[c0] * r[c0] / v[c0]))
  expect_equal(s$s_ty1, sum(tau_x[!c0] * r[!c0] / v[!c0]))
  expect_identical(
    nashtree:::bcf_coding_suff_stats(tau_x, r, Z, rep(1, n)),
    nashtree:::bcf_coding_suff_stats(tau_x, r, Z, NULL)
  )
})

# Runs `expr` with the tau_0 and adaptive-coding helpers instrumented to record
# their arguments and results on every call.
record_tau_0_and_coding <- function(expr) {
  rec <- new.env()
  rec$tau_0 <- list()
  rec$coding <- list()
  ns <- asNamespace("nashtree")
  record <- function(what) {
    # as.list() leaves out the dot-named binding
    bquote({
      .rec <- .(rec)
      .rec[[.(what)]][[length(.rec[[.(what)]]) + 1L]] <- c(
        as.list(environment()),
        list(result = returnValue())
      )
    })
  }
  suppressMessages({
    trace("bcf_tau_0_posterior", where = ns, print = FALSE,
      exit = record("tau_0"))
    trace("bcf_coding_suff_stats", where = ns, print = FALSE,
      exit = record("coding"))
  })
  on.exit(suppressMessages({
    untrace("bcf_tau_0_posterior", where = ns)
    untrace("bcf_coding_suff_stats", where = ns)
  }))
  result <- expr
  list(result = result, tau_0 = rec$tau_0, coding = rec$coding)
}

test_that("Gaussian tau_0 and coding updates condition on the weighted likelihood", {
  skip_on_cran()
  set.seed(9393)
  n <- 200
  X <- matrix(runif(n * 3), ncol = 3)
  pi_x <- 0.3 + 0.4 * X[, 1]
  Z <- rbinom(n, 1, pi_x)
  v <- ifelse(runif(n) < 0.2, 25, 1)
  y <- 1 + 2 * X[, 2] + Z + rnorm(n, 0, 0.5 * sqrt(v))
  fit <- function(...) {
    bcf(X, Z, y, propensity_train = pi_x, num_gfr = 2, num_burnin = 0,
      num_mcmc = 3, prognostic_forest_params = list(num_trees = 20),
      treatment_effect_forest_params = list(num_trees = 10), ...)
  }
  rec <- record_tau_0_and_coding(fit(observation_weights = v,
    general_params = list(random_seed = 5L, adaptive_coding = TRUE)))
  # One tau_0 and one coding update per GFR and MCMC iteration
  expect_length(rec$tau_0, 5)
  expect_length(rec$coding, 5)
  for (call in rec$tau_0) {
    expect_identical(call$precision, 1 / v)
    # Exact conditional mean of tau_0 given the forests: weighted least
    # squares of the partial residual on the basis, shrunk by the prior
    b <- as.numeric(call$basis)
    r <- as.numeric(call$partial_resid)
    wls <- sum(b * r / v) / (sum(b^2 / v) + call$sigma2 / call$prior_var)
    expect_equal(call$result$mean, wls)
    expect_equal(as.numeric(call$result$cov),
      call$sigma2 / (sum(b^2 / v) + call$sigma2 / call$prior_var))
  }
  for (call in rec$coding) {
    expect_identical(call$precision, 1 / v)
    t1 <- as.numeric(call$tau_x)[Z == 1]
    expect_equal(call$result$s_tt1, sum(t1^2 / v[Z == 1]))
  }

  # No weights, and probit with weights, keep the unweighted sums
  rec <- record_tau_0_and_coding(fit(
    general_params = list(random_seed = 5L, adaptive_coding = TRUE)))
  expect_true(all(vapply(rec$tau_0, function(x) is.null(x$precision), NA)))
  expect_true(all(vapply(rec$coding, function(x) is.null(x$precision), NA)))
  rec <- record_tau_0_and_coding(suppressWarnings(bcf(X, Z, as.numeric(y > 2),
    propensity_train = pi_x, observation_weights = v, num_gfr = 1,
    num_burnin = 0, num_mcmc = 2,
    general_params = list(random_seed = 5L, adaptive_coding = TRUE,
      outcome_model = OutcomeModel("binary", "probit")),
    prognostic_forest_params = list(num_trees = 20),
    treatment_effect_forest_params = list(num_trees = 10))))
  expect_length(rec$tau_0, 3)
  expect_true(all(vapply(rec$tau_0, function(x) is.null(x$precision), NA)))
})

test_that("Weighted Gaussian BCF recovers tau_0 under extreme heteroskedasticity", {
  skip_on_cran()
  # A fifth of the rows have 400 times the noise variance, passed as
  # observation_weights; the true effect is the constant 1. With the
  # unweighted tau_0 update the posterior mean of tau_0 was 38.7 here (the
  # weighted forest offsetting it); weighted, it stays near 1.
  set.seed(9494)
  n <- 500
  X <- matrix(runif(n * 3), ncol = 3)
  pi_x <- 0.3 + 0.4 * X[, 1]
  Z <- rbinom(n, 1, pi_x)
  v <- ifelse(runif(n) < 0.2, 400, 1)
  y <- 1 + 2 * X[, 2] + Z + rnorm(n, 0, 0.5 * sqrt(v))
  m <- bcf(X, Z, y, propensity_train = pi_x, observation_weights = v,
    num_gfr = 5, num_burnin = 50, num_mcmc = 100,
    general_params = list(random_seed = 3L))
  expect_lt(abs(mean(m$tau_0_samples) - 1), 1)
  expect_lt(abs(mean(m$tau_hat_train) - 1), 0.25)
})
