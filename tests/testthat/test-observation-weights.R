make_bart_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  X <- matrix(runif(n * p), ncol = p)
  y <- sin(X[, 1] * pi) + rnorm(n, 0, 0.1)
  n_train <- as.integer(0.8 * n)
  list(
    X_train = X[1:n_train, ], y_train = y[1:n_train],
    X_test = X[(n_train + 1):n, ], n_train = n_train, n_test = n - n_train
  )
}

make_bcf_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  X <- matrix(runif(n * p), ncol = p)
  pi_X <- 0.25 + 0.5 * X[, 1]
  Z <- rbinom(n, 1, pi_X)
  y <- pi_X * 5 + X[, 2] * 2 * Z + rnorm(n, 0, 1)
  n_train <- as.integer(0.8 * n)
  list(
    X_train = X[1:n_train, ], Z_train = Z[1:n_train], y_train = y[1:n_train],
    pi_train = pi_X[1:n_train], X_test = X[(n_train + 1):n, ],
    Z_test = Z[(n_train + 1):n], pi_test = pi_X[(n_train + 1):n],
    n_train = n_train, n_test = n - n_train
  )
}

test_that("BART: uniform weights produce identical predictions to no weights", {
  skip_on_cran()
  d <- make_bart_data()
  num_mcmc <- 10

  set.seed(1)
  m1 <- bart(
    X_train = d$X_train, y_train = d$y_train, X_test = d$X_test,
    num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc,
    general_params = list(random_seed = 1L)
  )

  set.seed(1)
  m2 <- bart(
    X_train = d$X_train, y_train = d$y_train, X_test = d$X_test,
    observation_weights = rep(1.0, d$n_train),
    num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc,
    general_params = list(random_seed = 1L)
  )

  expect_equal(m1$y_hat_train, m2$y_hat_train)
  expect_equal(m1$y_hat_test, m2$y_hat_test)
})

test_that("BART: non-uniform weights run and produce correct output shape", {
  skip_on_cran()
  d <- make_bart_data()
  num_mcmc <- 10
  weights <- runif(d$n_train, 0.5, 2.0)

  expect_no_error(
    m <- bart(
      X_train = d$X_train, y_train = d$y_train, X_test = d$X_test,
      observation_weights = weights,
      num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc
    )
  )
  expect_equal(dim(m$y_hat_train), c(d$n_train, num_mcmc))
  expect_equal(dim(m$y_hat_test), c(d$n_test, num_mcmc))
})

test_that("BART: all-zero weights (prior mode) run with num_gfr = 0", {
  skip_on_cran()
  d <- make_bart_data()
  num_mcmc <- 10

  expect_no_error(
    m <- bart(
      X_train = d$X_train, y_train = d$y_train,
      observation_weights = rep(0.0, d$n_train),
      num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc
    )
  )
  expect_equal(dim(m$y_hat_train), c(d$n_train, num_mcmc))
})

test_that("BART: non-numeric observation_weights raises error", {
  skip_on_cran()
  d <- make_bart_data()
  expect_error(
    bart(
      X_train = d$X_train, y_train = d$y_train,
      observation_weights = as.character(rep(1, d$n_train)),
      num_gfr = 0, num_burnin = 0, num_mcmc = 5
    ),
    "numeric"
  )
})

test_that("BART: wrong-length observation_weights raises error", {
  skip_on_cran()
  d <- make_bart_data()
  expect_error(
    bart(
      X_train = d$X_train, y_train = d$y_train,
      observation_weights = rep(1.0, d$n_train + 1),
      num_gfr = 0, num_burnin = 0, num_mcmc = 5
    ),
    "nrow"
  )
})

test_that("BART: negative observation_weights raises error", {
  skip_on_cran()
  d <- make_bart_data()
  weights <- rep(1.0, d$n_train)
  weights[1] <- -1.0
  expect_error(
    bart(
      X_train = d$X_train, y_train = d$y_train,
      observation_weights = weights,
      num_gfr = 0, num_burnin = 0, num_mcmc = 5
    ),
    "negative"
  )
})

test_that("BART: all-zero weights with num_gfr > 0 raises error", {
  skip_on_cran()
  d <- make_bart_data()
  expect_error(
    bart(
      X_train = d$X_train, y_train = d$y_train,
      observation_weights = rep(0.0, d$n_train),
      num_gfr = 5, num_burnin = 0, num_mcmc = 10
    ),
    "num_gfr"
  )
})

test_that("BART: observation_weights with cloglog outcome raises error", {
  skip_on_cran()
  d <- make_bart_data()
  y_ord <- sample(1:3, d$n_train, replace = TRUE)
  expect_error(
    bart(
      X_train = d$X_train, y_train = y_ord,
      observation_weights = rep(1.0, d$n_train),
      num_gfr = 0, num_burnin = 0, num_mcmc = 5,
      general_params = list(outcome_model = OutcomeModel(outcome = "ordinal", link = "cloglog"))
    ),
    "cloglog"
  )
})

test_that("BART: observation_weights with variance forest raises warning", {
  skip_on_cran()
  d <- make_bart_data()
  expect_warning(
    bart(
      X_train = d$X_train, y_train = d$y_train,
      observation_weights = rep(1.0, d$n_train),
      num_gfr = 0, num_burnin = 0, num_mcmc = 5,
      variance_forest_params = list(num_trees = 5)
    ),
    "variance forest"
  )
})

test_that("BCF: uniform weights produce identical predictions to no weights", {
  skip_on_cran()
  d <- make_bcf_data()
  num_mcmc <- 10

  set.seed(1)
  m1 <- bcf(
    X_train = d$X_train, Z_train = d$Z_train, y_train = d$y_train,
    propensity_train = d$pi_train, X_test = d$X_test,
    Z_test = d$Z_test, propensity_test = d$pi_test,
    num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc,
    general_params = list(random_seed = 1L)
  )

  set.seed(1)
  m2 <- bcf(
    X_train = d$X_train, Z_train = d$Z_train, y_train = d$y_train,
    propensity_train = d$pi_train, X_test = d$X_test,
    Z_test = d$Z_test, propensity_test = d$pi_test,
    observation_weights = rep(1.0, d$n_train),
    num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc,
    general_params = list(random_seed = 1L)
  )

  expect_equal(m1$y_hat_train, m2$y_hat_train)
  expect_equal(m1$tau_hat_train, m2$tau_hat_train)
})

test_that("BCF: non-uniform weights run and produce correct output shape", {
  skip_on_cran()
  d <- make_bcf_data()
  num_mcmc <- 10
  weights <- runif(d$n_train, 0.5, 2.0)

  expect_no_error(
    m <- bcf(
      X_train = d$X_train, Z_train = d$Z_train, y_train = d$y_train,
      propensity_train = d$pi_train, X_test = d$X_test,
      Z_test = d$Z_test, propensity_test = d$pi_test,
      observation_weights = weights,
      num_gfr = 0, num_burnin = 0, num_mcmc = num_mcmc
    )
  )
  expect_equal(dim(m$y_hat_train), c(d$n_train, num_mcmc))
  expect_equal(dim(m$tau_hat_train), c(d$n_train, num_mcmc))
  expect_equal(dim(m$y_hat_test), c(d$n_test, num_mcmc))
  expect_equal(dim(m$tau_hat_test), c(d$n_test, num_mcmc))
})

test_that("BCF: negative observation_weights raises error", {
  skip_on_cran()
  d <- make_bcf_data()
  weights <- rep(1.0, d$n_train)
  weights[1] <- -1.0
  expect_error(
    bcf(
      X_train = d$X_train, Z_train = d$Z_train, y_train = d$y_train,
      propensity_train = d$pi_train,
      observation_weights = weights,
      num_gfr = 0, num_burnin = 0, num_mcmc = 5
    ),
    "negative"
  )
})

test_that("BCF: all-zero weights with num_gfr > 0 raises error", {
  skip_on_cran()
  d <- make_bcf_data()
  expect_error(
    bcf(
      X_train = d$X_train, Z_train = d$Z_train, y_train = d$y_train,
      propensity_train = d$pi_train,
      observation_weights = rep(0.0, d$n_train),
      num_gfr = 5, num_burnin = 0, num_mcmc = 10
    ),
    "num_gfr"
  )
})

test_that("BCF: observation_weights with cloglog outcome raises error", {
  skip_on_cran()
  d <- make_bcf_data()
  y_bin <- rbinom(d$n_train, 1, 0.5)
  expect_error(
    bcf(
      X_train = d$X_train, Z_train = d$Z_train, y_train = y_bin,
      propensity_train = d$pi_train,
      observation_weights = rep(1.0, d$n_train),
      num_gfr = 0, num_burnin = 0, num_mcmc = 5,
      general_params = list(outcome_model = OutcomeModel(outcome = "binary", link = "cloglog"))
    ),
    "cloglog"
  )
})

# ---------------------------------------------------------------------------
# observation_weights are precision weights: y_i ~ N(mu(X_i), sigma^2 / w_i)
# ---------------------------------------------------------------------------

# X carries no signal; half the rows sit at y = 1, half at y = 0. Weighting one
# half pulls the fit toward it. Until 0.4.6.9013 the weights reached the
# sampler as variance multipliers and pulled the other way.
make_two_level_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  list(
    X = matrix(runif(n * 2), ncol = 2),
    X_test = matrix(runif(100 * 2), ncol = 2),
    high = rep(c(TRUE, FALSE), each = n / 2),
    y = rep(c(1, 0), each = n / 2) + rnorm(n, sd = 0.05)
  )
}

fit_two_level <- function(d, w, ...) {
  bart(
    X_train = d$X, y_train = d$y, X_test = d$X_test, observation_weights = w,
    num_gfr = 10, num_burnin = 0, num_mcmc = 50,
    general_params = list(random_seed = 7L, num_threads = 1L),
    mean_forest_params = list(num_trees = 20, min_samples_leaf = 50), ...
  )
}

test_that("BART: a larger weight gives a row more influence", {
  skip_on_cran()
  d <- make_two_level_data()
  unweighted <- mean(fit_two_level(d, NULL)$y_hat_test)
  toward_high <- mean(fit_two_level(d, ifelse(d$high, 100, 1))$y_hat_test)
  toward_low <- mean(fit_two_level(d, ifelse(d$high, 0.01, 1))$y_hat_test)
  expect_gt(toward_high, 0.9)
  expect_lt(toward_low, 0.1)
  expect_gt(toward_high, unweighted)
  expect_lt(toward_low, unweighted)
})

test_that("BART: a constant weight c multiplies the global variance by c", {
  skip_on_cran()
  d <- make_two_level_data()
  base <- mean(fit_two_level(d, NULL)$sigma2_global_samples)
  scaled <- mean(fit_two_level(d, rep(4, length(d$y)))$sigma2_global_samples)
  expect_gt(scaled / base, 3.2)
  expect_lt(scaled / base, 5)
})

test_that("BART: random effects see the observation weights", {
  skip_on_cran()
  # Two groups hold the same labels with opposite weights: group 1 weights its
  # y = 1 rows, group 2 its y = 0 rows. The forest cannot see the groups, so
  # only weighted random effects separate the two intercepts (before
  # 0.4.6.9013 a Gaussian fit gave its random effects no weights at all).
  d <- make_two_level_data()
  groups <- rep(1:2, length(d$y) / 2)
  w <- ifelse((groups == 1) == d$high, 100, 1)
  m <- bart(
    X_train = d$X, y_train = d$y, observation_weights = w,
    rfx_group_ids_train = groups,
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 0, num_burnin = 50, num_mcmc = 50,
    general_params = list(random_seed = 8L, num_threads = 1L),
    mean_forest_params = list(num_trees = 20, min_samples_leaf = 1000)
  )
  rfx <- rowMeans(m$rfx_preds_train)
  expect_gt(mean(rfx[groups == 1]) - mean(rfx[groups == 2]), 0.8)
})

test_that("BART: zero-weight rows count as absent from the likelihood", {
  skip_on_cran()
  # 300 of 400 rows sit at y = 1 with weight 0; the other 100 at y = 0. The
  # weighted fit must match a fit on the 100 rows alone, in its mean and in
  # its global variance (zero-weight rows count towards neither the residual
  # sum of squares nor the shape of sigma^2's posterior). The outcome is still
  # standardized over all rows, so the match is close, not bitwise.
  d <- make_two_level_data()
  absent <- seq_len(300)
  d$y <- rnorm(length(d$y), sd = 0.05)
  d$y[absent] <- 1
  w <- rep(1, length(d$y))
  w[absent] <- 0
  weighted <- fit_two_level(d, w)
  kept <- d
  kept$X <- d$X[-absent, , drop = FALSE]
  kept$y <- d$y[-absent]
  subset <- fit_two_level(kept, NULL)
  expect_lt(abs(mean(weighted$y_hat_test)), 0.02)
  expect_lt(
    abs(log(mean(weighted$sigma2_global_samples) /
      mean(subset$sigma2_global_samples))),
    log(1.5)
  )
})

test_that("BCF: a larger weight gives a row more influence", {
  skip_on_cran()
  d <- make_two_level_data()
  set.seed(3)
  Z <- rbinom(length(d$y), 1, 0.5)
  fit <- function(w) {
    m <- bcf(
      X_train = d$X, Z_train = Z, y_train = d$y,
      propensity_train = rep(0.5, length(d$y)), observation_weights = w,
      num_gfr = 10, num_burnin = 0, num_mcmc = 30,
      general_params = list(random_seed = 9L, num_threads = 1L),
      prognostic_forest_params = list(num_trees = 20, min_samples_leaf = 50),
      treatment_effect_forest_params = list(num_trees = 10, min_samples_leaf = 50)
    )
    mean(m$y_hat_train)
  }
  expect_gt(fit(ifelse(d$high, 100, 1)), 0.85)
  expect_lt(fit(ifelse(d$high, 0.01, 1)), 0.15)
})
