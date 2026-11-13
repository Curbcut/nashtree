# Generates bcf-fixes-reference.rds: seeded draws from bcf() captured on
# cma-rfx (84f5ccbb) before two fixes in R/bcf.R: the GFR loop updating the
# treatment forest's leaf scale, and observation weights in the Gaussian tau_0
# and adaptive-coding updates. test-bcf-fixes.R compares current draws with
# these: fits in settings the fixes do not touch must be bitwise identical,
# fits in the settings they target must differ.
#
# Regenerate only from a build of 84f5ccbb:
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-bcf-fixes-reference.R
library(nashtree)

bcf_fixes_reference_fits <- function() {
  set.seed(7171)
  n <- 160
  n_test <- 40
  p <- 4
  X <- matrix(runif((n + n_test) * p), ncol = p)
  pi_x <- 0.25 + 0.5 * X[, 1]
  Z <- rbinom(n + n_test, 1, pi_x)
  mu_x <- 2 * sin(pi * X[, 2]) + X[, 3]
  tau_x <- 1 + X[, 4]
  groups <- sample(1:4, n + n_test, replace = TRUE)
  y <- mu_x + tau_x * Z + c(-0.5, 0, 0.3, 0.6)[groups] +
    rnorm(n + n_test, 0, 0.5)
  y_bin <- as.numeric(mu_x - 1 + tau_x * Z + rnorm(n + n_test) > 1)
  w <- runif(n, 0.5, 2)
  tr <- seq_len(n)
  te <- n + seq_len(n_test)

  keep <- function(m) {
    list(
      y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
      mu_hat_train = m$mu_hat_train, tau_hat_train = m$tau_hat_train,
      tau_hat_test = m$tau_hat_test, sigma2 = m$sigma2_global_samples,
      sigma2_leaf_mu = m$sigma2_leaf_mu_samples,
      sigma2_leaf_tau = m$sigma2_leaf_tau_samples, tau_0 = m$tau_0_samples,
      b_0 = m$b_0_samples, b_1 = m$b_1_samples, rfx = m$rfx_preds_train,
      sigma2_x = m$sigma2_x_hat_train
    )
  }
  fit <- function(seed, y_train = y[tr], tau_params = list(), general = list(),
                  ...) {
    m <- bcf(X[tr, ], Z[tr], y_train, propensity_train = pi_x[tr],
      X_test = X[te, ], Z_test = Z[te], propensity_test = pi_x[te],
      num_gfr = 3, num_burnin = 3, num_mcmc = 6,
      general_params = c(list(random_seed = seed), general),
      prognostic_forest_params = list(num_trees = 30),
      treatment_effect_forest_params = c(list(num_trees = 15), tau_params),
      ...)
    keep(m)
  }
  out <- list()
  out$gaussian_intercept <- fit(41L)
  out$gaussian_adaptive <- fit(42L,
    tau_params = list(sample_intercept = FALSE),
    general = list(adaptive_coding = TRUE))
  out$gaussian_intercept_adaptive <- fit(43L,
    general = list(adaptive_coding = TRUE))
  m <- bcf(X[tr, ], Z[tr], y[tr], propensity_train = pi_x[tr],
    rfx_group_ids_train = groups[tr],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 2, num_burnin = 2, num_mcmc = 5,
    general_params = list(random_seed = 44L, num_chains = 2,
      adaptive_coding = TRUE),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15))
  out$gaussian_intercept_adaptive_rfx_chains <- keep(m)
  out$gaussian_variance_forest <- fit(45L,
    general = list(adaptive_coding = TRUE),
    variance_forest_params = list(num_trees = 10))
  out$gaussian_tau_leaf_root <- keep(bcf(X[tr, ], Z[tr], y[tr],
    propensity_train = pi_x[tr], num_gfr = 0, num_burnin = 3, num_mcmc = 6,
    general_params = list(random_seed = 46L),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15,
      sample_sigma2_leaf = TRUE)))
  # Weights are passed as 1 / w: the reference was captured when nashtree read
  # observation_weights as variance multipliers (0.4.6.9013 made them precision
  # weights), so 1 / w reproduces its draws.
  out$probit_weighted_intercept_adaptive <- fit(47L, y_train = y_bin[tr],
    general = list(adaptive_coding = TRUE,
      outcome_model = OutcomeModel("binary", "probit")),
    observation_weights = 1 / w)
  out$gaussian_tau_leaf_gfr <- fit(48L,
    tau_params = list(sample_sigma2_leaf = TRUE),
    general = list(keep_gfr = TRUE))
  out$gaussian_weighted_intercept <- fit(49L, observation_weights = 1 / w)
  out$gaussian_weighted_adaptive <- fit(50L,
    tau_params = list(sample_intercept = FALSE),
    general = list(adaptive_coding = TRUE), observation_weights = 1 / w)
  out$gaussian_weighted_intercept_adaptive <- fit(51L,
    general = list(adaptive_coding = TRUE), observation_weights = 1 / w)
  out
}

if (sys.nframe() == 0L) {
  ref <- suppressWarnings(bcf_fixes_reference_fits())
  ref$platform <- R.version$platform
  ref$generated_from <- "cma-rfx 84f5ccbb"
  saveRDS(ref, file.path("tests", "testthat", "fixtures",
    "bcf-fixes-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
