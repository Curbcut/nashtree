# Generates bcf-reference.rds: seeded draws from the Gaussian and probit paths
# of bcf(), captured on branch cma-zinb (242aae46, whose R/bcf.R is that of
# cma-nb) before count outcomes were added to bcf(). test-bcf-count.R compares
# current draws with these bitwise to prove those paths did not change.
#
# Regenerate only from a build of 242aae46:
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-bcf-reference.R
library(nashtree)

bcf_reference_fits <- function() {
  set.seed(6161)
  n <- 160
  n_test <- 40
  p <- 4
  X <- matrix(runif((n + n_test) * p), ncol = p)
  pi_x <- 0.25 + 0.5 * X[, 1]
  Z <- rbinom(n + n_test, 1, pi_x)
  dose <- runif(n + n_test)
  mu_x <- 2 * sin(pi * X[, 2]) + X[, 3]
  tau_x <- 1 + X[, 4]
  groups <- sample(1:4, n + n_test, replace = TRUE)
  y_gauss <- mu_x + tau_x * Z + c(-0.5, 0, 0.3, 0.6)[groups] + rnorm(n + n_test, 0, 0.5)
  y_dose <- mu_x + tau_x * dose + rnorm(n + n_test, 0, 0.5)
  y_bin <- as.numeric(mu_x - 1 + tau_x * Z + rnorm(n + n_test) > 1)
  w <- runif(n, 0.5, 2)
  tr <- seq_len(n)
  te <- n + seq_len(n_test)

  keep <- function(m, extra = list()) {
    c(list(
      y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
      mu_hat_train = m$mu_hat_train, tau_hat_train = m$tau_hat_train,
      tau_hat_test = m$tau_hat_test, sigma2 = m$sigma2_global_samples,
      sigma2_leaf_mu = m$sigma2_leaf_mu_samples, tau_0 = m$tau_0_samples,
      b_0 = m$b_0_samples, b_1 = m$b_1_samples, rfx = m$rfx_preds_train
    ), extra)
  }
  out <- list()
  a <- bcf(X[tr, ], Z[tr], y_gauss[tr], propensity_train = pi_x[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = pi_x[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_plus_treatment"),
    num_gfr = 3, num_burnin = 3, num_mcmc = 8,
    general_params = list(random_seed = 31L, num_chains = 2),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15))
  out$gaussian_rfx_chains <- keep(a, list(
    pred = predict(a, X[te, ], Z[te], pi_x[te], rfx_group_ids = groups[te]),
    contrast = computeContrastBCFModel(a, X[te, ], X[te, ], rep(0, n_test),
      rep(1, n_test), pi_x[te], pi_x[te], groups[te], groups[te],
      cbind(1, rep(0, n_test)), cbind(1, rep(1, n_test)))))
  b <- bcf(X[tr, ], Z[tr], y_gauss[tr], propensity_train = pi_x[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = pi_x[te],
    num_gfr = 0, num_burnin = 4, num_mcmc = 8,
    general_params = list(random_seed = 32L, adaptive_coding = TRUE),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15))
  out$gaussian_adaptive_root <- keep(b)
  set.seed(33)
  out$gaussian_adaptive_root$ppd <- sampleBCFPosteriorPredictive(b, X[te, ],
    Z[te], pi_x[te], num_draws_per_sample = 2)
  cd <- bcf(X[tr, ], dose[tr], y_dose[tr], propensity_train = pi_x[tr],
    observation_weights = 1 / w, num_gfr = 2, num_mcmc = 6,
    general_params = list(random_seed = 34L),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15))
  out$gaussian_dose_weighted <- keep(cd)
  pr <- bcf(X[tr, ], Z[tr], y_bin[tr], propensity_train = pi_x[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = pi_x[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 2, num_burnin = 2, num_mcmc = 8,
    general_params = list(random_seed = 35L, num_chains = 2,
      outcome_model = OutcomeModel("binary", "probit")),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15))
  out$probit_rfx_chains <- keep(pr, list(
    prob = predict(pr, X[te, ], Z[te], pi_x[te], rfx_group_ids = groups[te],
      terms = "y_hat", scale = "probability")))
  out
}

if (sys.nframe() == 0L) {
  ref <- suppressWarnings(bcf_reference_fits())
  ref$platform <- R.version$platform
  ref$generated_from <- "cma-zinb 242aae46 (R/bcf.R as on cma-nb)"
  saveRDS(ref, file.path("tests", "testthat", "fixtures", "bcf-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
