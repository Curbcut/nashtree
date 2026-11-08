# Generates zinb-rfx-reference.rds: seeded draws captured on branch
# cma-bcf-count (d8a0cb7b) before the zero part of the zero-inflated family got
# its own random intercept. test-zinb-rfx.R compares current draws with these
# bitwise to prove that the paths below did not change:
#   (a) ZINB bart() without group ids,
#   (b) ZINB bart() with group ids and general_params$zero_random_effects = FALSE,
#   (c) ZINB bcf() with group ids and general_params$zero_random_effects = FALSE,
#   (d) negative-binomial bart() and bcf() with group ids.
# d8a0cb7b ignores the unknown zero_random_effects entry of general_params, so
# the same calls run on both builds.
#
# Regenerate only from a build of d8a0cb7b:
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-zinb-rfx-reference.R
library(nashtree)

zinb_rfx_reference_fits <- function() {
  set.seed(7171)
  n <- 220
  n_test <- 40
  p <- 4
  X <- matrix(runif((n + n_test) * p), ncol = p)
  offset <- log(runif(n + n_test, 0.5, 2))
  groups <- sample(letters[1:6], n + n_test, replace = TRUE)
  u_count <- c(a = -0.4, b = 0, c = 0.3, d = 0.1, e = -0.2, f = 0.5)
  u_zero <- c(a = 0.5, b = -0.3, c = 0, d = 0.6, e = -0.6, f = 0.2)
  ps <- 0.3 + 0.4 * X[, 1]
  Z <- rbinom(n + n_test, 1, ps)
  log_mu <- offset + 1 + X[, 1] - X[, 2] + 0.3 * Z + u_count[groups]
  pi <- pnorm(-0.4 + X[, 3] - 0.3 * Z + u_zero[groups])
  y <- ifelse(runif(n + n_test) < pi, 0,
    rnbinom(n + n_test, mu = exp(log_mu), size = 4))
  tr <- seq_len(n)
  te <- n + seq_len(n_test)
  zinb <- OutcomeModel("count", "log", family = "zinb")
  negbin <- OutcomeModel("count", "log")

  keep_bart <- function(m, extra = list()) {
    c(list(
      y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
      zero_hat_train = m$zero_hat_train, zero_hat_test = m$zero_hat_test,
      dispersion = m$negbin_dispersion_samples,
      sigma2_leaf = m$sigma2_leaf_samples,
      rfx = m$rfx_preds_train, rfx_test = m$rfx_preds_test
    ), extra)
  }
  keep_bcf <- function(m, extra = list()) {
    c(list(
      y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
      mu_hat_train = m$mu_hat_train, tau_hat_train = m$tau_hat_train,
      tau_hat_test = m$tau_hat_test,
      zero_hat_train = m$zero_hat_train, zero_hat_test = m$zero_hat_test,
      zero_mu_hat_train = m$zero_mu_hat_train,
      zero_tau_hat_train = m$zero_tau_hat_train,
      dispersion = m$negbin_dispersion_samples,
      tau_0 = m$tau_0_samples, rfx = m$rfx_preds_train,
      rfx_test = m$rfx_preds_test
    ), extra)
  }
  out <- list()

  # (a) ZINB bart() without group ids: GFR warm start, two chains, test set
  a <- bart(X[tr, ], y[tr], X_test = X[te, ],
    offset_train = offset[tr], offset_test = offset[te],
    num_gfr = 3, num_burnin = 3, num_mcmc = 6,
    general_params = list(outcome_model = zinb, random_seed = 41L,
      num_chains = 2),
    mean_forest_params = list(num_trees = 30),
    zero_forest_params = list(num_trees = 15))
  set.seed(42)
  out$bart_zinb_no_rfx <- keep_bart(a, list(
    response = predict(a, X[te, ], offset = offset[te], terms = "y_hat",
      scale = "response"),
    zero_probability = predict(a, X[te, ], terms = "y_hat",
      scale = "zero_probability"),
    ppd = sampleBARTPosteriorPredictive(a, X[te, ], offset = offset[te],
      num_draws_per_sample = 2)))

  # (b) ZINB bart() with group ids, zero-part random effects switched off
  b <- bart(X[tr, ], y[tr], X_test = X[te, ],
    offset_train = offset[tr], offset_test = offset[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 3, num_burnin = 3, num_mcmc = 6,
    general_params = list(outcome_model = zinb, random_seed = 43L,
      num_chains = 2, zero_random_effects = FALSE),
    mean_forest_params = list(num_trees = 30),
    zero_forest_params = list(num_trees = 15))
  set.seed(44)
  out$bart_zinb_rfx_off <- keep_bart(b, list(
    response = predict(b, X[te, ], rfx_group_ids = groups[te],
      offset = offset[te], terms = "y_hat", scale = "response"),
    zero_probability = predict(b, X[te, ], rfx_group_ids = groups[te],
      terms = "y_hat", scale = "zero_probability"),
    ppd = sampleBARTPosteriorPredictive(b, X[te, ], rfx_group_ids = groups[te],
      offset = offset[te], num_draws_per_sample = 2)))

  # (c) ZINB bcf() with group ids, zero-part random effects switched off
  cz <- bcf(X[tr, ], Z[tr], y[tr], propensity_train = ps[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = ps[te],
    offset_train = offset[tr], offset_test = offset[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 3, num_burnin = 3, num_mcmc = 6,
    general_params = list(outcome_model = zinb, random_seed = 45L,
      num_chains = 2, zero_random_effects = FALSE),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15),
    zero_prognostic_forest_params = list(num_trees = 15),
    zero_treatment_effect_forest_params = list(num_trees = 10))
  set.seed(46)
  out$bcf_zinb_rfx_off <- keep_bcf(cz, list(
    response = predict(cz, X[te, ], Z[te], ps[te], rfx_group_ids = groups[te],
      offset = offset[te], terms = "y_hat", scale = "response"),
    zero_probability = predict(cz, X[te, ], Z[te], ps[te],
      rfx_group_ids = groups[te], terms = "y_hat", scale = "zero_probability"),
    contrast = computeContrastBCFModel(cz, X[te, ], X[te, ], rep(0, n_test),
      rep(1, n_test), ps[te], ps[te], groups[te], groups[te],
      matrix(1, n_test, 1), matrix(1, n_test, 1), scale = "response", offset_0 = offset[te], offset_1 = offset[te]),
    ppd = sampleBCFPosteriorPredictive(cz, X[te, ], Z[te], ps[te],
      rfx_group_ids = groups[te], offset = offset[te],
      num_draws_per_sample = 2)))

  # (d) negative-binomial bart() and bcf() with group ids
  d1 <- bart(X[tr, ], y[tr], X_test = X[te, ],
    offset_train = offset[tr], offset_test = offset[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 3, num_burnin = 3, num_mcmc = 6,
    general_params = list(outcome_model = negbin, random_seed = 47L,
      num_chains = 2),
    mean_forest_params = list(num_trees = 30))
  out$bart_negbin_rfx <- keep_bart(d1, list(
    rate = predict(d1, X[te, ], rfx_group_ids = groups[te],
      offset = offset[te], terms = "y_hat", scale = "rate")))
  d2 <- bcf(X[tr, ], Z[tr], y[tr], propensity_train = ps[tr],
    X_test = X[te, ], Z_test = Z[te], propensity_test = ps[te],
    offset_train = offset[tr], offset_test = offset[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 3, num_burnin = 3, num_mcmc = 6,
    general_params = list(outcome_model = negbin, random_seed = 48L,
      num_chains = 2),
    prognostic_forest_params = list(num_trees = 30),
    treatment_effect_forest_params = list(num_trees = 15))
  out$bcf_negbin_rfx <- keep_bcf(d2, list(
    rate = predict(d2, X[te, ], Z[te], ps[te], rfx_group_ids = groups[te],
      offset = offset[te], terms = "y_hat", scale = "rate")))
  out
}

if (sys.nframe() == 0L) {
  ref <- suppressWarnings(zinb_rfx_reference_fits())
  ref$platform <- R.version$platform
  ref$generated_from <- "cma-bcf-count d8a0cb7b"
  saveRDS(ref, file.path("tests", "testthat", "fixtures", "zinb-rfx-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
