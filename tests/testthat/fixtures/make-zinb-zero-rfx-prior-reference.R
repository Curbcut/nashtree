# Generates zinb-zero-rfx-prior-reference.rds: seeded draws captured on branch
# cma-rfx (1273f128) before the zero part of the zero-inflated family got a
# random-effects prior of its own (zero_random_effects_params).
# test-zinb-zero-rfx-prior.R compares current draws with these bitwise to show
# that the paths below did not change:
#   (a) Gaussian bart() and bcf() with group ids,
#   (b) probit bart() and bcf() with group ids, cloglog bart(),
#   (c) negative-binomial bart() and bcf() with group ids,
#   (d) ZINB bart() and bcf() with zero-part random effects and the default
#       (empty) zero_random_effects_params: GFR-started chains, root-started
#       chains and warm starts, so every place the zero-part random effects
#       are initialised or reset is exercised.
# Every call leaves zero_random_effects_params at its default, so the same
# calls run on both builds.
#
# Regenerate only from a build of 1273f128:
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-zinb-zero-rfx-prior-reference.R
library(nashtree)

zinb_zero_rfx_prior_reference_fits <- function() {
  set.seed(9191)
  n <- 200
  n_test <- 30
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
  y_gauss <- 1 + X[, 1] - X[, 2] + 0.5 * Z + u_count[groups] +
    rnorm(n + n_test, 0, 0.5)
  y_bin <- as.numeric(y_gauss > median(y_gauss))
  tr <- seq_len(n)
  te <- n + seq_len(n_test)
  zinb <- OutcomeModel("count", "log", family = "zinb")
  negbin <- OutcomeModel("count", "log")
  probit <- OutcomeModel("binary", "probit")
  cloglog <- OutcomeModel("binary", "cloglog")
  rfx_spec <- list(model_spec = "intercept_only")

  keep_bart <- function(m) {
    list(
      y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
      zero_hat_train = m$zero_hat_train, zero_hat_test = m$zero_hat_test,
      sigma2 = m$sigma2_global_samples, sigma2_leaf = m$sigma2_leaf_samples,
      dispersion = m$negbin_dispersion_samples,
      rfx = m$rfx_preds_train, rfx_test = m$rfx_preds_test,
      zero_rfx = m$zero_rfx_preds_train, zero_rfx_test = m$zero_rfx_preds_test
    )
  }
  keep_bcf <- function(m) {
    list(
      y_hat_train = m$y_hat_train, y_hat_test = m$y_hat_test,
      mu_hat_train = m$mu_hat_train, tau_hat_train = m$tau_hat_train,
      tau_hat_test = m$tau_hat_test, sigma2 = m$sigma2_global_samples,
      zero_hat_train = m$zero_hat_train, zero_hat_test = m$zero_hat_test,
      zero_tau_hat_train = m$zero_tau_hat_train,
      dispersion = m$negbin_dispersion_samples, tau_0 = m$tau_0_samples,
      rfx = m$rfx_preds_train, rfx_test = m$rfx_preds_test,
      zero_rfx = m$zero_rfx_preds_train, zero_rfx_test = m$zero_rfx_preds_test
    )
  }
  fit_bart <- function(yy, model, seed, num_gfr = 3, num_chains = 2,
                       offsets = FALSE, extra = list(), ...) {
    args <- list(X_train = X[tr, ], y_train = yy[tr], X_test = X[te, ],
      rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
      random_effects_params = rfx_spec,
      num_gfr = num_gfr, num_burnin = 3, num_mcmc = 6,
      general_params = c(list(outcome_model = model, random_seed = seed,
        num_chains = num_chains), extra),
      mean_forest_params = list(num_trees = 30), ...)
    if (offsets) {
      args$offset_train <- offset[tr]
      args$offset_test <- offset[te]
    }
    do.call(bart, args)
  }
  fit_bcf <- function(yy, model, seed, num_gfr = 3, num_chains = 2,
                      offsets = FALSE, extra = list(), ...) {
    args <- list(X_train = X[tr, ], Z_train = Z[tr], y_train = yy[tr],
      propensity_train = ps[tr],
      X_test = X[te, ], Z_test = Z[te], propensity_test = ps[te],
      rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
      random_effects_params = rfx_spec,
      num_gfr = num_gfr, num_burnin = 3, num_mcmc = 6,
      general_params = c(list(outcome_model = model, random_seed = seed,
        num_chains = num_chains), extra),
      prognostic_forest_params = list(num_trees = 30),
      treatment_effect_forest_params = list(num_trees = 15), ...)
    if (offsets) {
      args$offset_train <- offset[tr]
      args$offset_test <- offset[te]
    }
    do.call(bcf, args)
  }
  out <- list()

  # (a) Gaussian
  out$bart_gaussian_rfx <- keep_bart(fit_bart(y_gauss,
    OutcomeModel("continuous", "identity"), 51L))
  out$bcf_gaussian_rfx <- keep_bcf(fit_bcf(y_gauss,
    OutcomeModel("continuous", "identity"), 52L))

  # (b) probit and cloglog
  out$bart_probit_rfx <- keep_bart(fit_bart(y_bin, probit, 53L,
    extra = list(sample_sigma2_global = FALSE)))
  out$bcf_probit_rfx <- keep_bcf(fit_bcf(y_bin, probit, 54L,
    extra = list(sample_sigma2_global = FALSE)))
  cl <- bart(X[tr, ], y_bin[tr], X_test = X[te, ],
    num_gfr = 3, num_burnin = 3, num_mcmc = 6,
    general_params = list(outcome_model = cloglog, random_seed = 55L,
      sample_sigma2_global = FALSE, num_chains = 2),
    mean_forest_params = list(num_trees = 30))
  out$bart_cloglog <- keep_bart(cl)

  # (c) negative binomial
  out$bart_negbin_rfx <- keep_bart(fit_bart(y, negbin, 56L, offsets = TRUE))
  out$bcf_negbin_rfx <- keep_bcf(fit_bcf(y, negbin, 57L, offsets = TRUE))

  # (d) ZINB with zero-part random effects and the default prior
  zb_gfr <- fit_bart(y, zinb, 58L, offsets = TRUE,
    zero_forest_params = list(num_trees = 15))
  zb_root <- fit_bart(y, zinb, 59L, num_gfr = 0, offsets = TRUE,
    zero_forest_params = list(num_trees = 15))
  zb_warm <- fit_bart(y, zinb, 60L, num_gfr = 0, num_chains = 1,
    offsets = TRUE, zero_forest_params = list(num_trees = 15),
    previous_model_json = saveBARTModelToJsonString(zb_gfr),
    previous_model_warmstart_sample_num = 12)
  out$bart_zinb_zero_rfx_gfr <- keep_bart(zb_gfr)
  out$bart_zinb_zero_rfx_root <- keep_bart(zb_root)
  out$bart_zinb_zero_rfx_warm <- keep_bart(zb_warm)
  zc_args <- list(zero_prognostic_forest_params = list(num_trees = 15),
    zero_treatment_effect_forest_params = list(num_trees = 10))
  zc_gfr <- do.call(fit_bcf, c(list(y, zinb, 61L, offsets = TRUE), zc_args))
  zc_root <- do.call(fit_bcf, c(list(y, zinb, 62L, num_gfr = 0,
    offsets = TRUE), zc_args))
  zc_warm <- do.call(fit_bcf, c(list(y, zinb, 63L, num_gfr = 0,
    num_chains = 1, offsets = TRUE,
    previous_model_json = saveBCFModelToJsonString(zc_gfr),
    previous_model_warmstart_sample_num = 12), zc_args))
  out$bcf_zinb_zero_rfx_gfr <- keep_bcf(zc_gfr)
  out$bcf_zinb_zero_rfx_root <- keep_bcf(zc_root)
  out$bcf_zinb_zero_rfx_warm <- keep_bcf(zc_warm)
  out
}

if (sys.nframe() == 0L) {
  ref <- suppressWarnings(zinb_zero_rfx_prior_reference_fits())
  ref$platform <- R.version$platform
  ref$generated_from <- "cma-rfx 1273f128"
  saveRDS(ref, file.path("tests", "testthat", "fixtures",
    "zinb-zero-rfx-prior-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
