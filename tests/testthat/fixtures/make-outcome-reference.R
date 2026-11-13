# Generates outcome-paths-reference.rds: seeded draws from the Gaussian,
# probit and cloglog paths of bart(), captured before the negative-binomial
# outcome was added (bart.R as of 0ac0e7d0). The same draws are reproduced
# bitwise by bd5c5c48 (the nashtree rename), the base of the negative-binomial
# branch. test-negbin.R compares current draws with these bitwise to prove
# those paths did not change.
#
# Regenerate only from a build whose R/bart.R predates the count outcome:
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-outcome-reference.R
library(nashtree)

outcome_reference_fits <- function() {
  set.seed(4242)
  n <- 150
  p <- 4
  X <- matrix(runif(n * p), ncol = p)
  X_test <- matrix(runif(40 * p), ncol = p)
  f_X <- 2 * sin(pi * X[, 1]) + X[, 2]
  y_gauss <- f_X + rnorm(n, 0, 0.5)
  y_bin <- as.numeric(f_X + rnorm(n) > 1.5)
  y_ord <- cut(f_X + rnorm(n), c(-Inf, 0.8, 1.8, Inf), labels = FALSE)
  w <- runif(n, 0.5, 2)
  groups <- sample(1:4, n, replace = TRUE)
  groups_test <- sample(1:4, 40, replace = TRUE)

  keep <- function(m, extra = list()) {
    c(list(
      y_hat_train = m$y_hat_train,
      y_hat_test = m$y_hat_test,
      sigma2 = m$sigma2_global_samples,
      sigma2_leaf = m$sigma2_leaf_samples
    ), extra)
  }
  out <- list()
  g <- bart(X, y_gauss, X_test = X_test, num_gfr = 3, num_burnin = 3,
    num_mcmc = 8, general_params = list(random_seed = 1L, num_chains = 2))
  out$gaussian <- keep(g, list(pred = predict(g, X_test, type = "mean")))
  # Weights are passed as 1 / w: the reference was captured when nashtree read
  # observation_weights as variance multipliers (0.4.6.9013 made them precision
  # weights), so 1 / w reproduces its draws.
  gw <- bart(X, y_gauss, X_test = X_test, observation_weights = 1 / w,
    num_gfr = 0, num_burnin = 3, num_mcmc = 8,
    general_params = list(random_seed = 2L))
  out$gaussian_weighted <- keep(gw)
  gv <- bart(X, y_gauss, X_test = X_test, num_gfr = 2, num_mcmc = 6,
    general_params = list(random_seed = 3L),
    variance_forest_params = list(num_trees = 20))
  out$gaussian_variance_forest <- keep(gv, list(s2 = gv$sigma2_x_hat_train))
  pr <- bart(X, y_bin, X_test = X_test, num_gfr = 2, num_burnin = 2,
    num_mcmc = 8, rfx_group_ids_train = groups, rfx_group_ids_test = groups_test,
    random_effects_params = list(model_spec = "intercept_only"),
    general_params = list(random_seed = 4L,
      outcome_model = OutcomeModel("binary", "probit")))
  out$probit_rfx <- keep(pr, list(
    prob = predict(pr, X_test, rfx_group_ids = groups_test,
      scale = "probability", terms = "y_hat")))
  cb <- bart(X, y_bin, X_test = X_test, num_gfr = 2, num_mcmc = 6,
    general_params = list(random_seed = 5L,
      outcome_model = OutcomeModel("binary", "cloglog")))
  out$cloglog_binary <- keep(cb)
  co <- bart(X, y_ord, X_test = X_test, num_gfr = 2, num_mcmc = 6,
    general_params = list(random_seed = 6L,
      outcome_model = OutcomeModel("ordinal", "cloglog")))
  out$cloglog_ordinal <- keep(co, list(cut = co$cloglog_cutpoint_samples))
  out
}

if (sys.nframe() == 0L) {
  ref <- outcome_reference_fits()
  ref$platform <- R.version$platform
  ref$generated_from <- "0ac0e7d0 (R/bart.R)"
  saveRDS(ref, file.path("tests", "testthat", "fixtures", "outcome-paths-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
