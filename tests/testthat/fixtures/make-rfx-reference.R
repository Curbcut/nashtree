# Generates rfx-unweighted-reference.rds: seeded draws from Gaussian BART and
# BCF models with random effects, captured from the fork at commit 0ac0e7d0,
# before the random-effects conditional learned to honour variance weights.
# The same draws are reproduced bitwise by bd5c5c48 (the nashtree rename),
# the base of the negative-binomial branch, so the reference stands for both.
# test-negbin.R and test-random-effects-weights.R compare current draws with
# these bitwise to prove the unweighted path did not change.
#
# Regenerate only from a build of a pre-change commit (0ac0e7d0 under the
# package name stochtree, or bd5c5c48 as nashtree):
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-rfx-reference.R
library(nashtree)

rfx_reference_fits <- function() {
  set.seed(20261108)
  n <- 200
  p <- 4
  X <- matrix(runif(n * p), ncol = p)
  groups <- sample(1:6, n, replace = TRUE)
  rfx_basis <- cbind(1, runif(n))
  rfx_coefs <- matrix(rnorm(12, 0, 1.5), ncol = 2)
  f_X <- 3 * sin(pi * X[, 1]) + 2 * X[, 2]
  y <- f_X + rowSums(rfx_coefs[groups, ] * rfx_basis) + rnorm(n, 0, 0.7)
  Z <- rbinom(n, 1, 0.3 + 0.4 * X[, 3])
  pi_X <- 0.3 + 0.4 * X[, 3]
  y_bcf <- f_X + 1.5 * Z * X[, 4] + rfx_coefs[groups, 1] + rnorm(n, 0, 0.7)

  out <- list()
  bart_custom <- bart(
    X_train = X, y_train = y,
    rfx_group_ids_train = groups, rfx_basis_train = rfx_basis,
    num_gfr = 3, num_burnin = 5, num_mcmc = 10,
    general_params = list(random_seed = 1234L, num_chains = 2)
  )
  out$bart_custom <- list(
    y_hat_train = bart_custom$y_hat_train,
    rfx_preds_train = bart_custom$rfx_preds_train,
    sigma2 = bart_custom$sigma2_global_samples,
    rfx = bart_custom$rfx_samples$extract_parameter_samples()
  )
  bart_intercept <- bart(
    X_train = X, y_train = y, rfx_group_ids_train = groups,
    num_gfr = 0, num_burnin = 5, num_mcmc = 10,
    general_params = list(random_seed = 99L),
    random_effects_params = list(model_spec = "intercept_only")
  )
  out$bart_intercept <- list(
    y_hat_train = bart_intercept$y_hat_train,
    rfx_preds_train = bart_intercept$rfx_preds_train,
    sigma2 = bart_intercept$sigma2_global_samples,
    rfx = bart_intercept$rfx_samples$extract_parameter_samples()
  )
  bcf_intercept <- bcf(
    X_train = X, Z_train = Z, y_train = y_bcf, propensity_train = pi_X,
    rfx_group_ids_train = groups,
    num_gfr = 2, num_burnin = 3, num_mcmc = 8,
    general_params = list(random_seed = 7L),
    random_effects_params = list(model_spec = "intercept_only")
  )
  out$bcf_intercept <- list(
    y_hat_train = bcf_intercept$y_hat_train,
    tau_hat_train = bcf_intercept$tau_hat_train,
    rfx_preds_train = bcf_intercept$rfx_preds_train,
    rfx = bcf_intercept$rfx_samples$extract_parameter_samples()
  )
  out
}

if (sys.nframe() == 0L) {
  ref <- rfx_reference_fits()
  ref$platform <- R.version$platform
  ref$generated_from <- "0ac0e7d0"
  saveRDS(ref, file.path("tests", "testthat", "fixtures", "rfx-unweighted-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
