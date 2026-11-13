# Generates negbin-reference.rds: seeded draws from the plain negative-binomial
# (count, log link) outcome of bart(), captured on branch cma-nb (da675755)
# before the zero-inflated family was added. test-zinb.R compares current
# draws with these bitwise to prove that `family = "negbin"` (the default) did
# not change.
#
# Regenerate only from a build of da675755:
#   R_LIBS=<scratch lib> Rscript tests/testthat/fixtures/make-negbin-reference.R
library(nashtree)

negbin_reference_fits <- function() {
  set.seed(5151)
  n <- 200
  n_test <- 40
  p <- 4
  X <- matrix(runif((n + n_test) * p), ncol = p)
  offset <- log(runif(n + n_test, 0.5, 2))
  groups <- sample(1:5, n + n_test, replace = TRUE)
  log_mu <- offset + 1 + X[, 1] - X[, 2] + c(-0.3, 0, 0.2, 0.4, -0.1)[groups]
  y <- rnbinom(n + n_test, mu = exp(log_mu), size = 3)
  tr <- seq_len(n)
  te <- n + seq_len(n_test)
  count_model <- OutcomeModel(outcome = "count", link = "log")

  keep <- function(m, extra = list()) {
    c(list(
      y_hat_train = m$y_hat_train,
      y_hat_test = m$y_hat_test,
      dispersion = m$negbin_dispersion_samples,
      sigma2_leaf = m$sigma2_leaf_samples,
      outcome_mean = m$model_params$outcome_mean
    ), extra)
  }
  out <- list()
  a <- bart(X[tr, ], y[tr], X_test = X[te, ],
    offset_train = offset[tr], offset_test = offset[te],
    rfx_group_ids_train = groups[tr], rfx_group_ids_test = groups[te],
    random_effects_params = list(model_spec = "intercept_only"),
    num_gfr = 3, num_burnin = 3, num_mcmc = 8,
    general_params = list(outcome_model = count_model, random_seed = 21L,
      num_chains = 2),
    mean_forest_params = list(num_trees = 30))
  out$offset_rfx_chains <- keep(a, list(
    rate = predict(a, X[te, ], rfx_group_ids = groups[te], offset = offset[te],
      scale = "rate", terms = "y_hat"),
    rfx = a$rfx_preds_train))
  # Weights are passed as 1 / w: the reference was captured when nashtree read
  # observation_weights as variance multipliers (0.4.6.9013 made them precision
  # weights), so 1 / w reproduces its draws.
  b <- bart(X[tr, ], y[tr], X_test = X[te, ], num_gfr = 0, num_burnin = 5,
    num_mcmc = 8, observation_weights = 1 / runif(n, 0.5, 2),
    general_params = list(outcome_model = count_model, random_seed = 22L),
    mean_forest_params = list(num_trees = 20))
  out$root_weighted <- keep(b)
  set.seed(23)
  ppd <- sampleBARTPosteriorPredictive(b, X = X[te, ], num_draws_per_sample = 3)
  out$root_weighted$ppd <- ppd
  out
}

if (sys.nframe() == 0L) {
  ref <- negbin_reference_fits()
  ref$platform <- R.version$platform
  ref$generated_from <- "cma-nb da675755 (R/bart.R)"
  saveRDS(ref, file.path("tests", "testthat", "fixtures", "negbin-reference.rds"))
  cat("wrote reference on", ref$platform, "\n")
}
