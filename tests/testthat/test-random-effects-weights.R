# Variance weights in the random-effects conditional (src/random_effects.cpp).
#
# A RandomEffectsDataset row i with variance weight v_i has residual variance
# global_variance * v_i, so it enters the group and working-parameter
# conditionals with precision weight 1 / v_i. Without weights the sampler must
# be bitwise unchanged; the reference draws were captured from 0ac0e7d0 by
# fixtures/make-rfx-reference.R.

reference_platform_matches <- function(ref) {
  # std::mt19937 is portable but the draws also depend on the floating-point
  # environment (libm, Eigen code paths), so bitwise identity is only
  # meaningful on the platform that produced the reference.
  identical(ref$platform, R.version$platform)
}

test_that("Unweighted random-effects sampling is bitwise identical to 0ac0e7d0", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "rfx-unweighted-reference.rds"))
  skip_if_not(
    reference_platform_matches(ref),
    paste("reference draws were generated on", ref$platform)
  )
  source(test_path("fixtures", "make-rfx-reference.R"), local = TRUE)
  current <- rfx_reference_fits()
  for (name in c("bart_custom", "bart_intercept", "bcf_intercept")) {
    expect_identical(current[[name]], ref[[name]], label = name)
  }
})

run_rfx_sampler <- function(y, groups, variance_weights, num_iter, seed) {
  n <- length(y)
  num_groups <- length(unique(groups))
  basis <- matrix(1, n, 1)
  rfx_dataset <- createRandomEffectsDataset(groups, basis, variance_weights)
  outcome <- createOutcome(y)
  tracker <- createRandomEffectsTracker(groups)
  model <- createRandomEffectsModel(1, num_groups)
  model$set_working_parameter(1)
  model$set_group_parameters(matrix(0, 1, num_groups))
  model$set_working_parameter_cov(matrix(1))
  model$set_group_parameter_cov(matrix(1))
  model$set_variance_prior_shape(1)
  model$set_variance_prior_scale(1)
  samples <- createRandomEffectSamples(1, num_groups, tracker)
  rng <- createCppRNG(seed)
  for (i in seq_len(num_iter)) {
    model$sample_random_effect(rfx_dataset, outcome, tracker, samples, TRUE, 1.0, rng)
  }
  beta <- samples$extract_parameter_samples()$beta_samples
  beta[, -(1:200), drop = FALSE] # drop burn-in
}

test_that("Weighted random effects recover group means under heteroskedastic noise", {
  skip_on_cran()
  set.seed(11)
  num_groups <- 20
  per_group <- 60
  groups <- rep(seq_len(num_groups), each = per_group)
  n <- length(groups)
  group_means <- rnorm(num_groups, 0, 2)
  # Half the rows are precise (sd 0.1), half nearly uninformative (sd 25)
  noisy <- rep(rep(c(FALSE, TRUE), each = per_group / 2), num_groups)
  v <- ifelse(noisy, 25^2, 0.1^2)
  y <- group_means[groups] + rnorm(n, 0, sqrt(v))

  weighted <- run_rfx_sampler(y, groups, v, 1200, 1)
  post_mean <- rowMeans(weighted)
  post_sd <- apply(weighted, 1, sd)

  # Under a flat prior the group posterior is N(sum(y / v) / sum(1 / v), 1 / sum(1 / v));
  # the N(0, sigma_xi) prior is negligible next to a precision of about 3,000.
  precision <- tapply(1 / v, groups, sum)
  wmean <- tapply(y / v, groups, sum) / precision
  expect_lt(max(abs(post_mean - wmean)), 0.01)
  expect_lt(max(abs(post_mean - group_means)), 6 / sqrt(min(precision)))
  sd_ratio <- post_sd / sqrt(1 / precision)
  expect_true(all(sd_ratio > 0.7 & sd_ratio < 1.4))

  # Unit weights reproduce the unweighted conditional, which treats every row
  # alike and so is pulled about by the noisy half.
  unit <- run_rfx_sampler(y, groups, rep(1, n), 1200, 1)
  unweighted <- run_rfx_sampler(y, groups, NULL, 1200, 1)
  expect_equal(unit, unweighted, tolerance = 1e-12)
  rmse_weighted <- sqrt(mean((post_mean - group_means)^2))
  rmse_unweighted <- sqrt(mean((rowMeans(unweighted) - group_means)^2))
  expect_lt(rmse_weighted * 20, rmse_unweighted)
})
