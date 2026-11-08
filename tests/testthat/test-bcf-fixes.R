# Fixes to bcf() checked against seeded draws captured before them
# (fixtures/make-bcf-fixes-reference.R, cma-rfx 84f5ccbb).

# Fits whose settings no fix touches: must be bitwise identical to the reference.
bcf_fixes_unchanged <- c(
  "gaussian_intercept", "gaussian_adaptive", "gaussian_intercept_adaptive",
  "gaussian_intercept_adaptive_rfx_chains", "gaussian_variance_forest",
  "gaussian_tau_leaf_root", "probit_weighted_intercept_adaptive",
  "gaussian_weighted_intercept", "gaussian_weighted_adaptive",
  "gaussian_weighted_intercept_adaptive"
)
# Fits in the settings a fix targets: must differ from the reference.
bcf_fixes_changed <- c("gaussian_tau_leaf_gfr")

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
