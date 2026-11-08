# Internal helpers for the negative-binomial (count, log link) outcome of bart().
#
# Model: y_i ~ NB(mean mu_i, size r), log mu_i = o_i + c + f(x_i) [+ u_g(i)],
# with o_i an optional offset, c a fixed centering constant and f the forest.
#
# Each Gibbs iteration (see bart()):
#   1. Pólya-Gamma step (Polson, Scott and Windle 2013). With
#      psi_i = log mu_i - log r, the NB likelihood is proportional to
#      exp(kappa_i psi_i) / (1 + exp(psi_i))^(y_i + r), kappa_i = (y_i - r) / 2.
#      Drawing omega_i ~ PG(y_i + r, psi_i) makes the pseudo-outcome
#      z_i = kappa_i / omega_i conditionally N(psi_i, 1 / omega_i), so the forest
#      and random effects see a Gaussian regression with variance weights
#      1 / omega_i and the global error variance fixed at 1.
#   2. Forest, leaf scale and random effects: the existing samplers.
#   3. Dispersion: r | y, mu from its exact full conditional, by slice sampling
#      on log r (Neal 2003). See .negbinSampleDispersion() for why the
#      Zhou and Carin (2012) CRT-Gamma update is not used here.

# Log-likelihood of y ~ NB(mean mu, size r), dropping the constant -lgamma(y + 1).
.negbinLogLik <- function(r, y, mu) {
  sum(lgamma(y + r) - lgamma(r) - r * log1p(mu / r) - y * log1p(r / mu))
}

# Full conditional of r given y and mu under a Gamma(shape, rate) prior, drawn by
# a univariate slice sampler on theta = log r with stepping out and shrinkage
# (Neal 2003, Annals of Statistics 31(3), 705-767, Figures 3 and 5).
#
# The forest parametrizes log mu, so mu is held fixed while r is updated. The
# CRT-Gamma update of Zhou and Carin (2012), r | L ~ Gamma(a0 + sum L_i,
# b0 - sum log(1 - p_i)) with p_i = mu_i / (mu_i + r), is the conditional of r
# with p_i held fixed, i.e. for a forest on logit p_i = log mu_i - log r. Used
# with mu held fixed it does not leave p(r | y, mu) invariant: on simulated
# data (n = 2,000, r = 2, mu at the truth) it settled at a mean of 2.44 against
# the exact conditional's 2.05 (sd 0.12), with lag-one autocorrelation 0.92.
# In the log-mu parametrization, mean and dispersion are orthogonal, and an
# exact draw of r mixes well.
.negbinSampleDispersion <- function(
  r,
  y,
  mu,
  prior_shape,
  prior_rate,
  width = 1,
  max_steps = 32L
) {
  log_lower <- log(1e-4)
  log_upper <- log(1e6)
  log_post <- function(theta) {
    if (theta < log_lower || theta > log_upper) {
      return(-Inf)
    }
    r_val <- exp(theta)
    # Gamma(shape, rate) prior on r plus the Jacobian of r = exp(theta)
    value <- .negbinLogLik(r_val, y, mu) +
      prior_shape * theta -
      prior_rate * r_val
    if (is.nan(value)) -Inf else value
  }
  theta0 <- log(r)
  log_post0 <- log_post(theta0)
  if (!is.finite(log_post0)) {
    # theta0 must lie in the slice for the shrinkage loop to terminate
    stop(
      "The negative-binomial dispersion update met a non-finite likelihood ",
      "(r = ", signif(r, 4), ", range of mu = ",
      paste(signif(range(mu), 4), collapse = " to "),
      "); the log rates have probably overflowed"
    )
  }
  log_slice <- log_post0 - rexp(1)
  left <- theta0 - width * runif(1)
  right <- left + width
  steps_left <- floor(max_steps * runif(1))
  steps_right <- max_steps - 1 - steps_left
  while (steps_left > 0 && log_post(left) > log_slice) {
    left <- left - width
    steps_left <- steps_left - 1
  }
  while (steps_right > 0 && log_post(right) > log_slice) {
    right <- right + width
    steps_right <- steps_right - 1
  }
  repeat {
    theta1 <- runif(1, left, right)
    if (log_post(theta1) > log_slice) {
      return(exp(theta1))
    }
    if (theta1 < theta0) left <- theta1 else right <- theta1
  }
}

# Starting value for r: method of moments around the constant-rate fit
# m_i = exp(o_i + c), Var(y_i) = m_i + m_i^2 / r, clamped to [0.1, 1000]. It
# overstates the overdispersion when covariates matter (r too small); the
# sampler moves away from it within a few iterations.
.negbinDispersionInit <- function(y, m) {
  excess <- sum((y - m)^2 - m)
  if (!is.finite(excess) || excess <= 0) {
    return(1000)
  }
  min(max(sum(m^2) / excess, 0.1), 1000)
}

# Validate a log-scale offset vector against the expected length.
.validateNegbinOffset <- function(offset, n, name) {
  if (!is.numeric(offset) || (!is.null(dim(offset)) && NCOL(offset) != 1)) {
    stop(name, " must be a numeric vector")
  }
  offset <- as.numeric(offset)
  if (length(offset) != n) {
    stop("length(", name, ") must equal the number of observations (", n, ")")
  }
  if (any(!is.finite(offset))) {
    stop(name, " must contain only finite values (it is on the log scale)")
  }
  offset
}

# Read the negative-binomial fields of a list of BART JSON objects (one per
# chain or model) into model_params and output.
.loadBartNegbinFieldsFromJson <- function(json_object_list, model_params, output) {
  json_default <- json_object_list[[1]]
  model_params[["has_offset"]] <- json_default$get_boolean("has_offset")
  model_params[["negbin_dispersion_prior_shape"]] <- json_default$get_scalar(
    "negbin_dispersion_prior_shape"
  )
  model_params[["negbin_dispersion_prior_rate"]] <- json_default$get_scalar(
    "negbin_dispersion_prior_rate"
  )
  model_params[["sample_negbin_dispersion"]] <- json_default$get_boolean(
    "sample_negbin_dispersion"
  )
  model_params[["negbin_dispersion_init"]] <- json_default$get_scalar(
    "negbin_dispersion_init"
  )
  dispersion <- unlist(lapply(json_object_list, function(js) {
    js$get_vector("negbin_dispersion_samples", "parameters")
  }))
  output[["negbin_dispersion_samples"]] <- dispersion
  list(model_params = model_params, output = output)
}
