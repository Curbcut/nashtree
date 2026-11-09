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

# ---------------------------------------------------------------------------
# Zero-inflated negative binomial (family = "zinb")
#
# Model: with probability pi_i = Phi(eta0_i), eta0_i = c0 + g(x_i), y_i is a
# structural zero; otherwise y_i ~ NB(mu_i, r) as above. With the structural-
# zero indicator s_i as a latent variable, one Gibbs iteration of bart() is:
#   1. s_i | y_i = 0, mu_i, pi_i, r ~ Bernoulli(pi_i / (pi_i + (1 - pi_i) NB(0 | mu_i, r)));
#      s_i = 0 when y_i > 0 (.zinbSampleStructuralZeros()).
#   2. The zero forest g on s by the Albert and Chib (1993) probit step.
#   3. The count forest on the rows with s_i = 0 by the Polya-Gamma step above;
#      rows with s_i = 1 carry an infinite variance weight, so 1 / w_i = 0 and
#      they drop out of every sufficient statistic of the forest and of the
#      random effects.
#   4. r | y, mu over the rows with s_i = 0 (.negbinSampleDispersion()).
# ---------------------------------------------------------------------------

# log NB(0 | mean mu, size r) = -r log(1 + mu / r)
.negbinLogZeroProb <- function(mu, r) {
  -r * log1p(mu / r)
}

# Draw the structural-zero indicators. eta_count is log(mu) and eta_zero the
# probit index of pi, both full linear predictors. The Bernoulli probability is
# formed on the log-odds scale so that pi near 0 or 1 and large mu stay finite.
.zinbSampleStructuralZeros <- function(y, eta_count, eta_zero, r) {
  s <- integer(length(y))
  zero <- which(y == 0)
  if (length(zero) == 0) {
    return(s)
  }
  log_odds <- pnorm(eta_zero[zero], log.p = TRUE) -
    pnorm(eta_zero[zero], lower.tail = FALSE, log.p = TRUE) -
    .negbinLogZeroProb(exp(eta_count[zero]), r)
  s[zero] <- as.integer(runif(length(zero)) < plogis(log_odds))
  s
}

# Constant-rate ZINB fit used to initialise the sampler: the excess-zero
# fraction over the constant-rate NB fit (method of moments, as for "negbin")
# starts a maximum-likelihood fit of (pi, log rate, r) with no covariates. The
# NB fit alone understates the zero mass (the extra zeros inflate the variance,
# so the moment estimate of r is too small and the NB zero probability too
# large), which is why the three parameters are refined jointly. Returns the
# zero probability, the log rate per unit exposure, the dispersion and the
# posterior probability that each zero is structural under that fit.
.zinbConstantFit <- function(y, offset) {
  exposure <- exp(offset)
  log_rate_nb <- log(sum(y) / sum(exposure))
  m <- exp(offset + log_rate_nb)
  r_nb <- .negbinDispersionInit(y, m)
  p0_obs <- mean(y == 0)
  p0_nb <- mean(exp(.negbinLogZeroProb(m, r_nb)))
  excess <- if (p0_nb < 1) max((p0_obs - p0_nb) / (1 - p0_nb), 0) else 0
  pi_start <- min(max(excess, 0.01), 0.9)
  is_zero <- y == 0
  y_pos <- y[!is_zero]
  negloglik <- function(par) {
    pi <- plogis(par[1])
    r <- exp(par[3])
    mu <- exp(offset + par[2])
    log_p0 <- .negbinLogZeroProb(mu[is_zero], r)
    # log(pi + (1 - pi) p0) via log-sum-exp
    a <- log(pi)
    b <- log1p(-pi) + log_p0
    hi <- pmax(a, b)
    ll_zero <- hi + log(exp(a - hi) + exp(b - hi))
    mu_pos <- mu[!is_zero]
    ll_pos <- log1p(-pi) + lgamma(y_pos + r) - lgamma(r) - lgamma(y_pos + 1) -
      r * log1p(mu_pos / r) - y_pos * log1p(r / mu_pos)
    value <- -(sum(ll_zero) + sum(ll_pos))
    if (is.finite(value)) value else 1e300
  }
  start <- c(qlogis(pi_start), log_rate_nb - log1p(-pi_start), log(r_nb))
  fit <- tryCatch(
    optim(start, negloglik, method = "BFGS", control = list(maxit = 200)),
    error = function(e) NULL
  )
  par <- if (!is.null(fit) && all(is.finite(fit$par))) fit$par else start
  zero_probability <- plogis(par[1])
  log_rate <- par[2]
  dispersion <- min(max(exp(par[3]), 0.1), 1000)
  log_p0 <- .negbinLogZeroProb(exp(offset + log_rate), dispersion)
  structural_probability <- numeric(length(y))
  structural_probability[is_zero] <- plogis(
    qlogis(zero_probability) - log_p0[is_zero]
  )
  list(
    zero_probability = zero_probability,
    log_rate = log_rate,
    dispersion = dispersion,
    excess_zero_fraction = excess,
    structural_probability = structural_probability
  )
}

# JSON label of the zero forest: forests are written in the order mean,
# variance, zero, and a count model never has a variance forest.
.zeroForestJsonLabel <- function(include_mean_forest, include_variance_forest) {
  paste0(
    "forest_",
    as.integer(include_mean_forest) + as.integer(include_variance_forest)
  )
}

# Read the ZINB fields of a list of BART JSON objects (one per chain or model)
# into model_params and output.
.loadBartZinbFieldsFromJson <- function(json_object_list, model_params, output) {
  json_default <- json_object_list[[1]]
  model_params[["zero_outcome_mean"]] <- json_default$get_scalar("zero_outcome_mean")
  model_params[["zinb_zero_probability_init"]] <- json_default$get_scalar(
    "zinb_zero_probability_init"
  )
  model_params[["num_trees_zero"]] <- json_default$get_scalar("num_trees_zero")
  model_params[["sample_sigma2_leaf_zero"]] <- json_default$get_boolean(
    "sample_sigma2_leaf_zero"
  )
  model_params[["sigma2_leaf_zero_init"]] <- json_default$get_scalar(
    "sigma2_leaf_zero_init"
  )
  sampled <- vapply(json_object_list, function(js) {
    js$get_boolean("sample_sigma2_leaf_zero")
  }, logical(1))
  if (length(unique(sampled)) > 1) {
    stop(
      "Cannot combine zero-inflated models that differ in zero_forest_params$sample_sigma2_leaf"
    )
  }
  if (isTRUE(model_params[["sample_sigma2_leaf_zero"]])) {
    output[["sigma2_leaf_zero_samples"]] <- unlist(lapply(
      json_object_list,
      function(js) js$get_vector("sigma2_leaf_zero_samples", "parameters")
    ))
  }
  zero_rfx <- .readZeroRfxFlags(json_object_list)
  model_params[["has_zero_rfx"]] <- zero_rfx
  if (zero_rfx) {
    output[["zero_rfx_samples"]] <- .loadZeroRfxSamplesFromJson(json_object_list)
  }
  list(model_params = model_params, output = output)
}

# Whether the zero part of a list of zero-inflated JSON objects (one per chain
# or model) carries random effects. Models saved before zero-part random
# effects existed have no has_zero_rfx field and none.
.readZeroRfxFlags <- function(json_object_list) {
  flags <- vapply(json_object_list, function(js) {
    json_contains_field_cpp(js$json_ptr, "has_zero_rfx") &&
      js$get_boolean("has_zero_rfx")
  }, logical(1))
  if (length(unique(flags)) > 1) {
    stop(
      "Cannot combine zero-inflated models with and without zero-part random effects"
    )
  }
  flags[1]
}

# Zero-part random effects of a zero-inflated BART or BCF model: written after
# the count part's random effects, so they are random_effect_container_1
.loadZeroRfxSamplesFromJson <- function(json_object_list) {
  loadRandomEffectSamplesCombinedJson(json_object_list, 1)
}

# Prior of the zero part's random effects in a zero-inflated model. The keys
# of `zero_params` (zero_random_effects_params) override the count part's
# `rfx_params` (random_effects_params after defaults); unset keys keep the
# count part's values, so an empty list gives the count part's prior exactly.
# The groups, basis (model_spec) and unseen-group policy are the count part's
# and cannot be set here. Returns the same initial values and priors the count
# part computes in bart() and bcf().
.zinbZeroRfxPriorKeys <- c(
  "working_parameter_prior_mean",
  "group_parameter_prior_mean",
  "working_parameter_prior_cov",
  "group_parameter_prior_cov",
  "variance_prior_shape",
  "variance_prior_scale"
)

.zinbCheckZeroRfxParams <- function(zero_params) {
  if (is.null(zero_params)) {
    return(invisible(NULL))
  }
  if (!is.list(zero_params)) {
    stop("zero_random_effects_params must be a list")
  }
  if (length(zero_params) > 0) {
    keys <- names(zero_params)
    if (is.null(keys) || any(!nzchar(keys))) {
      stop("zero_random_effects_params must be a named list")
    }
    bad <- setdiff(keys, .zinbZeroRfxPriorKeys)
    if (length(bad) > 0) {
      stop(
        "zero_random_effects_params accepts only ",
        paste(.zinbZeroRfxPriorKeys, collapse = ", "),
        " (the zero part shares the count part's groups, basis and unseen-group policy); got: ",
        paste(bad, collapse = ", ")
      )
    }
  }
  for (key in c("variance_prior_shape", "variance_prior_scale")) {
    val <- zero_params[[key]]
    if (
      !is.null(val) &&
        (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0)
    ) {
      stop(sprintf(
        "zero_random_effects_params$%s must be a single positive number",
        key
      ))
    }
  }
  invisible(NULL)
}

.zinbZeroRfxPrior <- function(
  rfx_params,
  zero_params,
  num_components,
  num_groups
) {
  p <- preprocessParams(rfx_params[.zinbZeroRfxPriorKeys], zero_params)
  if (is.null(p$working_parameter_prior_mean)) {
    alpha_init <- rep(0, num_components)
  } else {
    alpha_init <- expand_dims_1d(p$working_parameter_prior_mean, num_components)
  }
  if (is.null(p$group_parameter_prior_mean)) {
    xi_init <- matrix(rep(alpha_init, num_groups), num_components, num_groups)
  } else {
    xi_init <- expand_dims_2d(
      p$group_parameter_prior_mean,
      num_components,
      num_groups
    )
  }
  if (is.null(p$working_parameter_prior_cov)) {
    sigma_alpha_init <- diag(1, num_components, num_components)
  } else {
    sigma_alpha_init <- expand_dims_2d_diag(
      p$working_parameter_prior_cov,
      num_components
    )
  }
  if (is.null(p$group_parameter_prior_cov)) {
    sigma_xi_init <- diag(1, num_components, num_components)
  } else {
    sigma_xi_init <- expand_dims_2d_diag(
      p$group_parameter_prior_cov,
      num_components
    )
  }
  list(
    alpha_init = alpha_init,
    xi_init = xi_init,
    sigma_alpha_init = sigma_alpha_init,
    sigma_xi_init = sigma_xi_init,
    sigma_xi_shape = p$variance_prior_shape,
    sigma_xi_scale = p$variance_prior_scale
  )
}

# Column indices of X_train a forest may split on, from keep / drop lists of
# names or indices (keep takes precedence), as bart() does for its mean and
# variance forests.
.resolveForestVariableSubset <- function(X_train, keep_vars, drop_vars, keep_label, drop_label) {
  p <- ncol(X_train)
  resolve <- function(vars, label) {
    if (is.character(vars)) {
      if (!all(vars %in% names(X_train))) {
        stop(label, " includes some variable names that are not in X_train")
      }
      return(unname(which(names(X_train) %in% vars)))
    }
    if (any(vars > p)) {
      stop(label, " includes some variable indices that exceed the number of columns in X_train")
    }
    if (any(vars < 0)) {
      stop(label, " includes some negative variable indices")
    }
    vars
  }
  if (!is.null(keep_vars)) {
    return(resolve(keep_vars, keep_label))
  }
  if (!is.null(drop_vars)) {
    drop <- resolve(drop_vars, drop_label)
    return((1:p)[!(1:p %in% drop)])
  }
  1:p
}

# Albert and Chib (1993) latent draw: z_i ~ N(eta_i, 1) truncated to z_i > 0
# when s_i = 1 and z_i <= 0 when s_i = 0, by inversion on the log scale so that
# |eta_i| in the tens stays finite.
.probitLatentDraw <- function(eta, s) {
  u <- runif(length(eta))
  e <- numeric(length(eta))
  pos <- s == 1L
  # s = 0: e = z - eta < -eta
  e[!pos] <- qnorm(
    log(u[!pos]) + pnorm(-eta[!pos], log.p = TRUE),
    log.p = TRUE
  )
  # s = 1: e > -eta, i.e. -e < eta
  e[pos] <- -qnorm(log(u[pos]) + pnorm(eta[pos], log.p = TRUE), log.p = TRUE)
  eta + e
}

# Count family stored under outcome_model in a BART / BCF JSON object. Models
# saved before the zero-inflated family existed carry no family field and are
# plain negative-binomial models.
.readOutcomeFamilyFromJson <- function(json_object, outcome) {
  if (!identical(outcome, "count")) {
    return(NULL)
  }
  if (json_contains_field_subfolder_cpp(json_object$json_ptr, "outcome_model", "family")) {
    json_object$get_string("family", "outcome_model")
  } else {
    "negbin"
  }
}
