# Internal helpers for count (negative-binomial and zero-inflated
# negative-binomial) outcomes in bcf().
#
# Count part: y_i ~ NB(mu_i, r) with
#   log mu_i = o_i + c + mu_c(x_i, pi_hat_i) + [tau_0 + tau_c(x_i)] b(z_i) [+ rfx],
# b(z) the treatment basis (z, or the adaptive coding b_0 (1 - z) + b_1 z), so
# tau_0 + tau_c(x) is a log rate ratio for a binary treatment and a log rate
# ratio per unit of dose for a continuous one. The Polya-Gamma step of R/negbin.R
# turns each iteration into the weighted Gaussian regression bcf() already
# fits: pseudo-outcome (y_i - r) / (2 omega_i) + log r - o_i with variance
# weights 1 / omega_i and unit error variance. The prognostic forest, the
# basis-regression treatment forest, tau_0, the adaptive coding, the random
# effects, GFR and multiple chains then run on it; tau_0 and the coding
# parameters, which are sampled in R, use the precision weights omega_i.
#
# Zero part (family = "zinb"): pi_i = Phi(c0 + mu_z(x_i, pi_hat_i) + tau_z(x_i) z_i),
# a probit prognostic forest and a probit treatment forest on the structural-
# zero indicators s_i, drawn as in bart() (.zinbSampleStructuralZeros()).

.validateBcfCountInputs <- function(
  y_train,
  X_train,
  X_test,
  offset_train,
  offset_test,
  observation_weights,
  outcome_is_zinb,
  num_trees_variance,
  general_params,
  negbin_dispersion_prior_shape,
  negbin_dispersion_prior_rate,
  negbin_dispersion_init,
  sample_negbin_dispersion,
  log_rate_ratio_prior_sd,
  num_trees_zero_mu,
  num_trees_zero_tau,
  zero_tau_delta_max
) {
  if (num_trees_variance > 0) {
    stop(
      "A variance forest cannot be combined with a count outcome: the negative-binomial variance is mu + mu^2 / r"
    )
  }
  if (isTRUE(general_params[["sample_sigma2_global"]])) {
    stop(
      "sample_sigma2_global = TRUE is not supported with a count outcome: the global error variance of the Polya-Gamma augmented model is fixed at 1"
    )
  }
  y <- as.numeric(y_train)
  if (anyNA(y) || any(!is.finite(y))) {
    stop("A count outcome must not contain missing or infinite values")
  }
  if (any(y < 0) || any(y != round(y))) {
    stop("A count outcome must contain non-negative integers only")
  }
  if (all(y == 0)) {
    stop("A count outcome must contain at least one positive count")
  }
  if (
    !is.numeric(negbin_dispersion_prior_shape) ||
      length(negbin_dispersion_prior_shape) != 1 ||
      !is.finite(negbin_dispersion_prior_shape) ||
      negbin_dispersion_prior_shape <= 0 ||
      !is.numeric(negbin_dispersion_prior_rate) ||
      length(negbin_dispersion_prior_rate) != 1 ||
      !is.finite(negbin_dispersion_prior_rate) ||
      negbin_dispersion_prior_rate < 0
  ) {
    stop(
      "negbin_dispersion_prior_shape must be a positive number and negbin_dispersion_prior_rate a non-negative number"
    )
  }
  if (
    !is.null(negbin_dispersion_init) &&
      (!is.numeric(negbin_dispersion_init) ||
        length(negbin_dispersion_init) != 1 ||
        !is.finite(negbin_dispersion_init) ||
        negbin_dispersion_init <= 0)
  ) {
    stop("negbin_dispersion_init must be a single positive number")
  }
  if (is.null(negbin_dispersion_init) && !isTRUE(sample_negbin_dispersion)) {
    stop("negbin_dispersion_init must be supplied when sample_negbin_dispersion = FALSE")
  }
  if (
    !is.numeric(log_rate_ratio_prior_sd) ||
      length(log_rate_ratio_prior_sd) != 1 ||
      !is.finite(log_rate_ratio_prior_sd) ||
      log_rate_ratio_prior_sd <= 0
  ) {
    stop("treatment_effect_forest_params$log_rate_ratio_prior_sd must be a single positive number")
  }
  if (outcome_is_zinb) {
    if (!any(y == 0)) {
      stop(
        "A zero-inflated count outcome needs at least one zero in y_train; use family = 'negbin'"
      )
    }
    if (!is.null(observation_weights)) {
      stop(
        "observation_weights are not supported with a zero-inflated count outcome"
      )
    }
    for (num_trees in c(num_trees_zero_mu, num_trees_zero_tau)) {
      if (!is.numeric(num_trees) || length(num_trees) != 1 || num_trees < 1) {
        stop(
          "A zero-inflated count outcome requires zero-part forests with num_trees >= 1"
        )
      }
    }
    if (
      !is.numeric(zero_tau_delta_max) ||
        length(zero_tau_delta_max) != 1 ||
        !(zero_tau_delta_max > 0)
    ) {
      stop("zero_treatment_effect_forest_params$delta_max must be a single positive number")
    }
  }
  offset_train <- if (is.null(offset_train)) {
    rep(0, NROW(X_train))
  } else {
    .validateNegbinOffset(offset_train, NROW(X_train), "offset_train")
  }
  if (!is.null(offset_test)) {
    if (is.null(X_test)) {
      stop("offset_test was supplied without X_test")
    }
    offset_test <- .validateNegbinOffset(offset_test, NROW(X_test), "offset_test")
  } else if (!is.null(X_test)) {
    if (any(offset_train != 0)) {
      warning(
        "offset_train was supplied but offset_test was not; test-set predictions use an offset of 0 (unit exposure)"
      )
    }
    offset_test <- rep(0, NROW(X_test))
  }
  list(offset_train = offset_train, offset_test = offset_test)
}

# Read the zero-inflated fields of a list of BCF JSON objects (one per chain or
# model) into model_params.
.loadBcfZinbFieldsFromJson <- function(json_object_list, model_params) {
  json_default <- json_object_list[[1]]
  for (field in c(
    "zero_outcome_mean",
    "zinb_zero_probability_init",
    "num_trees_zero_mu",
    "num_trees_zero_tau",
    "sigma2_leaf_zero_mu",
    "sigma2_leaf_zero_tau"
  )) {
    model_params[[field]] <- json_default$get_scalar(field)
  }
  model_params
}

# Count fields of a BCF model: dispersion draws and priors, offset flag, the
# family and, for "zinb", the zero-part forests and fields. Shared by the three
# JSON loaders; `combined` selects the multi-object forest loader.
.loadBcfCountFieldsFromJson <- function(json_object_list, model_params, output, combined) {
  negbin_fields <- .loadBartNegbinFieldsFromJson(
    json_object_list,
    model_params,
    output
  )
  model_params <- negbin_fields$model_params
  output <- negbin_fields$output
  if (identical(model_params[["outcome_model"]]$family, "zinb")) {
    model_params <- .loadBcfZinbFieldsFromJson(json_object_list, model_params)
    load_forest <- function(label) {
      if (combined) {
        loadForestContainerCombinedJson(json_object_list, label)
      } else {
        loadForestContainerJson(json_object_list[[1]], label)
      }
    }
    # Forests are written mu, tau, [variance,] zero mu, zero tau; a count model
    # has no variance forest
    output[["forests_zero_mu"]] <- load_forest("forest_2")
    output[["forests_zero_tau"]] <- load_forest("forest_3")
  }
  list(model_params = model_params, output = output)
}
