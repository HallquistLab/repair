# Shared operations used by the three M-MOM estimation steps.

# Backend and correction-control normalization.
mmom_backend <- function(backend, supported = c("auto", "r", "rcpp")) {
  if (length(backend) > 1L) backend <- backend[[1L]]
  backend <- tolower(backend)
  backend <- rlang::arg_match(backend, supported)
  if (identical(backend, "auto")) {
    return(if (exists("mmom_normalize_err_x_cpp")) "rcpp" else "r")
  }
  backend
}

# Parse a control object (for preliminary and final moment corrections)
mmom_control <- function(control, npredictors, label) {
  if (!is.list(control)) {
    repair_error("{label} must be a list.")
  }
  correction <- control$correction %||% "minimum"
  if (!correction %in% c("none", "minimum", "step_up", "custom")) {
    repair_error("{label}$correction must be one of `none`, `minimum`, `step_up`, or `custom`.")
  }
  
  if (identical(correction, "none") || identical(correction, "step_up")) {
    alpha <- NULL
  } else if (identical(correction, "minimum")) {
    alpha <- npredictors + 2L
  } else if (identical(correction, "custom")) {
    alpha <- control$alpha %||% (npredictors + 2L)
    if (!is.numeric(alpha) || length(alpha) != 1L ||
        !is.finite(alpha) || alpha <= 0) {
      repair_error("{label}$alpha must be a finite, positive scalar.")
    }
    if (alpha < npredictors + 2L) {
      repair_error("{label}$alpha must be greater than or equal to `npredictors + 2`.")
    }
  }
  list(
    correction = correction,
    alpha = alpha
  )
}

mmom_prepare <- function(
  data_list, preliminary_moment_control, final_moment_control, backend
) {
  # Coerce the supported row-wise and matrix-style input forms.
  data <- mmom_coerce_data(data_list)
  outcome <- as.numeric(data$outcome)
  predictors <- as.matrix(data$predictors)
  storage.mode(predictors) <- "double"
  nobs <- length(outcome)
  npredictors <- ncol(predictors)
  if (nrow(predictors) != nobs || is.null(npredictors) ||
      npredictors < 1L) {
    repair_error(
      "`data_list` must provide an nobs by npredictors predictor matrix."
    )
  }
  if (nobs < max(8L, npredictors + 2L)) {
    repair_error(
      "M-MOM requires at least max(8, npredictors + 2) observations."
    )
  }
  if (any(!is.finite(outcome)) || any(!is.finite(predictors))) {
    repair_error("M-MOM inputs must be finite.")
  }
  if (is.null(colnames(predictors))) {
    colnames(predictors) <- mmom_predictor_names(npredictors)
  }
  if (anyNA(colnames(predictors)) || any(!nzchar(colnames(predictors))) ||
    anyDuplicated(c("(Intercept)", colnames(predictors)))) {
    repair_error(
      paste(
        "Predictor names must be nonempty, unique,",
        "and distinct from `(Intercept)`."
      )
    )
  }
  err_x <- mmom_normalize_err_x(data$err_x, nobs, npredictors, backend)
  err_y <- data$err_y
  if (is.null(err_y)) err_y <- rep(0, nobs)
  if (!is.numeric(err_y) || length(err_y) != nobs ||
      any(!is.finite(err_y)) || any(err_y < 0)) {
    repair_error(
      "`err_y` must be a finite, nonnegative vector of length nobs."
    )
  }
  err_xy <- mmom_normalize_err_xy(data$err_xy, nobs, npredictors)

  # Validate the predictor and joint outcome/predictor error blocks.
  mmom_validate_joint_cov(err_x, err_y, err_xy, backend)
  design <- cbind(`(Intercept)` = 1, predictors)
  list(
    outcome = outcome,
    predictors = predictors,
    err_x = err_x,
    err_y = err_y,
    err_xy = err_xy,
    design = design,
    design_index = seq_len(npredictors + 1L) + 1L,
    predictor_index = seq.int(2L, npredictors + 1L),
    nobs = nobs,
    npredictors = npredictors,
    ncoef = npredictors + 1L,
    backend = backend,
    preliminary_moment_control = mmom_control(
      preliminary_moment_control, npredictors, "preliminary_moment"
    ),
    final_moment_control = mmom_control(
      final_moment_control, npredictors, "final_moment"
    )
  )
}

mmom_predictor_names <- function(npredictors) {
  letters <- c(letters, as.vector(outer(letters, letters, paste0)))
  paste0("predictor_", letters[seq_len(npredictors)])
}

mmom_coerce_data <- function(data_list) {
  if (!is.list(data_list)) repair_error("`data_list` must be a list.")
  if (!is.null(data_list$outcome) &&
      !is.null(data_list$predictors)) {
    return(data_list)
  }
  if (!length(data_list) || !all(vapply(data_list, is.list, logical(1)))) {
    repair_error("`data_list` must contain observation lists.")
  }
  get_field <- function(observation, fields, default = NULL) {
    matches <- fields[fields %in% names(observation)]
    if (!length(matches)) default else observation[[matches[[1L]]]]
  }
  predictors <- lapply(
    data_list,
    function(observation) get_field(observation, c("predictors", "x"))
  )
  if (any(vapply(predictors, is.null, logical(1)))) {
    repair_error("Each observation needs `predictors`.")
  }
  list(
    outcome = vapply(
      data_list,
      function(observation) {
        as.numeric(get_field(observation, c("outcome", "y")))
      },
      numeric(1)
    ),
    predictors = do.call(rbind, predictors),
    err_x = lapply(
      data_list,
      function(observation) get_field(observation, "err_x")
    ),
    err_y = vapply(
      data_list,
      function(observation) {
        as.numeric(get_field(
          observation, c("err_y"), 0
        ))
      },
      numeric(1)
    ),
    err_xy = {
      values <- lapply(
        data_list,
        function(observation) {
          get_field(observation, "err_xy")
        }
      )
      if (all(vapply(values, is.null, logical(1)))) NULL else values
    }
  )
}

# Covariance coercion and positive-semidefinite checks.
mmom_normalize_err_x <- function(
  err_x, nobs, npredictors, backend, tolerance = 1e-8
) {
  out <- mmom_coerce_err_x(err_x, nobs, npredictors)
  if (backend == "rcpp") {
    normalized <- mmom_cpp_call(
      mmom_normalize_err_x_cpp, out, nobs, npredictors, tolerance
    )
    if (normalized$bad_row > 0L) {
      repair_error(
        "`err_x` is not positive semidefinite for row {normalized$bad_row}.",
      )
    }
    return(normalized$err_x)
  }
  for (row_index in seq_len(nobs)) {
    out[, , row_index] <- mmom_symmetrize(out[, , row_index])
    diagnostics <- mmom_matrix_diagnostics(out[, , row_index])
    if (!is.finite(diagnostics$min_eigen) ||
        diagnostics$min_eigen < -tolerance * max(
          1, abs(diagnostics$max_eigen)
        )) {
      repair_error(
        "`err_x` is not positive semidefinite for row {row_index}."
      )
    }
  }
  out
}

mmom_coerce_err_x <- function(err_x, nobs, npredictors) {
  if (is.matrix(err_x)) {
    if (!identical(dim(err_x), c(npredictors, npredictors))) {
      repair_error(
        "A common `err_x` matrix must be square with one row per predictor."
      )
    }
    out <- array(
      rep(as.numeric(err_x), nobs),
      c(npredictors, npredictors, nobs)
    )
  } else if (is.array(err_x) && length(dim(err_x)) == 3L) {
    dimensions <- dim(err_x)
    if (identical(dimensions, c(npredictors, npredictors, nobs))) {
      out <- err_x
    } else if (identical(dimensions, c(nobs, npredictors, npredictors))) {
      out <- aperm(err_x, c(2, 3, 1))
    } else {
      repair_error(
        paste(
          "`err_x` must be npredictors by npredictors by nobs or",
          "nobs by npredictors by npredictors."
        )
      )
    }
  } else if (is.list(err_x)) {
    if (length(err_x) != nobs) {
      repair_error("An `err_x` list must have length nobs.")
    }
    out <- array(NA_real_, c(npredictors, npredictors, nobs))
    for (row_index in seq_len(nobs)) {
      block <- as.matrix(err_x[[row_index]])
      if (!identical(dim(block), c(npredictors, npredictors))) {
        repair_error(
          paste(
            "Every `err_x` block must be square with",
            "one row per predictor."
          )
        )
      }
      out[, , row_index] <- block
    }
  } else {
    repair_error(
      paste(
        "`err_x` must be a matrix, three-dimensional array,",
        "or list of matrices."
      )
    )
  }
  if (any(!is.finite(out))) {
    repair_error("`err_x` contains non-finite values.")
  }
  out
}

mmom_normalize_err_xy <- function(err_xy, nobs, npredictors) {
  if (is.null(err_xy)) {
    return(matrix(
      0,
      nrow = nobs,
      ncol = npredictors
    ))
  }
  mmom_coerce_err_xy(
    err_xy = err_xy,
    nobs = nobs,
    npredictors = npredictors
  )
}

mmom_coerce_err_xy <- function(err_xy, nobs, npredictors) {
  if (is.numeric(err_xy) && is.null(dim(err_xy))) {
    if (length(err_xy) != npredictors) {
      repair_error(
        "`err_xy` must have length npredictors."
      )
    }
    out <- matrix(
      rep(err_xy, each = nobs),
      nrow = nobs,
      ncol = npredictors
    )
  } else if (is.matrix(err_xy)) {
    dimensions <- dim(err_xy)
    if (identical(dimensions, c(npredictors, 1L))) {
      out <- matrix(
        rep(as.numeric(err_xy), each = nobs),
        nrow = nobs,
        ncol = npredictors
      )
    } else if (identical(dimensions, c(nobs, npredictors))) {
      out <- err_xy
    } else {
      repair_error(
        paste(
          "`err_xy` must be length npredictors,",
          "npredictors by 1, or nobs by npredictors."
        )
      )
    }
  } else if (is.list(err_xy)) {
    if (length(err_xy) != nobs) {
      repair_error(
        "An `err_xy` list must have length nobs."
      )
    }
    out <- matrix(
      NA_real_,
      nrow = nobs,
      ncol = npredictors
    )
    for (row_index in seq_len(nobs)) {
      block <- err_xy[[row_index]]
      if (is.matrix(block)) {
        if (!identical(dim(block), c(npredictors, 1L))) {
          repair_error(
            paste(
              "Every `err_xy` matrix must have dimensions",
              "npredictors by 1."
            )
          )
        }
        block <- as.numeric(block)
      } else {
        block <- as.numeric(block)
        if (length(block) != npredictors) {
          repair_error(
            paste(
              "Every `err_xy` element must have length",
              "npredictors."
            )
          )
        }
      }
      out[row_index, ] <- block
    }
  } else {
    repair_error(
      paste(
        "`err_xy` must be a numeric vector, matrix,",
        "or list of vectors or matrices."
      )
    )
  }
  if (any(!is.finite(out))) {
    repair_error("`err_xy` contains non-finite values.")
  }
  out
}


mmom_validate_joint_cov <- function(
  err_x, err_y, err_xy, backend, tolerance = 1e-8
) {
  if (backend == "rcpp") {
    bad_row <- mmom_cpp_call(
      mmom_joint_err_cov_bad_row_cpp,
      err_x,
      err_y,
      err_xy,
      tolerance
    )
    if (bad_row > 0L) {
      repair_error("The joint measurement-error covariance is not PSD for row {bad_row}.")
    }
    return(invisible(NULL))
  }
  for (row_index in seq_along(err_y)) {
    joint <- rbind(
      c(err_y[[row_index]], err_xy[row_index, ]),
      cbind(err_xy[row_index, ], err_x[, , row_index])
    )
    diagnostics <- mmom_matrix_diagnostics(joint)
    if (!is.finite(diagnostics$min_eigen) ||
        diagnostics$min_eigen < -tolerance * max(
          1, abs(diagnostics$max_eigen)
        )) {
      repair_error(
        "The joint measurement-error covariance is not PSD for row {row_index}.")
    }
  }
}

# Small matrix helpers and native-call forwarding.
mmom_cpp_call <- function(fun, ...) {
  tryCatch(
    fun(...),
    error = function(error) repair_error(conditionMessage(error))
  )
}

mmom_symmetrize <- function(matrix_value) (matrix_value + t(matrix_value)) / 2

mmom_matrix_diagnostics <- function(matrix_value) {
  values <- tryCatch(
    eigen(
      mmom_symmetrize(as.matrix(matrix_value)),
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    error = function(error) NA_real_
  )
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(list(
      min_eigen = NA_real_,
      max_eigen = NA_real_,
      condition_number = Inf
    ))
  }
  min_eigen <- min(values)
  max_eigen <- max(values)
  list(
    min_eigen = min_eigen,
    max_eigen = max_eigen,
    condition_number = if (min_eigen > 0) max_eigen / min_eigen else Inf
  )
}

mmom_smallest_det_root <- function(left, right) {
  if (!requireNamespace("geigen", quietly = TRUE)) {
    repair_error("The `geigen` package is required for M-MOM estimation.")
  }
  # Roots solve det(a - lambda * b) = 0. The error-free intercept makes b
  # singular, so the solver must permit positive-semidefinite b; its symmetric
  # mode requires positive definiteness. Ignore infinite roots and numerical
  # imaginary noise. NA means no usable finite root and selects the alpha-only
  # correction branch in the caller, matching the established R estimator.
  roots <- tryCatch(
    geigen::geigen(
      left, right, symmetric = FALSE, only.values = TRUE
    )$values,
    error = function(error) NULL
  )
  if (is.null(roots)) {
    return(NA_real_)
  }

  # Keep numerically real generalized-eigenvalue roots only.
  if (is.complex(roots)) {
    scale <- suppressWarnings(max(abs(roots), na.rm = TRUE))
    if (!is.finite(scale)) scale <- 1
    roots <- ifelse(
      abs(Im(roots)) <= sqrt(.Machine$double.eps) * max(1, scale),
      Re(roots),
      NA_real_
    )
  }
  roots <- roots[is.finite(roots)]
  if (length(roots)) min(roots) else NA_real_
}

mmom_weighted_err_cov_sum <- function(err_x, weights, backend) {
  if (backend == "rcpp") {
    return(mmom_symmetrize(mmom_cpp_call(
      mmom_weighted_err_cov_sum_cpp, err_x, weights
    )))
  }
  result <- matrix(0, dim(err_x)[1L], dim(err_x)[2L])
  for (row_index in seq_len(dim(err_x)[3L])) {
    result <- result + weights[[row_index]] * err_x[, , row_index]
  }
  mmom_symmetrize(result)
}

# Moment construction, residual variance, and sandwich components.
mmom_error_moment <- function(err_x_sum, err_y_sum, err_xy_sum) {
  nvars <- nrow(err_x_sum) + 1L
  moment <- matrix(0, nvars + 1L, nvars + 1L)
  moment[1, 1] <- err_y_sum
  moment[1, 2:(nvars + 1L)] <- err_xy_sum
  moment[2:(nvars + 1L), 1] <- err_xy_sum

  # The intercept is measured without error.
  moment[2:(nvars + 1L), 2:(nvars + 1L)] <- {
    predictor_block <- matrix(0, nvars, nvars)
    predictor_block[-1L, -1L] <- err_x_sum
    predictor_block
  }
  moment
}

mmom_composite_err_y <- function(inputs, coef_predictors) {
  if (inputs$backend == "rcpp") {
    # The native path replaces the row-wise covariance loop.
    cov_times_coef <- mmom_cpp_call(
      mmom_err_cov_matvec_cpp, inputs$err_x, coef_predictors
    )
    repeated_coef <- rep(coef_predictors, each = inputs$nobs)
    return(
      inputs$err_y -
        2 * rowSums(inputs$err_xy * repeated_coef) +
        rowSums(cov_times_coef * repeated_coef)
    )
  }
  vapply(seq_len(inputs$nobs), function(row_index) {
    inputs$err_y[[row_index]] -
      2 * sum(coef_predictors * inputs$err_xy[row_index, ]) +
      as.numeric(crossprod(
        coef_predictors,
        inputs$err_x[, , row_index] %*% coef_predictors
      ))
  }, numeric(1))
}

mmom_tilde_crossprod <- function(inputs, coef_predictors, case_weight) {
  if (inputs$backend == "rcpp") {
    # The native path replaces the row-wise variance loop.
    return(mmom_cpp_call(
      mmom_tilde_crossprod_cpp,
      inputs$err_x,
      inputs$err_xy,
      coef_predictors,
      case_weight
    ))
  }
  tilde <- matrix(0, inputs$nobs, inputs$ncoef)
  for (row_index in seq_len(inputs$nobs)) {
    tilde[row_index, inputs$predictor_index] <-
      inputs$err_xy[row_index, ] -
      inputs$err_x[, , row_index] %*% coef_predictors
  }
  crossprod(tilde * case_weight)
}

mmom_correction <- function(lambda, nobs, alpha) {
  if (is.finite(lambda) && lambda <= 1 + 1 / nobs) {
    lambda - 1 / nobs - alpha / nobs
  } else {
    1 - alpha / nobs
  }
}

mmom_is_failure <- function(value) {
  inherits(value, "S7_object") &&
    S7::prop_exists(value, "converged") &&
    !isTRUE(S7::prop(value, "converged"))
}

mmom_failure <- function(step, message, diagnostics = list()) {
  class_mmom_fit(
    converged = FALSE,
    status_code = 1L,
    message = paste(step, message, sep = ": "),
    coefficients = NULL,
    vcov = NULL,
    preliminary_coef = NULL,
    corrected_design_moment = NULL,
    analysis_weights = NULL,
    diagnostics = c(list(step = step, reason = message), diagnostics)
  )
}
