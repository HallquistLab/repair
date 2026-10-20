# Small, shared operations used by the three M-MOM estimation steps.

`%||%` <- function(value, fallback) if (is.null(value)) fallback else value

# Backend and correction-control normalization.
mmom_backend <- function(backend) {
  if (length(backend) != 1L || !backend %in% c("auto", "r", "rcpp")) {
    repair_error("`backend` must be one of `auto`, `r`, or `rcpp`.")
  }
  if (backend == "auto") {
    return(if (exists("mmom_normalize_err_cov_cpp")) "rcpp" else "r")
  }
  backend
}

mmom_control <- function(control, npredictors, label) {
  if (!is.list(control)) {
    repair_error(sprintf("`%s` must be a list.", label))
  }
  alpha <- control$alpha %||% (npredictors + 2L)
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0) {
    repair_error(sprintf(
      "`%s$alpha` must be a finite, positive scalar.", label
    ))
  }
  list(
    regularize = isTRUE(control$regularize),
    alpha = max(alpha, npredictors + 2L)
  )
}

mmom_prepare <- function(
  data_list, preliminary_moment, variance_bread, backend
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
  err_cov <- mmom_normalize_err_cov(data$err_cov, nobs, npredictors, backend)
  err_var <- data$err_var
  if (is.null(err_var)) err_var <- rep(0, nobs)
  if (!is.numeric(err_var) || length(err_var) != nobs ||
      any(!is.finite(err_var)) || any(err_var < 0)) {
    repair_error(
      "`err_var` must be a finite, nonnegative vector of length nobs."
    )
  }
  err_cross <- mmom_normalize_err_cross(data$err_cross, nobs, npredictors)

  # Validate the predictor and joint outcome/predictor error blocks.
  mmom_validate_joint_cov(err_cov, err_var, err_cross, backend)
  design <- cbind(`(Intercept)` = 1, predictors)
  list(
    outcome = outcome,
    predictors = predictors,
    err_cov = err_cov,
    err_var = err_var,
    err_cross = err_cross,
    design = design,
    design_index = seq_len(npredictors + 1L) + 1L,
    predictor_index = seq.int(2L, npredictors + 1L),
    nobs = nobs,
    npredictors = npredictors,
    ncoef = npredictors + 1L,
    backend = backend,
    preliminary = mmom_control(
      preliminary_moment, npredictors, "preliminary_moment"
    ),
    variance = mmom_control(
      variance_bread, npredictors, "variance_bread"
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
  if (!is.null(data_list$y) && !is.null(data_list$x)) {
    data_list$outcome <- data_list$y
    data_list$predictors <- data_list$x
    data_list$err_var <- data_list$err_var %||%
      data_list$err_var_y %||% data_list$omega_y
    data_list$err_cross <- data_list$err_cross %||%
      data_list$err_cov_xy %||% data_list$omega_xy
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
    err_cov = lapply(
      data_list,
      function(observation) get_field(observation, c("err_cov", "omega_x"))
    ),
    err_var = vapply(
      data_list,
      function(observation) {
        as.numeric(get_field(
          observation, c("err_var", "err_var_y", "omega_y"), 0
        ))
      },
      numeric(1)
    ),
    err_cross = {
      values <- lapply(
        data_list,
        function(observation) {
          get_field(observation, c("err_cross", "err_cov_xy", "omega_xy"))
        }
      )
      if (all(vapply(values, is.null, logical(1)))) NULL else values
    }
  )
}

# Covariance coercion and positive-semidefinite checks.
mmom_normalize_err_cov <- function(
  err_cov, nobs, npredictors, backend, tolerance = 1e-8
) {
  out <- mmom_coerce_err_cov(err_cov, nobs, npredictors)
  if (backend == "rcpp") {
    normalized <- mmom_cpp_call(
      mmom_normalize_err_cov_cpp, out, nobs, npredictors, tolerance
    )
    if (normalized$bad_row > 0L) {
      repair_error(sprintf(
        "`err_cov` is not positive semidefinite for row %d.",
        normalized$bad_row
      ))
    }
    return(normalized$err_cov)
  }
  for (row_index in seq_len(nobs)) {
    out[, , row_index] <- mmom_symmetrize(out[, , row_index])
    diagnostics <- mmom_matrix_diagnostics(out[, , row_index])
    if (!is.finite(diagnostics$min_eigen) ||
        diagnostics$min_eigen < -tolerance * max(
          1, abs(diagnostics$max_eigen)
        )) {
      repair_error(sprintf(
        "`err_cov` is not positive semidefinite for row %d.", row_index
      ))
    }
  }
  out
}

mmom_coerce_err_cov <- function(err_cov, nobs, npredictors) {
  if (is.matrix(err_cov)) {
    if (!identical(dim(err_cov), c(npredictors, npredictors))) {
      repair_error(
        "A common `err_cov` matrix must be square with one row per predictor."
      )
    }
    out <- array(
      rep(as.numeric(err_cov), nobs),
      c(npredictors, npredictors, nobs)
    )
  } else if (is.array(err_cov) && length(dim(err_cov)) == 3L) {
    dimensions <- dim(err_cov)
    if (identical(dimensions, c(npredictors, npredictors, nobs))) {
      out <- err_cov
    } else if (identical(dimensions, c(nobs, npredictors, npredictors))) {
      out <- aperm(err_cov, c(2, 3, 1))
    } else {
      repair_error(
        paste(
          "`err_cov` must be npredictors by npredictors by nobs or",
          "nobs by npredictors by npredictors."
        )
      )
    }
  } else if (is.list(err_cov)) {
    if (length(err_cov) != nobs) {
      repair_error("An `err_cov` list must have length nobs.")
    }
    out <- array(NA_real_, c(npredictors, npredictors, nobs))
    for (row_index in seq_len(nobs)) {
      block <- as.matrix(err_cov[[row_index]])
      if (!identical(dim(block), c(npredictors, npredictors))) {
        repair_error(
          paste(
            "Every `err_cov` block must be square with",
            "one row per predictor."
          )
        )
      }
      out[, , row_index] <- block
    }
  } else {
    repair_error(
      paste(
        "`err_cov` must be a matrix, three-dimensional array,",
        "or list of matrices."
      )
    )
  }
  if (any(!is.finite(out))) {
    repair_error("`err_cov` contains non-finite values.")
  }
  out
}

mmom_normalize_err_cross <- function(err_cross, nobs, npredictors) {
  if (is.null(err_cross)) {
    return(matrix(0, nobs, npredictors))
  }
  if (is.vector(err_cross) && length(err_cross) == npredictors) {
    return(matrix(rep(err_cross, each = nobs), nobs, npredictors))
  }
  err_cross <- as.matrix(err_cross)
  if (!identical(dim(err_cross), c(nobs, npredictors))) {
    repair_error(
      paste(
        "`err_cross` must be length npredictors or an",
        "nobs by npredictors matrix."
      )
    )
  }
  if (any(!is.finite(err_cross))) {
    repair_error("`err_cross` contains non-finite values.")
  }
  err_cross
}

mmom_validate_joint_cov <- function(
  err_cov, err_var, err_cross, backend, tolerance = 1e-8
) {
  if (backend == "rcpp") {
    bad_row <- mmom_cpp_call(
      mmom_joint_err_cov_bad_row_cpp,
      err_cov,
      err_var,
      err_cross,
      tolerance
    )
    if (bad_row > 0L) {
      repair_error(sprintf(
        "The joint measurement-error covariance is not PSD for row %d.",
        bad_row
      ))
    }
    return(invisible(NULL))
  }
  for (row_index in seq_along(err_var)) {
    joint <- rbind(
      c(err_var[[row_index]], err_cross[row_index, ]),
      cbind(err_cross[row_index, ], err_cov[, , row_index])
    )
    diagnostics <- mmom_matrix_diagnostics(joint)
    if (!is.finite(diagnostics$min_eigen) ||
        diagnostics$min_eigen < -tolerance * max(
          1, abs(diagnostics$max_eigen)
        )) {
      repair_error(sprintf(
        "The joint measurement-error covariance is not PSD for row %d.",
        row_index
      ))
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

mmom_weighted_err_cov_sum <- function(err_cov, weights, backend) {
  if (backend == "rcpp") {
    return(mmom_symmetrize(mmom_cpp_call(
      mmom_weighted_err_cov_sum_cpp, err_cov, weights
    )))
  }
  result <- matrix(0, dim(err_cov)[1L], dim(err_cov)[2L])
  for (row_index in seq_len(dim(err_cov)[3L])) {
    result <- result + weights[[row_index]] * err_cov[, , row_index]
  }
  mmom_symmetrize(result)
}

# Moment construction, residual variance, and sandwich components.
mmom_error_moment <- function(err_cov_sum, err_var_sum, err_cross_sum) {
  nvars <- nrow(err_cov_sum) + 1L
  moment <- matrix(0, nvars + 1L, nvars + 1L)
  moment[1, 1] <- err_var_sum
  moment[1, 2:(nvars + 1L)] <- err_cross_sum
  moment[2:(nvars + 1L), 1] <- err_cross_sum

  # The intercept is measured without error.
  moment[2:(nvars + 1L), 2:(nvars + 1L)] <- {
    predictor_block <- matrix(0, nvars, nvars)
    predictor_block[-1L, -1L] <- err_cov_sum
    predictor_block
  }
  moment
}

mmom_composite_err_var <- function(inputs, coef_predictors) {
  if (inputs$backend == "rcpp") {
    # The native path replaces the row-wise covariance loop.
    cov_times_coef <- mmom_cpp_call(
      mmom_err_cov_matvec_cpp, inputs$err_cov, coef_predictors
    )
    repeated_coef <- rep(coef_predictors, each = inputs$nobs)
    return(
      inputs$err_var -
        2 * rowSums(inputs$err_cross * repeated_coef) +
        rowSums(cov_times_coef * repeated_coef)
    )
  }
  vapply(seq_len(inputs$nobs), function(row_index) {
    inputs$err_var[[row_index]] -
      2 * sum(coef_predictors * inputs$err_cross[row_index, ]) +
      as.numeric(crossprod(
        coef_predictors,
        inputs$err_cov[, , row_index] %*% coef_predictors
      ))
  }, numeric(1))
}

mmom_tilde_crossprod <- function(inputs, coef_predictors, case_weight) {
  if (inputs$backend == "rcpp") {
    # The native path replaces the row-wise variance loop.
    return(mmom_cpp_call(
      mmom_tilde_crossprod_cpp,
      inputs$err_cov,
      inputs$err_cross,
      coef_predictors,
      case_weight
    ))
  }
  tilde <- matrix(0, inputs$nobs, inputs$ncoef)
  for (row_index in seq_len(inputs$nobs)) {
    tilde[row_index, inputs$predictor_index] <-
      inputs$err_cross[row_index, ] -
      inputs$err_cov[, , row_index] %*% coef_predictors
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

mmom_solve <- function(left, right, step, failure_message, tolerance = 1e-10) {
  left <- as.matrix(left)
  if (nrow(left) != ncol(left)) {
    return(mmom_failure(
      step,
      "The solver matrix was not square.",
      list(solver = "qr", dimensions = dim(left))
    ))
  }

  qr_left <- qr(left, tol = tolerance, LAPACK = FALSE)
  rank <- qr_left$rank
  condition_number <- tryCatch(
    kappa(left, exact = TRUE),
    error = function(e) Inf
  )
  if (rank < ncol(left)) {
    return(mmom_failure(
      step,
      failure_message,
      list(
        solver = "qr",
        rank = rank,
        dimension = ncol(left),
        condition_number = condition_number
      )
    ))
  }

  result <- tryCatch(
    qr.coef(qr_left, right),
    error = function(e) NULL
  )
  if (is.null(result) || any(!is.finite(result))) {
    return(mmom_failure(
      step,
      "The rank-aware solve returned a non-finite result.",
      list(
        solver = "qr",
        rank = rank,
        dimension = ncol(left),
        condition_number = condition_number
      )
    ))
  }
  attr(result, "mmom_solver_diagnostics") <- list(
    solver = "qr",
    rank = rank,
    dimension = ncol(left),
    condition_number = condition_number
  )
  result
}

mmom_chol_solve <- function(matrix_value, right) {
  chol_factor <- chol(matrix_value)
  solution <- backsolve(
    chol_factor,
    forwardsolve(t(chol_factor), right)
  )
  list(
    solution = solution,
    inverse = chol2inv(chol_factor),
    diagnostics = list(
      solver = "chol",
      rank = ncol(matrix_value),
      dimension = ncol(matrix_value),
      condition_number = kappa(matrix_value, exact = TRUE)
    )
  )
}

mmom_is_failure <- function(value) {
  inherits(value, "S7_object") &&
    S7::prop_exists(value, "converged") &&
    !isTRUE(S7::prop(value, "converged"))
}

mmom_failure <- function(step, message, diagnostics = list()) {
  mmom_fit(
    converged = FALSE,
    status_code = 1L,
    message = paste(step, message, sep = ": "),
    coefficients = NULL,
    vcov = NULL,
    preliminary_coefficients = NULL,
    corrected_design_moment = NULL,
    analysis_weights = NULL,
    diagnostics = c(list(step = step, reason = message), diagnostics)
  )
}
