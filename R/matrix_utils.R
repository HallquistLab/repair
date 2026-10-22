repair_solve <- function(
    lhs,
    rhs,
    method = c("auto", "qr", "chol"),
    tolerance = 1e-10,
    fail_on_error = FALSE
) {
  requested_method <- method[[1L]]
  method <- rlang::arg_match(method)

  if (!is.logical(fail_on_error) ||
      length(fail_on_error) != 1L ||
      is.na(fail_on_error)) {
    repair_error("`fail_on_error` must be a single TRUE or FALSE.")
  }

  lhs <- as.matrix(lhs)

  diagnostics <- list(
    requested_method = requested_method,
    solver = method,
    rank = NA_integer_,
    dimension = NA_integer_,
    condition_number = Inf
  )

  out <- class_mat_solve(
    success = FALSE,
    message = NULL,
    lhs = lhs,
    rhs = rhs,
    solution = NULL,
    inverse = NULL,
    solver = method,
    diagnostics = diagnostics
  )

  fail <- function(message) {
    if (isTRUE(fail_on_error)) {
      repair_error(message)
    }
    out@success <- FALSE
    out@message <- message
    out@solver <- out@diagnostics$solver
    out
  }

  if (length(dim(lhs)) != 2L || nrow(lhs) != ncol(lhs)) {
    return(fail("`lhs` must be a square matrix."))
  }

  if (any(!is.finite(lhs))) {
    return(fail("`lhs` must contain only finite values."))
  }

  if (is.numeric(rhs) && any(!is.finite(rhs))) {
    return(fail("`rhs` must contain only finite values."))
  }

  diagnostics$dimension <- ncol(lhs)
  diagnostics$condition_number <- tryCatch(
    kappa(lhs, exact = TRUE),
    error = function(error) Inf
  )

  symmetric <- isSymmetric(lhs, tol = tolerance)

  if (identical(method, "auto")) {
    method <- if (symmetric) "chol" else "qr"
    diagnostics$solver <- method
  }

  out@diagnostics <- diagnostics
  out@solver <- method

  solve_qr <- function() {
    qr_lhs <- qr(
      lhs,
      tol = tolerance,
      LAPACK = FALSE
    )
    solution <- qr.coef(qr_lhs, rhs)
    list(
      success = TRUE,
      rank = qr_lhs$rank,
      solution = solution
    )
  }

  solve_chol <- function() {
    chol_factor <- chol(lhs)
    solution <- backsolve(
      chol_factor,
      forwardsolve(t(chol_factor), rhs)
    )
    list(
      success = TRUE,
      solution = solution,
      inverse = chol2inv(chol_factor)
    )
  }

  if (identical(method, "chol")) {
    if (!symmetric) {
      if (identical(requested_method, "auto")) {
        method <- "qr"
        diagnostics$solver <- method
      } else {
        return(fail(
          "`lhs` must be symmetric for the Cholesky solver."
        ))
      }
    } else {
      chol_result <- tryCatch(
        solve_chol(),
        error = function(error) {
          list(
            success = FALSE,
            message = conditionMessage(error)
          )
        }
      )

      if (isTRUE(chol_result$success)) {
        diagnostics$rank <- ncol(lhs)

        if (any(!is.finite(chol_result$solution))) {
          return(fail(
            "The Cholesky solve returned a non-finite result."
          ))
        }

        out@success <- TRUE
        out@message <- "The Cholesky solve was successful."
        out@solution <- chol_result$solution
        out@inverse <- chol_result$inverse
        out@solver <- "chol"

        diagnostics$solver <- "chol"
        out@diagnostics <- diagnostics

        return(out)
      }

      if (identical(requested_method, "auto")) {
        method <- "qr"
        diagnostics$solver <- method
      } else {
        return(fail(
          paste(
            "The Cholesky solve failed:",
            chol_result$message
          )
        ))
      }
    }
  }

  out@solver <- method
  out@diagnostics <- diagnostics

  qr_result <- tryCatch(
    solve_qr(),
    error = function(error) {
      list(
        success = FALSE,
        message = conditionMessage(error)
      )
    }
  )

  if (!isTRUE(qr_result$success)) {
    return(fail(
      paste(
        "The QR solve failed:",
        qr_result$message
      )
    ))
  }

  diagnostics$rank <- qr_result$rank
  out@diagnostics <- diagnostics

  if (qr_result$rank < ncol(lhs)) {
    return(fail("The solver matrix was rank-deficient."))
  }

  if (any(!is.finite(qr_result$solution))) {
    return(fail(
      "The QR solve returned a non-finite result."
    ))
  }

  out@success <- TRUE
  out@message <- "The QR solve was successful."
  out@solution <- qr_result$solution
  out@inverse <- NULL
  out@solver <- "qr"
  out@diagnostics <- diagnostics

  out
}