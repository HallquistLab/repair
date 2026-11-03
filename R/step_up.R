# Sequential alpha selection for the manuscript M-MOMS estimator. All estimates
# and candidate availability decisions use the canonical M-MOM pipeline; this
# file only implements the search and its predictor-block conditioning guards.
mmom_step_up_control <- function(control) {
  defaults <- list(coarse_grid_size = 9L, search_tolerance = 0.005,
    max_refinements = 12L, min_eigenvalue = sqrt(.Machine$double.eps),
    min_relative_eigenvalue = 0.05, max_condition_number = 1e5,
    candidate_alphas = NULL)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control)) ||
      any(!names(control) %in% names(defaults)))))
    repair_error("`step_up_control` must be a named list of documented search settings.")
  settings <- utils::modifyList(defaults, control, keep.null = TRUE)
  for (name in setdiff(names(defaults), "candidate_alphas")) {
    value <- settings[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0)
      repair_error("Step-up setting {.val {name}} must be a finite nonnegative scalar.")
  }
  if (settings$coarse_grid_size < 3 || settings$coarse_grid_size %% 1 != 0 ||
      settings$coarse_grid_size > .Machine$integer.max || settings$max_refinements %% 1 != 0 ||
      settings$max_refinements > .Machine$integer.max)
    repair_error("Step-up grid size must be an integer >= 3 and max_refinements a nonnegative integer.")
  if (settings$search_tolerance <= 0 || settings$max_condition_number < 1)
    repair_error("Step-up search_tolerance must be positive and max_condition_number >= 1.")
  if (!is.null(settings$candidate_alphas) && (!is.numeric(settings$candidate_alphas) ||
      !length(settings$candidate_alphas) || any(!is.finite(settings$candidate_alphas))))
    repair_error("Step-up candidate_alphas must be NULL or a nonempty finite numeric vector.")
  settings
}

mmom_step_up_guard <- function(fit, inputs, stage, alpha, settings) {
  result <- data.frame(alpha = alpha, usable = isTRUE(fit@converged), guard_pass = FALSE,
    reason = if (isTRUE(fit@converged)) "" else fit@message,
    min_eigenvalue = NA_real_, max_eigenvalue = NA_real_,
    relative_min_eigenvalue = NA_real_, condition_number = NA_real_)
  if (!isTRUE(fit@converged)) return(result)
  if (stage == "preliminary") {
    observed <- crossprod(inputs$design) / inputs$nobs
    corrected <- (crossprod(inputs$design) -
      fit@steps$step1@correction_prelim * fit@steps$step1@err_cov_design) / inputs$nobs
  } else {
    observed <- crossprod(fit@steps$step3@weighted_design)
    corrected <- fit@corrected_design_moment
  }
  index <- inputs$predictor_index
  info <- mmom_matrix_diagnostics(corrected[index, index, drop = FALSE])
  observed_info <- mmom_matrix_diagnostics(observed[index, index, drop = FALSE])
  relative <- if (is.finite(observed_info$max_eigen) && observed_info$max_eigen > sqrt(.Machine$double.eps))
    info$min_eigen / observed_info$max_eigen else NA_real_
  result$min_eigenvalue <- info$min_eigen
  result$max_eigenvalue <- info$max_eigen
  result$relative_min_eigenvalue <- relative
  result$condition_number <- info$condition_number
  checks <- c(eigenvalue_floor = is.finite(info$min_eigen) && info$min_eigen >= settings$min_eigenvalue,
    relative_eigenvalue_floor = is.finite(relative) && relative >= settings$min_relative_eigenvalue,
    condition_cap = is.finite(info$condition_number) && info$condition_number <= settings$max_condition_number)
  result$guard_pass <- all(checks)
  result$reason <- if (all(checks)) "ok" else names(checks)[which(!checks)[1L]]
  result
}

mmom_step_up <- function(inputs) {
  requested <- list(preliminary = inputs$preliminary_moment_control,
                    final = inputs$final_moment_control)
  active <- vapply(requested, function(x) x$correction == "step_up", logical(1))
  working <- inputs
  for (stage in names(requested)[active]) {
    field <- paste0(stage, "_moment_control")
    working[[field]] <- list(correction = "custom", alpha = requested[[stage]]$alpha)
  }
  fit <- mmom_fit_prepared(working)
  searches <- list()
  for (stage in names(requested)[active]) {
    field <- paste0(stage, "_moment_control")
    settings <- requested[[stage]]$step_up
    lower <- requested[[stage]]$alpha
    lambda <- fit@diagnostics[[if (stage == "preliminary") "lambda_prelim" else "lambda_final"]]
    upper <- if (length(lambda) == 1L && is.finite(lambda)) {
      if (lambda <= 1 + 1 / inputs$nobs) lambda * inputs$nobs - 1 else inputs$nobs
    } else lower
    # Cache candidates to avoid refitting the baseline and selected midpoint.
    fits <- list(fit)
    history <- mmom_step_up_guard(fit, working, stage, lower, settings)
    evaluate <- function(alpha) {
      known <- match(alpha, history$alpha)
      if (!is.na(known)) return(known)
      candidate <- working
      candidate[[field]]$alpha <- alpha
      candidate_fit <- mmom_fit_prepared(candidate)
      fits[[length(fits) + 1L]] <<- candidate_fit
      history <<- rbind(history, mmom_step_up_guard(candidate_fit, candidate, stage, alpha, settings))
      nrow(history)
    }
    if (upper > lower) {
      grid <- settings$candidate_alphas
      if (is.null(grid)) grid <- seq(lower, upper, length.out = settings$coarse_grid_size)
      grid <- sort(unique(c(lower, pmin(upper, pmax(lower, grid)))))
      for (alpha in grid) evaluate(alpha)
      passing <- history$alpha[history$guard_pass]
      if (length(passing)) {
        hi <- min(passing)
        below <- history$alpha[!history$guard_pass & history$alpha < hi]
        if (length(below)) {
          lo <- max(below)
          for (iteration in seq_len(settings$max_refinements)) {
            if (hi - lo <= settings$search_tolerance) break
            mid <- (hi + lo) / 2
            row <- evaluate(mid)
            if (history$guard_pass[row]) hi <- mid else lo <- mid
          }
        }
      }
    }
    passing <- which(history$guard_pass)
    selected <- if (length(passing)) passing[which.min(history$alpha[passing])] else 1L
    fit <- fits[[selected]]
    working[[field]]$alpha <- history$alpha[selected]
    searches[[stage]] <- list(alpha = history$alpha[selected], lower = lower, upper = upper,
      met_criteria = isTRUE(history$guard_pass[selected]), fallback = !length(passing),
      reason = if (length(passing)) "ok" else if (upper <= lower) "no_search_interval" else "no_passing_candidate",
      selected = history[selected, , drop = FALSE],
      candidates = history[order(history$alpha), , drop = FALSE], settings = settings)
  }
  # Preserve the requested policy for replay. Selected alphas are outcomes of the
  # search, not frozen controls to reuse for subsequent bootstrap samples.
  fit@diagnostics$preliminary_moment_control <- requested$preliminary
  fit@diagnostics$final_moment_control <- requested$final
  fit@diagnostics$alpha_prelim <- working$preliminary_moment_control$alpha
  fit@diagnostics$alpha_final <- working$final_moment_control$alpha
  fit@diagnostics$step <- fit@diagnostics[["step", exact = TRUE]] %||% "complete"
  fit@diagnostics$step_up <- searches
  fit
}
