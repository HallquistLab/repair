repair_mlm_pkg <- function(mlm) {
  if (inherits(mlm, "lme")) {
    return("nlme")
  } else if (inherits(mlm, "lmerMod")) {
    return("lme4")
  } else if (inherits(mlm, "glmerMod")) {
    return("lme4")
  } else if (inherits(mlm, "brmsfit")) {
    return("brms")
  } else {
    repair_error("Unsupported mlm object class: ", class(mlm))
  }
}

repair_check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    repair_error("Package '", pkg, "' is required but not installed.")
  }
}

# get rows corresponding to each level of a grouping variable
repair_group_rows <- function(grouping, levels) {
  grouping <- as.character(grouping)
  rows <- split(
    seq_along(grouping),
    grouping,
    drop = TRUE
  )
  rows[levels]
}

# Functions for retrieving relevant elements from mixed-model fits.
repair_get_mlm_elements <- S7::new_generic(
  name = "repair_get_mlm_elements",
  dispatch_args = "model"
)

S7::method(repair_get_mlm_elements, class_merMod) <- function(model) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    repair_error(
      "The {.pkg lme4} package is required but not installed."
    )
  }

  # Fixed-effects model components
  model_frame <- stats::model.frame(model)
  fixed_design <- as.matrix(lme4::getME(model, "X"))
  fixed_coefs <- lme4::fixef(model)
  response <- as.numeric(
    stats::model.response(model_frame)
  )
  offset <- stats::model.offset(model_frame)
  if (is.null(offset)) {
    offset <- rep(0, length(response))
  }
  offset <- as.numeric(offset)
  observation_weights <- stats::model.weights(model_frame)
  if (is.null(observation_weights)) {
    observation_weights <- rep(1, length(response))
  }
  observation_weights <- as.numeric(observation_weights)
  if (any(!is.finite(observation_weights)) ||
    any(observation_weights <= 0)) {
    repair_error(
      "The fitted model has invalid prior weights."
    )
  }

  # Basic dimensions and checks
  nobs <- nrow(fixed_design)
  if (length(response) != nobs ||
    length(offset) != nobs ||
    length(observation_weights) != nobs) {
    repair_error(
      paste(
        "Extracted {.pkg lme4} model components have",
        "inconsistent row counts."
      )
    )
  }

  # Random-effects components
  ranef_object <- lme4::ranef(model, condVar = FALSE)
  random_terms <- lme4::getME(model, "cnms")
  grouping <- lme4::getME(model, "flist")
  group_names <- names(random_terms)
  if (!length(group_names)) {
    repair_error(
      "The fitted model has no random-effects grouping factors."
    )
  }
  if (!identical(
    group_names,
    names(grouping)
  )) {
    repair_error(
      "The random-effects terms and grouping factors have",
      "inconsistent names."
    )
  }
  group_levels <- setNames(
    lapply(
      group_names,
      function(group_name) {
        rownames(ranef_object[[group_name]])
      }
    ),
    group_names
  )
  group_rows <- setNames(
    lapply(
      group_names,
      function(group_name) {
        repair_group_rows(
          grouping = grouping[[group_name]],
          levels = group_levels[[group_name]]
        )
      }
    ),
    group_names
  )

  # Full random-effects design matrix
  random_design_full <- as.matrix(lme4::getME(model, "Z"))
  if (nrow(random_design_full) != nobs) {
    repair_error(
      "The random-effects design has an invalid number of rows."
    )
  }
  group_ptr <- lme4::getME(model, "Gp")
  if (length(group_ptr) != length(group_names) + 1L) {
    repair_error(
      "The random-effects column pointers have an invalid length."
    )
  }
  random_design <-
    random_covariance <-
    setNames(
      vector("list", length(group_names)),
      group_names
    )
  variance_components <- lme4::VarCorr(model)

  # Extract random-effects designs and covariance matrices
  for (group_index in seq_along(group_names)) {
    group_name <- group_names[[group_index]]
    terms <- random_terms[[group_name]]
    nterms <- length(terms)
    nlevels <- length(group_levels[[group_name]])
    if (!nterms || !nlevels) {
      repair_error(
        "Grouping factor {.val {group_name}} has no random-effect",
        "terms or levels."
      )
    }

    # Gp contains zero-based column offsets into Z.
    start <- group_ptr[[group_index]] + 1L
    end <- group_ptr[[group_index + 1L]]
    expected_columns <- nterms * nlevels
    if ((end - start + 1L) != expected_columns) {
      repair_error(
        paste(
          "The random-effects design for grouping factor",
          "{.val {group_name}} has an unexpected number of columns."
        )
      )
    }
    group_columns <- start:end
    group_design <- random_design_full[
      ,
      group_columns,
      drop = FALSE
    ]

    # lme4 orders columns by grouping level and random-effect term.
    design_names <- paste(
      rep(group_levels[[group_name]], each = nterms),
      rep(terms, times = nlevels),
      sep = "::"
    )
    colnames(group_design) <- design_names
    random_design[[group_name]] <- group_design

    # VarCorr() supplies the covariance matrix for one grouping level.
    group_covariance <- as.matrix(
      variance_components[[group_name]]
    )
    if (!identical(
      dim(group_covariance),
      c(nterms, nterms)
    )) {
      repair_error(
        paste(
          "The random-effects covariance for grouping factor",
          "{.val {group_name}} does not match its random-effect terms."
        )
      )
    }
    expanded_covariance <- kronecker(
      diag(nlevels),
      group_covariance
    )
    dimnames(expanded_covariance) <- list(
      design_names,
      design_names
    )
    random_covariance[[group_name]] <-
      expanded_covariance
  }

  # Residual covariance operator
  #
  # For a Gaussian lme4 model:
  #   R = diag(sigma^2 / weights)
  # Therefore:
  #   R^-1 A = A * weights / sigma^2
  #
  # Contributions from non-target random effects are added during
  # REPAIR scoring.
  sigma2 <- stats::sigma(model)^2
  if (!is.finite(sigma2) || sigma2 <= 0) {
    repair_error(
      "The fitted model must have a finite, positive residual variance."
    )
  }
  residual_covariance <- list(
    type = "independent",
    solve = function(group = NULL,
                     level = NULL,
                     rows,
                     value) {
      if (missing(rows)) {
        repair_error(
          "Rows must be supplied when solving with the residual",
          "covariance operator."
        )
      }
      value_matrix <- as.matrix(value)
      if (nrow(value_matrix) != length(rows)) {
        repair_error(
          "The value supplied to the residual covariance operator",
          "has incompatible dimensions."
        )
      }
      value * observation_weights[rows] / sigma2
    },
    matrix = function(group = NULL,
                      level = NULL,
                      rows) {
      if (missing(rows)) {
        repair_error(
          "Rows must be supplied when requesting a residual",
          "covariance matrix."
        )
      }
      diag(
        sigma2 / observation_weights[rows],
        nrow = length(rows),
        ncol = length(rows)
      )
    }
  )

  class_mlm_elements(
    backend = "lme4",
    model = model,
    ranef = ranef_object,
    fixed_design = fixed_design,
    fixed_coefs = fixed_coefs,
    response = response,
    offset = offset,
    observation_weights = observation_weights,
    random_design = random_design,
    random_terms = random_terms,
    grouping = grouping,
    group_levels = group_levels,
    group_rows = group_rows,
    group_ptr = group_ptr,
    random_covariance = random_covariance,
    residual_covariance = residual_covariance,
    diagnostics = list(
      nobs = nobs,
      nfixed = ncol(fixed_design),
      groups = group_names,
      residual_variance = sigma2,
      random_effects = stats::setNames(
        lapply(
          group_names,
          function(group_name) {
            list(
              terms = random_terms[[group_name]],
              levels = group_levels[[group_name]],
              nterms = length(random_terms[[group_name]]),
              nlevels = length(group_levels[[group_name]])
            )
          }
        ),
        group_names
      )
    )
  )
}


S7::method(repair_get_mlm_elements, class_lme) <- function(model) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    repair_error(
      "The {.pkg nlme} package is required but not installed."
    )
  }

  # Fixed-effects model components
  model_data <- nlme::getData(model)
  response <- as.numeric(nlme::getResponse(model))
  fixed_coefs <- nlme::fixef(model)
  fixed_design <- stats::model.matrix(
    formula(model, fixed = TRUE),
    data = model_data
  )

  # nlme::lme() does not use lme4-style offsets or prior observation
  # weights. Variance-function information is represented in the
  # residual covariance blocks below.
  offset <- numeric()
  observation_weights <- numeric()

  nobs <- nrow(fixed_design)
  if (length(response) != nobs) {
    repair_error(
      "The response and fixed-effects design have inconsistent row counts."
    )
  }

  # Random-effects grouping structure
  grouping_data <- model$groups
  if (is.null(grouping_data) ||
    !ncol(grouping_data)) {
    repair_error(
      "Could not identify the grouping structure of the {.pkg nlme} model."
    )
  }
  grouping_data <- as.data.frame(grouping_data)
  re_struct <- model$modelStruct$reStruct
  if (is.null(re_struct) ||
    !length(re_struct)) {
    repair_error(
      "The fitted {.pkg nlme} model has no random-effects structure."
    )
  }
  group_names <- names(re_struct)
  if (is.null(group_names) ||
    any(!nzchar(group_names))) {
    repair_error(
      "Could not identify the grouping factors in the random-effects",
      "structure."
    )
  }
  missing_groups <- setdiff(group_names, names(grouping_data))
  if (length(missing_groups)) {
    repair_error(
      "The following random-effects grouping factors were not found",
      "in the fitted model data: {.val {missing_groups}}"
    )
  }

  # Initialize normalized random-effects structures
  random_effects <-
    random_terms <-
    random_design <-
    random_covariance <-
    grouping <-
    group_levels <-
    group_rows <-
    setNames(
      vector("list", length(group_names)),
      group_names
    )

  # Extract each random-effects grouping factor
  for (group_index in seq_along(group_names)) {
    group_name <- group_names[[group_index]]
    grouping_factor <- as.character(grouping_data[[group_name]])
    group_levels_current <- unique(grouping_factor)
    group_levels_current <- group_levels_current[
      !is.na(group_levels_current)
    ]
    group_rows_current <- repair_group_rows(
      grouping = grouping_factor,
      levels = group_levels_current
    )
    if (any(lengths(group_rows_current) == 0L)) {
      repair_error(
        "One or more grouping levels have no observations."
      )
    }

    random_formula <- stats::formula(re_struct[[group_name]])
    random_design_all <- stats::model.matrix(
      random_formula,
      data = model_data
    )
    if (nrow(random_design_all) != nobs) {
      repair_error(
        "The reconstructed random-effects design for grouping factor",
        "{.val {group_name}} has an invalid number of rows."
      )
    }
    random_design_terms <- colnames(random_design_all)
    if (is.null(random_design_terms) ||
      !length(random_design_terms)) {
      repair_error(
        "Could not determine the random-effects terms for grouping factor",
        "{.val {group_name}}."
      )
    }

    # Obtain random-effect estimates for this grouping level.
    ranef_group <- tryCatch(
      {
        if (length(group_names) == 1L) {
          nlme::ranef(model)
        } else {
          nlme::ranef(
            model,
            level = group_index
          )
        }
      },
      error = function(error) {
        repair_error(
          "Could not extract random effects for grouping factor",
          "{.val {group_name}}: ",
          conditionMessage(error)
        )
      }
    )
    if (!is.data.frame(ranef_group)) {
      ranef_group <- as.data.frame(ranef_group)
    }
    fitted_terms <- colnames(ranef_group)
    if (is.null(fitted_terms) ||
      !length(fitted_terms)) {
      repair_error(
        "Could not determine the fitted random-effects terms for",
        "grouping factor {.val {group_name}}."
      )
    }
    if (ncol(random_design_all) != length(fitted_terms)) {
      repair_error(
        "The reconstructed random-effects design for grouping factor",
        "{.val {group_name}} has {ncol(random_design_all)} columns,",
        "but the fitted model has {length(fitted_terms)} random-effect",
        "terms."
      )
    }
    if (!identical(random_design_terms, fitted_terms)) {
      repair_error(
        "Could not unambiguously map the reconstructed random-effects",
        "design to the fitted terms for grouping factor",
        "{.val {group_name}}."
      )
    }

    # Use the fitted random-effects row names when available because
    # they define the order used by nlme::ranef().
    fitted_levels <- rownames(ranef_group)

    if (!is.null(fitted_levels)) {
      missing_levels <- setdiff(
        fitted_levels,
        group_levels_current
      )
      if (length(missing_levels)) {
        repair_error(
          "The fitted random-effects output contains levels not found",
          "in the grouping data for {.val {group_name}}."
        )
      }
      group_levels_current <- fitted_levels
      group_rows_current <- repair_group_rows(
        grouping = grouping_factor,
        levels = group_levels_current
      )
    }

    random_effects[[group_name]] <- ranef_group
    random_terms[[group_name]] <- fitted_terms
    random_design[[group_name]] <- random_design_all
    grouping[[group_name]] <- grouping_factor
    group_levels[[group_name]] <- group_levels_current
    group_rows[[group_name]] <- group_rows_current

    # Extract the random-effects covariance for one grouping level
    # and expand it across grouping levels.
    random_covariance0 <- tryCatch(
      nlme::getVarCov(
        model,
        type = "random.effects"
      ),
      error = function(error) {
        repair_error(
          "Could not extract the random-effects covariance for grouping",
          "factor {.val {group_name}}: ",
          conditionMessage(error)
        )
      }
    )

    # Depending on the nlme model, getVarCov() may return a matrix or
    # a named list. Extract the relevant matrix when necessary.
    if (is.list(random_covariance0) && !is.matrix(random_covariance0)) {
      if (group_name %in% names(random_covariance0)) {
        random_covariance0 <- random_covariance0[[group_name]]
      } else if (length(random_covariance0) == 1L) {
        random_covariance0 <- random_covariance0[[1L]]
      } else {
        repair_error(
          "Could not identify the random-effects covariance matrix for",
          "grouping factor {.val {group_name}}."
        )
      }
    }
    random_covariance0 <- as.matrix(random_covariance0)
    if (!identical(
      dim(random_covariance0),
      c(length(fitted_terms), length(fitted_terms))
    )) {
      repair_error(
        "The random-effects covariance for grouping factor",
        "{.val {group_name}} has dimensions",
        "{.val {paste(dim(random_covariance0), collapse = ' x ')}};",
        "expected {length(fitted_terms)} by {length(fitted_terms)}."
      )
    }
    nlevels <- length(group_levels_current)
    random_covariance_expanded <- kronecker(diag(nlevels), random_covariance0)
    random_design_names <- paste(
      rep(group_levels_current, each = length(fitted_terms)),
      rep(fitted_terms, times = nlevels),
      sep = "::"
    )
    dimnames(random_covariance_expanded) <- list(
      random_design_names,
      random_design_names
    )
    random_covariance[[group_name]] <-
      random_covariance_expanded
  }

  # Conditional residual covariance
  #
  # For each level of the residual grouping factor,
  # getVarCov(..., type = "conditional") supplies the fitted
  # ordinary residual covariance. Contributions from non-target
  # random effects are added later during REPAIR scoring.
  residual_group <- group_names[[length(group_names)]]
  residual_levels <- group_levels[[residual_group]]
  residual_rows <- group_rows[[residual_group]]
  residual_blocks <- setNames(
    vector("list", length(residual_levels)),
    residual_levels
  )
  for (level in residual_levels) {
    covariance <- tryCatch(
      nlme::getVarCov(
        model,
        individuals = level,
        type = "conditional"
      ),
      error = function(error) {
        repair_error(
          "Could not extract the conditional residual covariance",
          "for grouping level {.val {level}}: ",
          conditionMessage(error)
        )
      }
    )
    # getVarCov() may return a named list of covariance objects.
    if (is.list(covariance) &&
      !is.matrix(covariance)) {
      if (level %in% names(covariance)) {
        covariance <- covariance[[level]]
      } else if (length(covariance) == 1L) {
        covariance <- covariance[[1L]]
      } else {
        repair_error(
          "Could not identify the conditional residual covariance for",
          "grouping level {.val {level}}."
        )
      }
    }
    covariance <- as.matrix(covariance)
    expected_dimension <- length(residual_rows[[level]])
    if (!identical(
      dim(covariance),
      c(expected_dimension, expected_dimension)
    )) {
      repair_error(
        "The conditional residual covariance for grouping level",
        "{.val {level}} has dimensions",
        "{.val {paste(dim(covariance), collapse = ' x ')}};",
        "expected {expected_dimension} by {expected_dimension}."
      )
    }
    if (any(!is.finite(covariance))) {
      repair_error(
        "The conditional residual covariance for grouping level",
        "{.val {level}} contains non-finite values."
      )
    }
    residual_blocks[[level]] <- covariance
  }

  residual_operator <- list(
    type = "nlme_conditional",
    group = residual_group,
    blocks = residual_blocks,
    solve = function(group = NULL,
                     level,
                     rows = NULL,
                     value) {
      covariance <- residual_blocks[[level]]
      if (is.null(covariance)) {
        repair_error(
          "No residual covariance was found for grouping level",
          "{.val {level}}."
        )
      }
      if (!is.null(rows) &&
        length(rows) != nrow(covariance)) {
        repair_error(
          "The supplied rows do not match grouping level",
          "{.val {level}}."
        )
      }
      result <- repair_solve(
        lhs = covariance,
        rhs = value,
        method = "chol",
        fail_on_error = TRUE
      )
      result@solution
    },
    matrix = function(group = NULL,
                      level,
                      rows = NULL) {
      covariance <- residual_blocks[[level]]
      if (is.null(covariance)) {
        repair_error(
          "No residual covariance was found for grouping level",
          "{.val {level}}."
        )
      }
      if (!is.null(rows) &&
        length(rows) != nrow(covariance)) {
        repair_error(
          "The supplied rows do not match grouping level",
          "{.val {level}}."
        )
      }
      covariance
    }
  )

  # Return normalized model elements
  class_mlm_elements(
    backend = "nlme",
    model = model,
    ranef = random_effects,
    fixed_design = fixed_design,
    fixed_coefs = fixed_coefs,
    response = response,
    offset = offset,
    observation_weights = observation_weights,
    random_design = random_design,
    random_terms = random_terms,
    grouping = grouping,
    group_levels = group_levels,
    group_rows = group_rows,
    group_ptr = NULL,
    random_covariance = random_covariance,
    residual_covariance = residual_operator,
    diagnostics = list(
      nobs = nobs,
      nfixed = ncol(fixed_design),
      groups = group_names,
      residual_group = residual_group,
      residual_structure = list(
        correlation = model$modelStruct$corStruct,
        variance = model$modelStruct$varStruct
      )
    )
  )
}
