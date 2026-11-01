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
  block_terms <- lme4::getME(model, "cnms")
  grouping <- lme4::getME(model, "flist")
  group_names <- names(grouping)
  # lme4 can represent one grouping factor with several independent blocks
  # (for example, (Days || Subject)). Merge these in level-major order.
  random_terms <- stats::setNames(lapply(group_names, function(group) {
    unlist(block_terms[names(block_terms) == group], use.names = FALSE)
  }), group_names)
  if (!length(group_names)) {
    repair_error(
      "The fitted model has no random-effects grouping factors."
    )
  }
  if (!setequal(names(block_terms), group_names)) {
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
  if (length(group_ptr) != length(block_terms) + 1L) {
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

    if (anyDuplicated(terms)) {
      repair_error("Random coefficient names must be unique within each grouping factor for likelihood-only prediction extraction.")
    }
    blocks <- which(names(block_terms) == group_name)
    # Gp contains zero-based offsets for blocks, not distinct grouping factors.
    for (block in blocks) {
      if (group_ptr[block + 1L] - group_ptr[block] !=
          length(block_terms[[block]]) * nlevels) {
        repair_error("A random-effects block has incompatible dimensions.")
      }
    }
    group_columns <- unlist(lapply(seq_len(nlevels), function(level) {
      unlist(lapply(blocks, function(block) {
        width <- length(block_terms[[block]])
        group_ptr[block] + (level - 1L) * width + seq_len(width)
      }), use.names = FALSE)
    }), use.names = FALSE)
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

    # Independent blocks have zero cross-block covariance in the fitted model.
    group_covariance <- matrix(0, nterms, nterms)
    term_offset <- 0L
    for (block in blocks) {
      width <- length(block_terms[[block]])
      indices <- term_offset + seq_len(width)
      covariance <- as.matrix(variance_components[[block]])
      if (!identical(dim(covariance), c(width, width))) {
        repair_error("A random-effects covariance block has incompatible dimensions.")
      }
      group_covariance[indices, indices] <- covariance
      term_offset <- term_offset + width
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
  # REPAIR prediction extraction.
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
