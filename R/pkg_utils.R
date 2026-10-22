# this file will contain misc utility functions for the repair package

set_caller_context <- function() {
  rlang::local_error_call(rlang::caller_env())
}

repair_warn <- function(msg, .class = "repair_warning", .call, .envir = rlang::caller_env()) {
  if (missing(.call)) {
    cli::cli_warn(msg, class = .class, .envir = .envir)
  } else {
    cli::cli_warn(msg, class = .class, call = .call, .envir = .envir)
  }
}

repair_error <- function(msg, .class = "repair_error", .call, .envir = rlang::caller_env()) {
  if (missing(.call)) {
    cli::cli_abort(msg, class = .class, .envir = .envir)
  } else {
    cli::cli_abort(msg, class = .class, call = .call, .envir = .envir)
  }
}

repair_message <- function(msg, .class = "repair_message", .call, .envir = rlang::caller_env()) {
  if (missing(.call)) {
    cli::cli_inform(msg, class = .class, .envir = .envir)
  } else {
    cli::cli_inform(msg, class = .class, call = .call, .envir = .envir)
  }
}

# Fall back to value if first one is NULL
`%||%` <- function(value, fallback) if (is.null(value)) fallback else value