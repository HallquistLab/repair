# this file will contain misc utility functions for the repair package

set_caller_context <- function() {
  rlang::local_error_call(rlang::caller_env())
}

repair_warn <- function(msg, .class = "repair_warning", .call) {
  if (missing(.call)) {
    cli::cli_warn(msg, class = .class)
  } else {
    cli::cli_warn(msg, class = .class, call = .call)
  }
}

repair_error <- function(msg, .class = "repair_error", .call) {
  if (missing(.call)) {
    cli::cli_abort(msg, class = .class)
  } else {
    cli::cli_abort(msg, class = .class, call = .call)
  }
}

repair_message <- function(msg, .class = "repair_message", .call) {
  if (missing(.call)) {
    cli::cli_inform(msg, class = .class)
  } else {
    cli::cli_inform(msg, class = .class, call = .call)
  }
}