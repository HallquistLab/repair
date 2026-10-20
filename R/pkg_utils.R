# this file will contain misc utility functions for the repair package

repair_warn <- function(msg) {
  cli::cli_warn(msg)
}

repair_error <- function(msg) {
  cli::cli_abort(msg)
}