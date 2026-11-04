# Load the frozen, pre-optimization estimator in its own lexical environment.
# It resolves unchanged shared covariance/solver helpers in the namespace, but
# its fit pipeline, moment steps, SEs, guards and exhaustive search are frozen.
step_up_legacy <- new.env(parent = asNamespace("repair"))
sys.source(test_path("fixtures", "step-up-60c46ec.R"), step_up_legacy)

expect_step_up_unchanged <- function(actual, expected) {
  expect_identical(actual@converged, expected@converged)
  expect_identical(actual@status_code, expected@status_code)
  expect_identical(actual@message, expected@message)
  for (field in c("coefficients", "vcov", "preliminary_coef", "corrected_design_moment",
                  "analysis_weights", "steps"))
    expect_identical(S7::prop(actual, field), S7::prop(expected, field), info = field)
  keep <- setdiff(names(expected@diagnostics), "step_up")
  expect_identical(actual@diagnostics[keep], expected@diagnostics[keep])
  expect_identical(names(actual@diagnostics$step_up), names(expected@diagnostics$step_up))
  for (stage in names(expected@diagnostics$step_up)) {
    old <- expected@diagnostics$step_up[[stage]]
    new <- actual@diagnostics$step_up[[stage]]
    keep <- setdiff(names(old), c("selected", "candidates"))
    expect_identical(new[keep], old[keep])
    rownames(new$selected) <- rownames(old$selected) <- NULL
    expect_identical(new$selected, old$selected)
    # Unneeded coarse candidates are omitted; every evaluated point has exactly
    # the same eligibility, guard measurements and reason as the old estimator.
    previous <- old$candidates[match(new$candidates$alpha, old$candidates$alpha), , drop = FALSE]
    rownames(previous) <- rownames(new$candidates) <- NULL
    expect_identical(new$candidates, previous)
  }
}
