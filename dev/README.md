# Development plans and audits

These documents describe design decisions, evidence, and future work. Historical
milestone audits retain their original scope; later follow-ups supersede earlier
support boundaries. README and executable vignettes describe current user-facing
behavior. A design draft does not make its proposed functions available.

## Current planning

- [First-release readiness audit](release-readiness-audit.md): accomplished work,
  dated CI evidence, merged versus pending website status, concrete release gaps,
  short-term user improvements, and a separate CRAN preparation path.
- [Person-level bootstrap and refitting design](bootstrap-refit-design.md): a
  future-version proposal covering resampling identity, replay/cohort contracts,
  cross-model covariance, failure handling, inference, and validation. Not a
  first-release commitment and not an implemented API.
- [Package design plan](package_plan.Rmd): implemented interfaces, statistical
  boundaries, and the development/acceptance workflow.

## Implementation evidence

- [Design audit](design-audit.md): original workflow corrections and follow-ups.
- [nlme validation](nlme-validation.md): effective rows, residual covariance,
  independent scoring and end-to-end reference checks.
- [Multiple-model validation](multimodel-validation.md): effect identities,
  cross-model covariance contracts, paired processes, and reproduction.
- [Reliability validation](reliability-validation.md): manuscript definitions,
  cohort/undefined-value behavior, and independent references.
- [Extractor contracts](extractor-contracts.md): public table schemas, ordering,
  keys, covariance reconstruction, and failure behavior.
- [CI checks](ci-checks.md): platform/R-version matrix, artifacts, local
  reproduction, and warning policy.
- [Fixture provenance](../tests/testthat/fixtures/README.md): independent
  numerical/reliability reference sources and deliberate regeneration policy.

The website work is tracked in [PR #8](https://github.com/zacharyvig/repair/pull/8).
Its build is verified but publication is still pending as of the release audit;
consult that audit's timestamp rather than treating this as a live status page.
