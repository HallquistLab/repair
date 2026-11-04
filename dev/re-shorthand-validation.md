# Concise random-coefficient references

November 4, 2026. `re(day)` now accepts an unquoted scalar coefficient name as
shorthand for `re("day")`. This is literal formula parsing, not evaluation of a
variable named `day`. No estimator or inference change accompanies the syntax.

## Identity and ambiguity

Quoted unqualified references already work when the catalog identifies one
coefficient. This includes a single named or unnamed MLM, and can include several
MLMs when only one has that coefficient name. A single MLM can still contain the
same coefficient in multiple grouping factors, requiring `group`. Multiple MLMs
with shared coefficient names require model qualification. Errors now identify
matching model/group/term entries and the relevant qualifier, rather than always
suggesting a model name even for within-model group ambiguity.

The parser converts only bare symbols to coefficient-name strings and then uses
the same alias handling, catalog lookup and canonical labels as quoted input.
Model/group qualifiers remain literal strings. Expressions, function calls,
numeric term arguments and the bare dot remain rejected. Backticked nonsyntactic
names are supported. A programmatically chosen name must be inserted into the
formula explicitly; `re(variable)` never looks up the value of `variable`.

## Evidence

`test-re-shorthand.R` verifies identical quoted/bare structural fits, predictions,
covariance and diagnostics; named and unqualified forms; nonsyntactic names and
intercept aliases; an active binding that must not be evaluated; rejected calls
and qualifier variables; model/group ambiguity including a single MLM; integrated
calls for both native engines and both computation backends; and reference replay
and paired bootstrap draws, estimates, alphas, covariance and intervals.
The existing grammar rejection test now uses `re(get("Days"))` instead of the
newly supported `re(Days)` form. Other existing fixtures and calculations remain
unchanged.

README and the main/nlme guides now show concise one-call fitting with unquoted
coefficient names and no unnecessary model qualifier. Installed help and the
workflow explain the literal-name and ambiguity contracts. The source-archive
hash, dated checks and hosted results are retained in `cran-comments.md` and the
ignored local `check/cran/` bundle. Earlier Win-builder receipts apply to their
own source archive; they do not validate a subsequent candidate. PRs 21 and 22
remain excluded, and no CRAN submission is made.

## Current archive evidence

The exact archive `repair_1.0.0.tar.gz` has SHA-256
`3cf98f3d57b5df3e64838b1183a7ef9523b01362ed799ea68b27c16450fcb171`.
Full `R CMD check --as-cran` reports zero errors, zero warnings and the two
NOTEs New submission and unavailable external current-time verification. All
4,076 assertions pass without failures, warnings or skips; examples, vignettes,
PDF and HTML manuals pass. The older dependency archive check passes the same
4,076 assertions with no errors/warnings and only the clock-verification NOTE.
The README runs against the checked installation and all 48 site pages validate.

One initial source-install load check segfaulted while source-loading tests and
installation were running against shared native build files. A clean isolated
installation passed temporary and final load checks; the crash did not recur.
Subsequent archive checks used separate unpacked build directories and passed.
The initial and isolated logs are retained for transparency; compiler diagnostics
were not suppressed and no runtime code change was used to bypass the load check.

The previous candidate and its Win-builder receipt are retained under
`check/cran-bb0ed01/`; the current source/check bundle is `check/cran/`. No new
Win-builder upload is implied by the earlier receipt. Check the exact final
archive before submission, once release preparation is frozen.
