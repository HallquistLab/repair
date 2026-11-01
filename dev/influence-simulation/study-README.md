# REPAIR influence / bootstrap / stacked-sandwich comparison

Read [report.html](report.html) or its [source](report.Rmd). The immutable executed
protocol is [design.md](design.md); postprocessing/report files are separate.
All earlier studies are retained unchanged. This reuses the robust screen's
600 measurement datasets, with two paired structural designs and two truths,
rather than claiming independent confirmatory evidence.

## Reproduce

The executed numerical package is commit
558bdcead1c15cedd1e332aed6f01e86d49f3721, parent robust-refit head 7770874.
Its preserved installed library is `/private/tmp/repair-if-lib`, R 4.5.2 /
lme4 2.0.6. Do not overwrite that library during PR builds. Later scope-only
guards and documentation are checked separately against the recorded outputs.

From this directory:

```sh
REPAIR_STUDY_LIB=/private/tmp/repair-if-lib REPAIR_STUDY_WORKERS=1 Rscript run.R benchmark
REPAIR_STUDY_LIB=/private/tmp/repair-if-lib REPAIR_STUDY_WORKERS=12 Rscript run.R main
REPAIR_STUDY_LIB=/private/tmp/repair-if-lib REPAIR_STUDY_WORKERS=4 Rscript run.R sensitivity
Rscript summarize.R
Rscript validate.R
Rscript -e 'rmarkdown::render("report.Rmd")'
```

Requires the previous `../robust_refit` source archive and actual native-bootstrap
checkpoints, plus the research repository's `R/stacked_sandwich_helpers.R` and
`R/derivative_backends.R`. These GPL research helpers are sourced unchanged,
not copied into the MIT package. The driver snapshots library, sources,
comparator and original checkpoint hashes. Resuming refuses incompatible
hashes; do not mix a different package binary or source revision into a phase.
Cross-platform bitwise reproducibility is not promised.

`results/main` contains 600 synthetic checkpoints with influences, covariance
blocks, diagnostics, interval records, cached focal bootstrap arrays/draws and
legacy checks. `results/sensitivity` contains 18 B1999 checkpoints, exact B499
prefixes and 72 half-step derivative checks. Both include provenance/session
files. The six excluded timing/validation pilots are in results/benchmark;
an initial pilot adapter failure is separately retained at
`/private/tmp/repair-if-benchmark-initial` and was corrected before the main run.
No main/sensitivity driver errors occurred.

`results/summary` contains all 17,200 interval records, 172 performance cells,
paired common-eligibility comparisons, availability diagnostics, signed
covariance contributions, HC0 identity checks, public bootstrap checks,
normal-critical legacy diagnostics, B/derivative sensitivity, plots, manifests
and independent validation output. Coverage/zero-exclusion denominators are
available intervals; delivery measures also use all planned datasets. No missing
interval is counted as a valid non-rejection.

The compact PR bundle in repair's `dev/influence-simulation` stores the endpoint
records, summaries, plots, executed drivers and provenance plus a base-R
verifier. It permits review of rates without fitting thousands of models. The
full cached/native histories remain in this research archive. Arm-B focal
bootstrap reuse is justified only by this design's common-centering invariance;
the report does not reconstruct or claim intercept intervals.
