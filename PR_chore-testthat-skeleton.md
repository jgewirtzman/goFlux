# chore: testthat skeleton

**Branch:** `chore/testthat-skeleton`
**Files:** `tests/testthat.R`, `tests/testthat/test-smoke.R` (new), `DESCRIPTION` (`testthat` in Suggests)

goFlux has no `tests/` directory. This adds the minimal testthat skeleton and one
smoke test (`goFlux()` + `best.flux()` on `manID.UGGA`), so that the feature
branches `feat/empirical-precision-mdf`, `feat/qc-flags`,
`feat/click-peak2-stacked` and `feat/aqua-diffusive-deebulliated` (all based on
this branch) only add their own test files. `R CMD check` stays as on master
apart from a new "checking tests ... OK".
