# feat: empirical.prec() and a confidence level for MDF

**Branch:** `feat/empirical-precision-mdf`
**Files:** `R/empirical.prec.R` (new), `R/MDF.R`, `R/goFlux.R`, `R/best.flux.R`, `man/`, `NAMESPACE`, `tests/testthat/test-empirical.prec.R` (new; adds `testthat` to Suggests)

Datasheet precision underestimates the noise of field-worn analyzers (2.5x on
a UGGA in our data), so the MDF built on it is too optimistic. This adds

- `empirical.prec(dataframe, gastype, method = c("mad", "allan"), by = NULL)`:
  `"mad"` = MAD(dx) x 1.4826 / sqrt(2) over the whole record or per `by` group
  (robust to moves, breaths and chamber changes, so no quiet-period selection
  is needed); `"allan"` = SD(dx)/sqrt(2) per `UniqueID` on the flagged rows.
  The docs say plainly that feeding the result through the existing `prec`
  argument at import (or the `*_prec` columns) already makes `MDF` use it.
- `conf = NULL` in `MDF()`, `goFlux()` and `best.flux()`: when given, MDF = z
  x prec / t x flux.term with z = qnorm(1 - (1 - conf)/2) (1.96 for 0.95, the
  2-sigma criterion), following Wassmann et al. (2018). `best.flux(conf =)`
  scales the `MDF` column before the MDF criterion and records `MDF.conf`.
  With `conf` absent nothing changes (z = 1): all existing outputs are
  byte-identical.

Tests: `tests/testthat/test-empirical.prec.R` (recovery of known noise levels,
default-unchanged and scaling checks).
