# feat: empirical.prec() and a confidence level for MDF

**Branch:** `feat/empirical-precision-mdf`
**Files:** `R/empirical.prec.R` (new), `R/MDF.R`, `R/goFlux.R`, `R/best.flux.R`, `man/`, `NAMESPACE`, `tests/testthat/test-empirical.prec.R` (new)

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

Based on `chore/testthat-skeleton` (tests/testthat.R, `testthat` in Suggests).

## Follow-up (candidates for the upstream PR)

- `empirical.prec(method = "mad")` now computes the MAD of first differences
  per run of constant logging interval when `POSIX.time` is present: a run
  breaks where the spacing changes by > `tol` (20 %), duplicate/backwards
  timestamps are ignored, runs shorter than 3 differences are dropped (a lone
  gap between files is not a run), differences are pooled per interval
  (rounded to 0.5 s). The interval is returned in a trailing `dt_s` column; a
  record (or `by` group) with more than one interval gets one row per interval
  and a warning. Single-interval output keeps its shape (`prec`, `n`) plus
  `dt_s`. Without `POSIX.time`: pooled MAD as before, `dt_s = NA`.
- `spec.at.interval(spec_1s, dt_s)` (new, exported): datasheet precision
  rescaled to a logging interval, `spec_1s / sqrt(dt_s)`; compare
  `empirical.prec()` of a 5-10 s record with this, not with the 1 s figure.
- `MDF()` docs: `t` is `max(Etime) + 1` in seconds from `POSIX.time` at any
  logging interval (never `nb.obs`); test added that a 5-s-subsampled trace
  gives the same `MDF * t / prec` as the 1 Hz trace. No code path used
  `nb.obs` as seconds.
- `MDF()` doc corrections: Christiansen et al. (2015) define MDF =
  precision / enclosure time (k = 1); their "SD x 3 x t99" is a GC method
  quantification limit (Corley 2003), not a per-closure MDF, so no 3 x t
  variant is attributed to them. Wassmann et al. (2018) used a fixed k = 3
  ("99 %") on the SD of replicate ambient analyses; `conf` is the z-scaled
  generalisation; Parkin, Venterea & Hargreaves (2012, J Environ Qual
  41:705-715) cited for alpha = 0.05.
- Tests: `tests/testthat/test-empirical.prec.R` (runs guard, mixed-interval
  warning, gap/duplicate handling, 5-s MDF invariance, `spec.at.interval`).
