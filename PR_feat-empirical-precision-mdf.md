# feat: empirical.prec() and a confidence level for MDF

**Branch:** `feat/empirical-precision-mdf`
**Files:** `R/empirical.prec.R` (new), `R/MDF.R`, `R/goFlux.R`, `R/best.flux.R`, `man/`, `NAMESPACE`, `tests/testthat/test-empirical.prec.R` (new)

Datasheet precision underestimates the noise of field-worn analyzers (2.5x on
a UGGA in our data), so the MDF built on it is too optimistic. This adds

- `empirical.prec(dataframe, gastype, method = c("hadamard", "allan",
  "mad_diff1"), by = NULL, tol = 0.2, warn = TRUE)`. The default
  `"hadamard"` works closure by closure on the flagged rows:
  sigma_i = MAD(second differences) / sqrt(6), MAD the normal-consistent
  `stats::mad` (1.4826, about the median). A straight line (the flux)
  cancels exactly in a second difference; white noise gives Var = 6 sigma^2.
  This is the one-sample Hadamard deviation (Riley 2008, NIST SP 1065,
  doi:10.6028/NIST.SP.1065). A group (`by`, e.g. analyzer x field day; all
  closures when `NULL`) gets the **median** of sigma_i. `"allan"` =
  SD(dx)/sqrt(2) per `UniqueID` on the flagged rows (unchanged).
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
  (rounded to 0.5 s); an interval with fewer than 3 pooled differences is
  dropped. The interval is returned in a trailing `dt_s` column; a
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

## Revision: second-difference estimator (replaces the whole-record first-difference MAD)

- **Default estimator changed** from `"mad"` (MAD of first differences over
  the whole record, centred on one median) to `"hadamard"` (above). The old
  estimator counts the spread of slopes *between* closures as noise: two
  noise-free closures rising 0 and 10 ppb per step give ~5-7 ppb. The
  inflation is ~1.0-1.2x at 1 Hz but 1.4-2.3x at 5-10 s logging
  (https://github.com/jgewirtzman/fluxqc/issues/1). In a bake-off on 12
  simulated field-day types and real records from 8 analyzer setups, the
  closure-median second-difference sigma was within 5 % of the truth except
  for red noise and alternate-row loggers, both caught by its checks.
- **Old estimator kept as `method = "mad_diff1"`**, documented as superseded
  and kept for comparison with earlier results; `method = "mad"` still
  reaches it by partial matching. Kept rather than dropped because it is the
  only option that runs on a whole imported record without `UniqueID`/`flag`,
  and keeping it costs no API complexity (one more `method` value, same code).
- **Logging interval**: second differences never span a change of interval,
  a gap or a duplicate/backwards timestamp (same run logic as before, now per
  closure); each closure uses differences at its own interval (median
  spacing); a group whose closures were logged at several intervals gets one
  row per interval and a warning.
- **Checks returned as columns, with warnings** (`warn = FALSE` silences):
  `ac1`, lag-1 autocorrelation of the second differences (white noise -2/3;
  warns above -0.5: red or drift-dominated noise); `zero.frac`, fraction of
  exactly-zero first differences (warns at >= 0.3: one gas per logged row,
  e.g. Picarro G4301, or quantized data); `prec.d1c` and `d1c.ratio`,
  sigma from first differences centred on each closure's own median and its
  ratio to `prec` (warns outside 0.8-1.2: trend leakage / non-white noise, or
  alternate-row data). Per-closure values in `attr(, "closures")`;
  `by = "UniqueID"` gives one row per closure.
- **Output columns** for `"hadamard"`: [`by`], `prec`, `n.closures`, `dt_s`,
  `n.diff`, `ac1`, `zero.frac`, `prec.d1c`, `d1c.ratio`. (`"allan"` and
  `"mad_diff1"` outputs unchanged.)
- **Docs**: Riley (2008) reference; the alternate-row warning; a paragraph
  that any difference-based sigma is blind to drift and slow leaks, which look
  like a flux, with periodic blank closures on an inert surface as the
  remedy (ambient-shoulder trends only flag); a paragraph on comparing with
  the datasheet at the logging interval (spec / sqrt(dt), `spec.at.interval()`).
- `manID.UGGA` (1 Hz, one closure): CH4 0.556 (`mad_diff1`, whole record
  incl. shoulders) -> 0.445 ppb (`hadamard`, window); CO2 0.264 -> 0.155 ppm.
- Tests: white noise recovered within 5 % (40 closures x 180 s; also 10 s
  logging); the two-slope toy case (0 vs `mad_diff1` > 5) and its noisy
  10-s version (~1 ppb vs > 3); a gap with a level jump inside a closure is
  not differenced across; alternate-row synthetic triggers the zero-difference
  warning; integrated-AR(1) red noise triggers the lag-1 warning; mixed
  intervals give one row per interval.
