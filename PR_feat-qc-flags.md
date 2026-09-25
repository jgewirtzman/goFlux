# feat: qc.flags() post-hoc quality flags and co2.tracer()

**Branch:** `feat/qc-flags`
**Files:** `R/qc.flags.R` (new: `qc.flags()`, `co2.tracer()`), `R/empirical.prec.R` (from `feat/empirical-precision-mdf`, needed for the record precision; rebase after that branch), `man/`, `NAMESPACE`, `tests/testthat/test-qc.flags.R`

`qc.flags(flux.result, dataframe = NULL, gastype = NULL, c0.mult = 1.5, by =
NULL, min.obs = 60, ambient.sigma = 3, ambient.secs = 10, noisy.mult = 1.5)`
appends independent logical flags to a `best.flux()` table and never modifies
or removes anything: `qc.c0` (C0 above c0.mult x the group median: chamber not
at ambient at closure), `qc.convex` (HM.k < 0), `qc.min.obs`, `qc.ambient`
(window start already enriched vs the pre-closure rows), `qc.noisy` (closure
MAD precision above noisy.mult x the record precision) and `qc.any`. Each is
a toggle (`NULL` skips it); the two that need the concentration data are
skipped when `dataframe` is absent.

`co2.tracer(co2.flux.result, flux.result, p.val = 0.05)` is deliberately a
one-liner: a logical aligned with `flux.result` saying whether the CO2
`best.flux` of the same UniqueID is positive with `LM.p.val < p.val`, i.e.
whether the chamber was sealed on a respiring surface. It is a comparison on
outputs, not a screen, and the docs say not to use it on open water, lit
foliage or dead wood.

Based on `chore/testthat-skeleton` (tests/testthat.R, `testthat` in Suggests).

## Follow-up (candidates for the upstream PR)

- `R/empirical.prec.R` synced from `feat/empirical-precision-mdf`:
  `method = "mad"` now returns one row per logging interval (`dt_s`) with a
  warning when a record changed interval.
- `qc.flags(noisy.mult =)` accordingly compares each closure's MAD precision
  with the record (or group) precision *at the closure's own logging interval*
  (median spacing of its flagged rows, rounded to 0.5 s), falling back to the
  interval with the most differences; it no longer assumes `empirical.prec()`
  returns a single row. Test added with a 1 Hz + 5 s mixed record.
- Fix: `qc.flags(by =)` always failed with "'dataframe' must contain a column
  that matches 'by'" because the internal grouping column was named `.by`,
  which the `\\<...\\>` word-boundary check in `empirical.prec()` cannot
  match (`.` is not a word character). Renamed to `qc_by`; covered by the new
  test.
