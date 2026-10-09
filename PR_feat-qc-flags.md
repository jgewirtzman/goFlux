# feat: qc.flags() post-hoc quality flags and co2.tracer()

**Branch:** `feat/qc-flags`
**Files:** `R/qc.flags.R` (new: `qc.flags()`, `co2.tracer()`), `R/empirical.prec.R` and `R/spec.at.interval.R` (identical copies from `feat/empirical-precision-mdf`, needed for the group precision; drop them when rebasing onto that branch), `man/`, `NAMESPACE`, `tests/testthat/test-qc.flags.R`

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

## Revision: ambient start judged at sealing; second-difference precision

- `R/empirical.prec.R` (and `R/spec.at.interval.R`, which its docs link to)
  synced from `feat/empirical-precision-mdf`: default precision is now the
  per-closure second-difference (one-sample Hadamard) sigma, group = median
  over closures, with the `ac1` / `zero.frac` / `d1c.ratio` checks. Not
  rebased onto that branch because the copy was added in this branch's first
  commit (add/add conflict); the files are byte-identical instead.
- **`qc.ambient` is judged at the moment of sealing**, not at the window
  start (which follows the dead-band transient and is legitimately enriched):
  median of the first `ambient.secs` (10) s after the recorded closure start
  (`seal.time`, default `"start.time"`) vs the median of the `ambient.pre`
  (60) s before it; fires when |difference| > `ambient.sigma` (3) x
  max(group sigma, normal-consistent MAD of the pre-closure record).
  `qc.ambient.dev` is the signed difference in units of that tolerance. The
  pre-closure trend (`qc.ambient.drift`, same units) is reported but never
  raises the flag. New arguments `ambient.pre = 60`, `seal.time =
  "start.time"`; the docs warn that `start.time_corr` (window start) and the
  `start.time` returned by `crop.meas()` (moved by the dead band) are not the
  sealing time.
- **`qc.clock`** (new): the fitting window starts more than one logging
  interval before the recorded start (field and analyzer clocks disagree);
  `qc.ambient` is then `NA` (not evaluated). Counted in `qc.any`.
- **`qc.noisy`** compares the closure's second-difference sigma with its
  group's median sigma at the closure's own logging interval (`qc.prec`, new
  column, is the group value used).
- **Docs**: a paragraph that none of these flags measures drift or slow leaks
  (invisible to a difference-based sigma; they look like a flux), with
  periodic blank closures on an inert surface as the remedy and the
  ambient-shoulder trend only as a flag.
- `manID.UGGA` CH4: `qc.ambient.dev` 2.36 (window start vs mean ambient, in closure Allan sigma) -> -1.81 (at sealing, in max(0.445 ppb, MAD of pre-closure)); `qc.noisy.ratio` 0.89 -> 1 (one closure: its own group).
- Tests: synthetic closures (clean with a 30 s dead band, 40 ppb enriched at
  sealing, ramping ambient, window before the recorded start, steep rise with
  `ambient.secs` 10 vs 2), the old window-start rule shown to fire on the
  clean closure, argument checks; `qc.noisy` on a closure with 4x noise;
  the `empirical.prec()` tests copied from the other branch.
