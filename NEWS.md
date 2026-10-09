# goFlux 0.5.0.9002 (fork additions, jgewirtzman/goFlux)

GitHub release v0.5.0.9002, Zenodo doi:10.5281/zenodo.23256675.

* `qc.flags()`: `qc.convex` could not fire in 0.5.0.9001. It tested only the
  Hutchinson-Mosier curvature (`HM.k < 0`), which `goFlux()` keeps at 0 or
  above unless `k.min < 0`. With `dataframe`, it now fits a quadratic to the
  rows with `flag == 1` of each closure and fires when the quadratic term has
  the sign of the overall slope and p < `convex.p` (default 0.05), as fluxqc
  did; `HM.k < 0` remains the fallback without `dataframe`. Fluxes, MDFs and
  detection classes are unaffected.
* `qc.flags(min.secs = )`: new `qc.min.secs` flag for closures shorter than a
  duration in seconds (`closure.time()`), whatever the logging interval.

# goFlux 0.5.0.9001 (fork additions, jgewirtzman/goFlux)

Additions to goFlux 0.5.0 (Qepanna/goFlux master, 3ab9c04). Cite this release
(GitHub release v0.5.0.9001, Zenodo doi:10.5281/zenodo.23254791) together with Rheault et al.
(2024); see `citation("goFlux")` and CITATION.cff.

## Changed behaviour

* `goFlux()`: the closure duration `t` used in the MDF and in the kappa-max
  bound of the HM model is the span of the retained window plus one logging
  interval (`closure.time()`), instead of `max(Etime) + 1`. Identical at 1 Hz;
  at 5 s and 10 s logging the MDF is smaller by (dt - 1)/t (2.2 % and
  5.0 % for a 180 s closure) and HM fits bounded by kappa-max can change.
  The same `t` is used by `goAquaFlux()`.
* `goAquaFlux()`: `bubble.window.size` default 15 (also upstream).
* `qc.flags()`: `seal.time` defaults to `cham.close` when present.

## New functions

* `empirical.prec()`: per-closure second-difference precision
  (MAD/sqrt(6), a robust median-based form of the one-sample Hadamard
  deviation) and its median over the closures of a group, per logging
  interval, with checks (`ac1`, `zero.frac`, `d1c.ratio`). `"allan"` and
  `"mad_diff1"` are kept as documented comparison methods.
* `spec.at.interval()`: datasheet precision at a logging interval (valid only
  when each logged value is an average).
* `closure.time()`: closure duration from `Etime`.
* `flux.class()`: emission / uptake / below MDF, with the precision,
  duration and MDF used; replaces fluxqc's `flag_detection()`. Like `MDF()`,
  it defaults to z = 1 (`conf = NULL`); pass `conf = 0.95` for z = 1.96.
* `qc.flags()`: post-hoc flags (`qc.c0`, `qc.convex`, `qc.min.obs`,
  `qc.ambient` at sealing, `qc.clock`, `qc.noisy`, `qc.leak`) and
  `co2.tracer()`.
* `find.clock.offset()`, `find.rise()`, `auto.id.rise()`,
  `windows.from.table()` (one `crop.meas()` call): window and clock utilities.
  `find.rise()` scales its gap limit with the logging interval (fluxqc's fixed
  5 s found no rise at 10 s logging); its thresholds default to CO2 in ppm.
* `process.fluxes()`, `write.outputs()`: goFlux -> best.flux -> flux.class ->
  qc.flags in one call, with every option recorded.
* `example.closures()`: eleven real, relabelled closures with their geometry
  (`inst/extdata/example_closures/`), including two floating-chamber bubbles,
  one with a transient overshoot; `goAquaFlux()`, `find.bubbles()` and
  `flux.plot.aqua()` have runnable examples on them.
* `flux.plot.aqua()`: with `diffusion.window = "deebulliated"`, the measured
  series in grey, the series minus the bubble steps in sky blue, the excluded
  rise samples as crosses and the step times as dashed lines; the bubble
  model is labelled "step" or "step + overshoot".
* `import2RData.flat()`: folder import that reports a failing file and
  continues.
* `MDF()` and `flux.term()` are exported; `MDF(conf = )` labels the
  multiplier z (1.96 for 0.95: a benchmark, not a calibrated 95 % test).
* `click.peak2(gases = )`: stacked read-only panels of other gases.
* `goAquaFlux(diffusion.window = "deebulliated")`, `find.bubbles(second.pass,
  settle.mult)`: diffusive flux on the de-ebulliated trace (validation in
  `inst/validation/`).
* `import.LI7810(dates = )`: subset a campaign-long `.data` file by date.

## Fixes

* `crop.meas(max.obs.length = "aux")` joined a NULL object.

## Documentation

* MDF, `best.flux()` p-value note, `empirical.prec()`, `flux.term()` (total
  volume, analyzer internal volume) and the website pages on precision and
  quality flags state what the MDF and the p-value measure (Cowan et al. 2025,
  doi:10.1111/ejss.70104).

## Provenance

Replaces fluxqc 0.2.4 (commit ecc7710), now retired: `flag_detection()` ->
`flux.class()`, `qc_screens()` -> `qc.flags()`, `find_clock_offset()` ->
`find.clock.offset()`, `find_rise()`/`auto_id_rise()` -> `find.rise()`/
`auto.id.rise()`, `windows_from_table()` -> `windows.from.table()`,
`process_fluxes()`/`write_outputs()` -> `process.fluxes()`/`write.outputs()`,
`fluxqc_examples()` -> `example.closures()`, `import2RData_flat()` ->
`import2RData.flat()`, `flux_term()` -> `flux.term()`, `closure_seconds()` ->
`closure.time()`, `precision_hadamard()`/`precision_campaign()` ->
`empirical.prec()`, `spec_at_interval()` -> `spec.at.interval()`,
`click_peak2_stacked()` -> `click.peak2(gases = )`,
`import_li7810_subset()` -> `import.LI7810(dates = )`. Not ported:
`mdf_multi()` and the `extra` comparison columns of `flag_detection()`,
`precision_rolling()`, `example_closure()`.
