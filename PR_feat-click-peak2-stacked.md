# feat: click.peak2(gases =) stacked read-only panels

**Branch:** `feat/click-peak2-stacked`
**Files:** `R/click.peak2.R`, `man/click.peak2.Rd`, `tests/testthat/test-click.peak2.stacked.R` (new)

When picking a CH4 window it helps to see CO2 (a rise confirms the seal) and
H2O (a step marks a chamber change) at the same time. `click.peak2()` gains
`gases = NULL`: a character vector of additional gas columns that are drawn as
read-only panels above the `gastype` panel, all sharing the time axis and the
blue nominal-window lines. The `gastype` panel is drawn last so
`graphics::identify()` still acts on it alone, and the returned list is
identical to the unstacked call (test stubs `identify` and the pop-up devices
and compares both runs). With `gases` absent nothing changes. Any number of
panels is accepted; columns missing from `ow.list` raise an error naming them.

Based on `chore/testthat-skeleton` (tests/testthat.R, `testthat` in Suggests).
