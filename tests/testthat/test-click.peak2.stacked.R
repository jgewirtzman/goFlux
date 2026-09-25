test_that("click.peak2 with gases= gives identical output to the single-panel call", {
  data(manID.UGGA)
  ow <- unname(split(manID.UGGA, manID.UGGA$UniqueID))
  expect_error(click.peak2(ow, gastype = "CH4dry_ppb", gases = "nope"), "gases")
  # identify() is interactive: stub it (and the pop-up devices) to return a
  # fixed window in both runs, drawing into a null pdf device instead
  run <- function(...) {
    pdf(NULL); on.exit(dev.off())
    testthat::with_mocked_bindings(
      suppressWarnings(click.peak2(ow, gastype = "CH4dry_ppb", sleep = 0,
                                   plot.lim = c(0, 1e6), ...)),
      identify = function(x, y, ...) c(20L, 200L),
      dev.new = function(...) invisible(NULL),
      dev.off = function(...) invisible(NULL),
      dev.flush = function(...) invisible(NULL),
      .package = "goFlux")
  }
  a <- run()
  b <- run(gases = c("CO2dry_ppm", "H2O_ppm"))
  expect_identical(a, b)
  expect_true(is.data.frame(a))
  expect_equal(sum(a$flag), 181 * length(ow))
})
