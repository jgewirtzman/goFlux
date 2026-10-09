test_that("closure.time() is the span plus one logging interval", {
  expect_equal(closure.time(0:179), 180)
  expect_equal(closure.time(seq(0, 175, by = 5)), 180)
  expect_equal(closure.time(seq(0, 170, by = 10)), 180)
  # duplicate timestamps and one gap do not change the interval
  expect_equal(closure.time(c(0, 5, 5, 10, 15, 25, 30)), 35)
  # unordered and NA input
  expect_equal(closure.time(c(10, 0, NA, 5)), 15)
  # a single timestamp falls back to a 1 s interval
  expect_equal(closure.time(0), 1)
  expect_error(closure.time("a"))
})

# Thin the 1 Hz example to a coarser logging interval, keeping the window.
thin_to <- function(d, dt) {
  keep <- d$flag == 0 | (d$flag == 1 & (d$Etime %% dt) == 0)
  d[keep, ]
}

test_that("goFlux() MDF uses span + interval at 1, 5 and 10 s logging", {
  data(manID.UGGA)
  for (dt in c(1, 5, 10)) {
    d <- thin_to(manID.UGGA, dt)
    flux <- suppressWarnings(goFlux(d, "CH4dry_ppb"))
    et <- d$Etime[d$flag == 1]
    # the example has one duplicate timestamp: 180 rows span 178 s at 1 Hz
    t_exp <- diff(range(et)) + dt
    expect_equal(flux$MDF, flux$prec / t_exp * flux$flux.term, tolerance = 1e-10,
                 info = paste("dt =", dt))
  }
})

test_that("at 1 Hz the duration equals the earlier max(Etime) + 1", {
  data(manID.UGGA)
  et <- manID.UGGA$Etime[manID.UGGA$flag == 1]
  expect_equal(closure.time(et - min(et)), max(et - min(et)) + 1)
})
