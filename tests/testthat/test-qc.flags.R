test_that("qc.flags appends flags only and co2.tracer joins on UniqueID", {
  data(manID.UGGA)
  best <- best.flux(goFlux(manID.UGGA, "CH4dry_ppb"))
  q <- qc.flags(best, dataframe = manID.UGGA, gastype = "CH4dry_ppb")
  expect_equal(nrow(q), nrow(best))
  for (nm in names(best)) expect_identical(q[[nm]], best[[nm]], info = nm)
  expect_true(all(c("qc.c0", "qc.convex", "qc.min.obs", "qc.ambient", "qc.noisy", "qc.any") %in% names(q)))
  expect_true(is.logical(q$qc.any))
  # synthetic contaminated / noisy closures
  fx <- data.frame(UniqueID = c("a", "b", "c"), C0 = c(2000, 2100, 11000), HM.k = c(0.01, -0.01, 0.01), nb.obs = c(300, 30, 300))
  q2 <- qc.flags(fx)
  expect_equal(q2$qc.c0, c(FALSE, FALSE, TRUE)); expect_equal(q2$qc.convex, c(FALSE, TRUE, FALSE))
  expect_equal(q2$qc.min.obs, c(FALSE, TRUE, FALSE)); expect_equal(q2$qc.any, c(FALSE, TRUE, TRUE))
  co2 <- best.flux(goFlux(manID.UGGA, "CO2dry_ppm"))
  tr <- co2.tracer(co2, best)
  expect_length(tr, nrow(best)); expect_true(is.logical(tr))
  expect_true(is.na(co2.tracer(co2, data.frame(UniqueID = "zzz"))))
})

test_that("qc.flags: qc.noisy compares each closure with the record precision of its own logging interval", {
  data(manID.UGGA)
  u <- manID.UGGA
  u5 <- u[seq(1, nrow(u), by = 5), ]; u5$UniqueID <- paste0(u5$UniqueID, "_5s")   # same closure logged every 5 s
  both <- rbind(u, u5)
  best <- suppressWarnings(best.flux(goFlux(both, "CH4dry_ppb", warn.length = 10)))
  # the record has two logging intervals: empirical.prec warns, qc.flags still works
  expect_warning(q <- qc.flags(best, dataframe = both, gastype = "CH4dry_ppb"), "2 logging intervals")
  expect_equal(nrow(q), 2L)
  expect_true(is.logical(q$qc.noisy)); expect_true(all(is.finite(q$qc.noisy.ratio)))
  # each closure vs the record precision at its own interval: neither is "noisy"
  expect_equal(q$qc.noisy, c(FALSE, FALSE))
  # 'by' path: one group per interval, no warning
  best$camp <- ifelse(grepl("_5s$", best$UniqueID), "b", "a")
  expect_silent(qb <- qc.flags(best, dataframe = both, gastype = "CH4dry_ppb", by = "camp"))
  expect_equal(qb$qc.noisy, c(FALSE, FALSE))
})
