test_that("empirical.prec: MAD over the record and Allan per UniqueID", {
  data(manID.UGGA)
  set.seed(1)
  x <- data.frame(UniqueID = rep(c("a", "b"), each = 600), flag = 1, camp = rep(c("x", "y"), each = 600),
                  CH4dry_ppb = c(2000 + 0.5 * (1:600) + rnorm(600, 0, 1.5), 2000 + rnorm(600, 0, 3)))
  m <- empirical.prec(x, "CH4dry_ppb")
  expect_equal(names(m), c("prec", "n", "dt_s"))
  expect_equal(m$n, 1199L)
  expect_true(is.na(m$dt_s))   # no POSIX.time: pooled MAD, interval unknown
  g <- empirical.prec(x, "CH4dry_ppb", by = "camp")
  expect_equal(g$camp, c("x", "y"))
  expect_equal(g$prec, c(1.5, 3), tolerance = 0.15)
  expect_equal(nrow(g), 2L)
  a <- empirical.prec(x, "CH4dry_ppb", method = "allan")
  expect_equal(a$UniqueID, c("a", "b"))
  expect_equal(a$prec, c(1.5, 3), tolerance = 0.1)
  expect_s3_class(empirical.prec(manID.UGGA, "CH4dry_ppb", method = "allan"), "data.frame")
  expect_error(empirical.prec(x, "nope"), "gastype")
})

test_that("MDF and best.flux accept conf; default behaviour unchanged", {
  expect_equal(MDF(0.9, 300, 0.7), 0.9 / 300 * 0.7)
  expect_equal(MDF(0.9, 300, 0.7, conf = 0.95), qnorm(0.975) * 0.9 / 300 * 0.7)
  expect_error(MDF(0.9, 300, 0.7, conf = 2), "conf")
  data(manID.UGGA)
  r0 <- goFlux(manID.UGGA, "CH4dry_ppb")
  r1 <- goFlux(manID.UGGA, "CH4dry_ppb", conf = 0.95)
  expect_equal(r1$MDF, r0$MDF * qnorm(0.975))
  b0 <- best.flux(r0); b1 <- best.flux(r0, conf = 0.95)
  expect_equal(b1$MDF, b0$MDF * qnorm(0.975)); expect_equal(unique(b1$MDF.conf), 0.95)
  expect_false("MDF.conf" %in% names(b0))
})

test_that("empirical.prec: MAD per run of constant logging interval", {
  set.seed(3)
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  # single interval (1 Hz): one row, dt_s = 1, shape as before plus dt_s
  one <- data.frame(POSIX.time = t0 + 0:599, CH4dry_ppb = 2000 + 0.1 * (0:599) + rnorm(600, 0, 1))
  r1 <- empirical.prec(one, "CH4dry_ppb")
  expect_equal(nrow(r1), 1L)
  expect_equal(r1$dt_s, 1)
  expect_equal(r1$prec, 1, tolerance = 0.15)
  expect_equal(r1$n, 599L)
  # interval changes 1 s -> 10 s mid-record: one row per interval + warning,
  # each with its own noise level
  mixed <- data.frame(POSIX.time = c(t0 + 0:599, t0 + 600 + 10 * (0:299)),
                      CH4dry_ppb = 2000 + c(rnorm(600, 0, 1), rnorm(300, 0, 0.4)))
  expect_warning(r2 <- empirical.prec(mixed, "CH4dry_ppb"), "2 logging intervals")
  expect_equal(nrow(r2), 2L)
  expect_equal(r2$dt_s, c(1, 10))
  expect_equal(r2$prec, c(1, 0.4), tolerance = 0.15)
  # the same, per group: the mixed group gets two rows, the clean one one row
  mixed$camp <- c(rep("a", 600), rep("a", 300))
  both <- rbind(cbind(one, camp = "b"), mixed)
  expect_warning(r3 <- empirical.prec(both, "CH4dry_ppb", by = "camp"), "group 'a'")
  expect_equal(r3$camp, c("a", "a", "b"))
  expect_equal(r3$dt_s, c(1, 10, 1))
  # a lone gap between files (one long difference) is not a run: still one row
  gap <- data.frame(POSIX.time = c(t0 + 0:299, t0 + 3600 + 0:299),
                    CH4dry_ppb = 2000 + rnorm(600, 0, 1))
  expect_silent(r4 <- empirical.prec(gap, "CH4dry_ppb"))
  expect_equal(nrow(r4), 1L); expect_equal(r4$dt_s, 1)
  expect_equal(r4$n, 598L)   # the 3300-s jump is dropped
  # duplicate and backwards timestamps are ignored, not treated as intervals
  dup <- one; dup$POSIX.time[c(100, 200)] <- dup$POSIX.time[c(99, 199)]
  expect_silent(r5 <- empirical.prec(dup, "CH4dry_ppb"))
  expect_equal(r5$dt_s, 1)
  # a single off-interval difference inside a run (a 4-s gap left by a duplicate
  # timestamp in a 5-s record) is not an interval of its own: still one row
  data(manID.UGGA)
  u5 <- manID.UGGA[seq(1, nrow(manID.UGGA), by = 5), ]
  expect_equal(sort(unique(diff(as.numeric(u5$POSIX.time)))), c(4, 5))
  expect_silent(r6 <- empirical.prec(u5, "CH4dry_ppb"))
  expect_equal(nrow(r6), 1L); expect_equal(r6$dt_s, 5)
  # fewer than 3 usable differences: NA
  expect_true(is.na(empirical.prec(one[1:3, ], "CH4dry_ppb")$prec))
  # tol
  expect_error(empirical.prec(one, "CH4dry_ppb", tol = -1), "tol")
})

test_that("MDF: t is elapsed seconds, so a 5-s-logged trace gives the same MDF*t/prec", {
  data(manID.UGGA)
  u <- manID.UGGA[manID.UGGA$UniqueID == unique(manID.UGGA$UniqueID)[1], ]
  u5 <- u[seq(1, nrow(u), by = 5), ]            # same closure, logged every 5 s
  r1 <- suppressWarnings(goFlux(u, "CH4dry_ppb"))
  r5 <- suppressWarnings(goFlux(u5, "CH4dry_ppb", warn.length = 10))
  expect_equal(r1$nb.obs, 180L); expect_equal(r5$nb.obs, 36L)
  # goFlux uses t = max(Etime) + 1 on the re-anchored flagged rows
  tt <- function(d) { e <- d$Etime[d$flag == 1]; max(e) - min(e) + 1 }
  expect_equal(r1$MDF * tt(u) / r1$prec, r1$flux.term)
  expect_equal(r5$MDF * tt(u5) / r5$prec, r5$flux.term)
  # identical up to the subsampled mean of Tcham/Pcham in flux.term
  expect_equal(r5$MDF * tt(u5) / r5$prec, r1$MDF * tt(u) / r1$prec, tolerance = 1e-3)
  # and NOT nb.obs: with t = nb.obs the 5-s MDF would be ~5x larger
  expect_equal(r5$MDF, r5$prec * r5$flux.term / tt(u5))
  expect_false(isTRUE(all.equal(r5$MDF, r5$prec * r5$flux.term / r5$nb.obs)))
  expect_lt(r5$MDF / r1$MDF, 1.1)
})

test_that("spec.at.interval rescales a 1-s datasheet precision", {
  expect_equal(spec.at.interval(0.9, c(1, 4, 9)), c(0.9, 0.45, 0.3))
  expect_error(spec.at.interval(0.9, 0), "dt_s")
  expect_error(spec.at.interval("a", 1), "spec_1s")
})
