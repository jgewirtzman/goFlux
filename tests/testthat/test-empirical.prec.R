test_that("empirical.prec: superseded first-difference MAD and Allan per UniqueID", {
  data(manID.UGGA)
  set.seed(1)
  x <- data.frame(UniqueID = rep(c("a", "b"), each = 600), flag = 1, camp = rep(c("x", "y"), each = 600),
                  CH4dry_ppb = c(2000 + 0.5 * (1:600) + rnorm(600, 0, 1.5), 2000 + rnorm(600, 0, 3)))
  m <- empirical.prec(x, "CH4dry_ppb", method = "mad_diff1")
  expect_equal(names(m), c("prec", "n", "dt_s"))
  expect_equal(m$n, 1199L)
  expect_true(is.na(m$dt_s))   # no POSIX.time: pooled MAD, interval unknown
  # the old name still reaches it (partial matching)
  expect_identical(empirical.prec(x, "CH4dry_ppb", method = "mad"), m)
  g <- empirical.prec(x, "CH4dry_ppb", method = "mad_diff1", by = "camp")
  expect_equal(g$camp, c("x", "y"))
  expect_equal(g$prec, c(1.5, 3), tolerance = 0.15)
  expect_equal(nrow(g), 2L)
  a <- empirical.prec(x, "CH4dry_ppb", method = "allan")
  expect_equal(a$UniqueID, c("a", "b"))
  expect_equal(a$prec, c(1.5, 3), tolerance = 0.1)
  expect_s3_class(empirical.prec(manID.UGGA, "CH4dry_ppb", method = "allan"), "data.frame")
  expect_error(empirical.prec(x, "nope"), "gastype")
})

# synthetic closures: nc closures of n rows logged every dt s, slopes per step,
# noise from noise(n)
sim.closures <- function(nc, n, dt = 1, slopes = stats::runif(nc, -1, 5),
                         noise = function(n) stats::rnorm(n, 0, 1)) {
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  do.call(rbind, lapply(seq_len(nc), function(i)
    data.frame(UniqueID = paste0("c", i), flag = 1, POSIX.time = t0 + 3600 * i + dt * (0:(n - 1)),
               CH4dry_ppb = 2000 + slopes[i] * (0:(n - 1)) + noise(n))))
}

test_that("empirical.prec (hadamard): white noise is recovered within 5 %, checks near their white-noise values", {
  set.seed(11)
  w <- sim.closures(40, 180)
  expect_silent(e <- empirical.prec(w, "CH4dry_ppb"))
  expect_equal(names(e), c("prec", "n.closures", "dt_s", "n.diff", "ac1", "zero.frac", "prec.d1c", "d1c.ratio"))
  expect_equal(e$prec, 1, tolerance = 0.05)
  expect_equal(e$n.closures, 40L); expect_equal(e$dt_s, 1); expect_equal(e$n.diff, 40L * 178L)
  expect_equal(e$ac1, -2/3, tolerance = 0.05)
  expect_equal(e$zero.frac, 0)
  expect_equal(e$d1c.ratio, 1, tolerance = 0.1)
  # the group value is the median of the per-closure values
  cl <- attr(e, "closures")
  expect_equal(nrow(cl), 40L)
  expect_equal(e$prec, stats::median(cl$prec))
  # by: one row per group; by = "UniqueID": one row per closure
  w$day <- ifelse(as.integer(sub("c", "", w$UniqueID)) <= 20, "d1", "d2")
  b <- empirical.prec(w, "CH4dry_ppb", by = "day")
  expect_equal(b$day, c("d1", "d2")); expect_equal(b$n.closures, c(20L, 20L))
  expect_equal(b$prec, c(1, 1), tolerance = 0.08)
  u <- empirical.prec(w, "CH4dry_ppb", by = "UniqueID")
  expect_equal(nrow(u), 40L)
  # a 10-s record: same noise per sample
  set.seed(12)
  expect_equal(empirical.prec(sim.closures(40, 60, dt = 10), "CH4dry_ppb")$prec, 1, tolerance = 0.05)
  # needs the fitted windows
  expect_error(empirical.prec(w[, c("POSIX.time", "CH4dry_ppb")], "CH4dry_ppb"), "UniqueID")
})

test_that("empirical.prec (hadamard): two closures with very different slopes give the noise, not the slope spread (fluxqc #1)", {
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  # noise-free toy case: 0 and 10 ppb per step
  toy <- data.frame(UniqueID = rep(c("a", "b"), c(31, 30)), flag = 1,
                    POSIX.time = t0 + 10 * (0:60), CH4dry_ppb = 2000 + c(0 * (1:31), 10 * (1:30)))
  expect_equal(empirical.prec(toy, "CH4dry_ppb", warn = FALSE)$prec, 0)
  expect_gt(empirical.prec(toy, "CH4dry_ppb", method = "mad_diff1")$prec, 5)
  # with 1 ppb noise, logged every 10 s
  set.seed(13)
  two <- sim.closures(2, 60, dt = 10, slopes = c(0, 10))
  e <- empirical.prec(two, "CH4dry_ppb")
  expect_equal(e$prec, 1, tolerance = 0.2)
  expect_gt(empirical.prec(two, "CH4dry_ppb", method = "mad_diff1")$prec, 3)
  # second differences do not span a gap or a change of interval within a closure
  g <- sim.closures(1, 120, slopes = 1)
  g$POSIX.time[61:120] <- g$POSIX.time[61:120] + 600     # a 10-min gap
  g$CH4dry_ppb[61:120] <- g$CH4dry_ppb[61:120] + 500     # and a level jump across it
  cl <- attr(empirical.prec(g, "CH4dry_ppb"), "closures")
  expect_equal(cl$n.diff, 116L)                           # 2 x (60 - 2)
  expect_equal(cl$prec, 1, tolerance = 0.25)
})

test_that("empirical.prec (hadamard): alternate-row logging and red noise trigger their warnings", {
  set.seed(14)
  a <- sim.closures(20, 180)
  i <- seq(2, nrow(a), by = 2)
  a$CH4dry_ppb[i] <- a$CH4dry_ppb[i - 1]                 # one gas per row: carried forward
  expect_warning(e <- empirical.prec(a, "CH4dry_ppb"), "alternate-row")
  expect_equal(e$zero.frac, 0.5, tolerance = 0.02)
  set.seed(15)
  red <- sim.closures(20, 180, noise = function(n)
    cumsum(as.numeric(stats::arima.sim(list(ar = 0.5), n))) + stats::rnorm(n, 0, 0.3))
  expect_warning(empirical.prec(red, "CH4dry_ppb"), "lag-1 autocorrelation")
  r <- suppressWarnings(empirical.prec(red, "CH4dry_ppb"))
  expect_gt(r$ac1, -0.5)
  expect_gt(r$d1c.ratio, 1.2)
  # warn = FALSE: same columns, no warnings
  expect_silent(empirical.prec(red, "CH4dry_ppb", warn = FALSE))
})

test_that("empirical.prec (hadamard): closures at two logging intervals give one row per interval", {
  data(manID.UGGA)
  u <- manID.UGGA
  u5 <- u[seq(1, nrow(u), by = 5), ]; u5$UniqueID <- paste0(u5$UniqueID, "_5s")
  both <- rbind(u, u5)
  expect_warning(e <- empirical.prec(both, "CH4dry_ppb"), "2 logging intervals")
  expect_equal(e$dt_s, c(1, 5))
  both$camp <- ifelse(grepl("_5s$", both$UniqueID), "b", "a")
  expect_silent(eb <- empirical.prec(both, "CH4dry_ppb", by = "camp"))
  expect_equal(eb$camp, c("a", "b")); expect_equal(eb$dt_s, c(1, 5))
})

test_that("empirical.prec (mad_diff1): MAD per run of constant logging interval", {
  set.seed(3)
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  # single interval (1 Hz): one row, dt_s = 1, shape as before plus dt_s
  one <- data.frame(POSIX.time = t0 + 0:599, CH4dry_ppb = 2000 + 0.1 * (0:599) + rnorm(600, 0, 1))
  r1 <- empirical.prec(one, "CH4dry_ppb", method = "mad_diff1")
  expect_equal(nrow(r1), 1L)
  expect_equal(r1$dt_s, 1)
  expect_equal(r1$prec, 1, tolerance = 0.15)
  expect_equal(r1$n, 599L)
  # interval changes 1 s -> 10 s mid-record: one row per interval + warning,
  # each with its own noise level
  mixed <- data.frame(POSIX.time = c(t0 + 0:599, t0 + 600 + 10 * (0:299)),
                      CH4dry_ppb = 2000 + c(rnorm(600, 0, 1), rnorm(300, 0, 0.4)))
  expect_warning(r2 <- empirical.prec(mixed, "CH4dry_ppb", method = "mad_diff1"), "2 logging intervals")
  expect_equal(nrow(r2), 2L)
  expect_equal(r2$dt_s, c(1, 10))
  expect_equal(r2$prec, c(1, 0.4), tolerance = 0.15)
  # the same, per group: the mixed group gets two rows, the clean one one row
  mixed$camp <- c(rep("a", 600), rep("a", 300))
  both <- rbind(cbind(one, camp = "b"), mixed)
  expect_warning(r3 <- empirical.prec(both, "CH4dry_ppb", method = "mad_diff1", by = "camp"), "group 'a'")
  expect_equal(r3$camp, c("a", "a", "b"))
  expect_equal(r3$dt_s, c(1, 10, 1))
  # a lone gap between files (one long difference) is not a run: still one row
  gap <- data.frame(POSIX.time = c(t0 + 0:299, t0 + 3600 + 0:299),
                    CH4dry_ppb = 2000 + rnorm(600, 0, 1))
  expect_silent(r4 <- empirical.prec(gap, "CH4dry_ppb", method = "mad_diff1"))
  expect_equal(nrow(r4), 1L); expect_equal(r4$dt_s, 1)
  expect_equal(r4$n, 598L)   # the 3300-s jump is dropped
  # duplicate and backwards timestamps are ignored, not treated as intervals
  dup <- one; dup$POSIX.time[c(100, 200)] <- dup$POSIX.time[c(99, 199)]
  expect_silent(r5 <- empirical.prec(dup, "CH4dry_ppb", method = "mad_diff1"))
  expect_equal(r5$dt_s, 1)
  # a single off-interval difference inside a run (a 4-s gap left by a duplicate
  # timestamp in a 5-s record) is not an interval of its own: still one row
  data(manID.UGGA)
  u5 <- manID.UGGA[seq(1, nrow(manID.UGGA), by = 5), ]
  expect_equal(sort(unique(diff(as.numeric(u5$POSIX.time)))), c(4, 5))
  expect_silent(r6 <- empirical.prec(u5, "CH4dry_ppb", method = "mad_diff1"))
  expect_equal(nrow(r6), 1L); expect_equal(r6$dt_s, 5)
  # fewer than 3 usable differences: NA
  expect_true(is.na(empirical.prec(one[1:3, ], "CH4dry_ppb", method = "mad_diff1")$prec))
  # tol
  expect_error(empirical.prec(one, "CH4dry_ppb", method = "mad_diff1", tol = -1), "tol")
})

test_that("spec.at.interval rescales a 1-s datasheet precision", {
  expect_equal(spec.at.interval(0.9, c(1, 4, 9)), c(0.9, 0.45, 0.3))
  expect_error(spec.at.interval(0.9, 0), "dt_s")
  expect_error(spec.at.interval("a", 1), "spec_1s")
})
