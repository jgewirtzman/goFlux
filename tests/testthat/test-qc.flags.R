test_that("qc.flags appends flags only and co2.tracer joins on UniqueID", {
  data(manID.UGGA)
  best <- best.flux(goFlux(manID.UGGA, "CH4dry_ppb"))
  q <- qc.flags(best, dataframe = manID.UGGA, gastype = "CH4dry_ppb")
  expect_equal(nrow(q), nrow(best))
  for (nm in names(best)) expect_identical(q[[nm]], best[[nm]], info = nm)
  expect_true(all(c("qc.c0", "qc.convex", "qc.min.obs", "qc.ambient", "qc.clock", "qc.noisy", "qc.any") %in% names(q)))
  expect_false(q$qc.clock)   # the window starts at the recorded start
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

test_that("qc.flags: ambient start is judged at sealing, pre-closure drift is reported only, clock mismatch gives NA", {
  t0 <- as.POSIXct("2024-06-01 12:00:00", tz = "UTC")
  # 90 s of ambient, sealed at t0 + 90, dead band 30 s, window 120 s; 1 ppb noise
  mk <- function(id, seal.excess = 0, pre.slope = 0, win.start = 30, rise = 0.5) {
    tt <- 0:239
    amb <- 2000 + pre.slope * (tt - 90) * (tt < 90)
    x <- ifelse(tt < 90, amb, 2000 + seal.excess + rise * (tt - 90)) + stats::rnorm(240, 0, 1)
    data.frame(UniqueID = id, POSIX.time = t0 + 3600 * match(id, letters) + tt,
               start.time = t0 + 3600 * match(id, letters) + 90,
               flag = as.numeric(tt >= 90 + win.start), Etime = tt - 90 - win.start, CH4dry_ppb = x)
  }
  set.seed(21)
  d <- rbind(mk("a"),                         # clean: enriched by the window start, not at sealing
             mk("b", seal.excess = 40),       # headspace already 40 ppb above ambient at sealing
             mk("c", pre.slope = 0.3),        # ambient ramps up 18 ppb over the minute before
             mk("d", win.start = -20),        # window begins 20 s before the recorded start
             mk("e", rise = 2))               # steep: rises ~9 ppb within the first 10 s
  fx <- data.frame(UniqueID = c("a", "b", "c", "d", "e"))
  q <- qc.flags(fx, dataframe = d, gastype = "CH4dry_ppb")
  expect_equal(q$qc.ambient, c(FALSE, TRUE, FALSE, NA, TRUE))
  expect_equal(q$qc.clock, c(FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_true(is.na(q$qc.ambient.dev[4]))
  expect_gt(q$qc.ambient.dev[2], 3)
  # drift is reported, not flagged
  expect_gt(abs(q$qc.ambient.drift[3]), 1)
  # the tolerance uses the group precision (~1 ppb here)
  expect_equal(q$qc.prec[1], 1, tolerance = 0.2)
  # the old rule (window start vs ambient) would have fired on the clean and steep closures
  w0 <- tapply(d$CH4dry_ppb[d$flag == 1 & d$Etime <= 10], d$UniqueID[d$flag == 1 & d$Etime <= 10], mean)
  expect_gt(unname(w0["a"]) - 2000, 3 * 1.5)
  # arguments
  q2 <- qc.flags(fx, dataframe = d, gastype = "CH4dry_ppb", ambient.secs = 5, ambient.pre = 30)
  expect_equal(q2$qc.ambient[2], TRUE)
  # a shorter sealing interval for the steep closure
  expect_false(qc.flags(fx, dataframe = d, gastype = "CH4dry_ppb", ambient.secs = 2)$qc.ambient[5])
  expect_error(qc.flags(fx, dataframe = d, gastype = "CH4dry_ppb", ambient.pre = 0), "ambient.pre")
  # without the seal-time column: not evaluated
  q3 <- qc.flags(fx, dataframe = d, gastype = "CH4dry_ppb", seal.time = "nope")
  expect_false("qc.ambient" %in% names(q3))
})

test_that("qc.flags: qc.noisy uses the per-closure second-difference precision", {
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  set.seed(22)
  d <- do.call(rbind, lapply(1:6, function(i) data.frame(UniqueID = paste0("m", i), flag = 1,
    POSIX.time = t0 + 3600 * i + 0:119, Etime = 0:119,
    CH4dry_ppb = 2000 + i * (0:119) + stats::rnorm(120, 0, if (i == 6) 4 else 1))))
  q <- qc.flags(data.frame(UniqueID = paste0("m", 1:6)), dataframe = d, gastype = "CH4dry_ppb")
  expect_equal(q$qc.noisy, c(rep(FALSE, 5), TRUE))
  expect_equal(q$qc.noisy.ratio[6], 4, tolerance = 0.3)
})

test_that("qc.flags: qc.leak from a measured leak rate or a blank-closure slope", {
  fx <- data.frame(UniqueID = c("a", "b"), flux.term = 0.7, MDF = c(0.001, 0.1), C0 = 2000, Ct = 2100)
  q <- qc.flags(fx, blank.slope = 0.01)
  expect_equal(q$qc.leak.flux, c(0.007, 0.007))
  expect_equal(q$qc.leak, c(TRUE, FALSE))
  q2 <- qc.flags(fx, leak.rate = 1e-3, blank.slope = 0.01)
  expect_equal(q2$qc.leak.flux, rep(max(1e-3 * 100 * 0.7, 0.007), 2))
  # det.MDF from flux.class() takes precedence over goFlux's MDF
  fx$det.MDF <- c(0.1, 0.001)
  expect_equal(qc.flags(fx, blank.slope = 0.01)$qc.leak, c(FALSE, TRUE))
  expect_true(qc.flags(fx, blank.slope = 0.01)$qc.any[2])
  expect_false("qc.leak" %in% names(qc.flags(fx)))
  # with the concentration data, the range of the flagged rows is used
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  d <- data.frame(UniqueID = rep(c("a", "b"), each = 100), flag = 1, Etime = rep(0:99, 2),
                  POSIX.time = t0 + c(0:99, 1000 + 0:99), CH4dry_ppb = c(2000 + 0:99, 2000 + 2 * (0:99)))
  q3 <- qc.flags(fx, dataframe = d, gastype = "CH4dry_ppb", leak.rate = 1e-3)
  expect_equal(q3$qc.leak.flux, 1e-3 * c(99, 198) * 0.7)
})

test_that("qc.flags: seal.time defaults to cham.close when present", {
  t0 <- as.POSIXct("2024-06-01 12:00:00", tz = "UTC")
  set.seed(3)
  tt <- 0:239
  x <- ifelse(tt < 90, 2000, 2040 + 0.5 * (tt - 90)) + stats::rnorm(240, 0, 1)
  d <- data.frame(UniqueID = "a", POSIX.time = t0 + tt, cham.close = t0 + 90,
                  start.time = t0 + 120, flag = as.numeric(tt >= 120), Etime = tt - 120, CH4dry_ppb = x)
  q <- qc.flags(data.frame(UniqueID = "a"), dataframe = d, gastype = "CH4dry_ppb")
  expect_true(q$qc.ambient)           # judged at cham.close (90 s), where the 40 ppb excess appears
})
