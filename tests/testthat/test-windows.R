# ported from fluxqc 0.2.4 tests/testthat/test-windows.R
t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")

make_window <- function(uid = "A", n = 900, start = 180) {
  set.seed(1)
  data.frame(UniqueID = uid, POSIX.time = t0 + 0:(n - 1),
             CO2dry_ppm = c(rep(410, 200), 410 + 0.3 * (0:400), rep(530, n - 601)) + rnorm(n, 0, 0.3),
             CH4dry_ppb = 2000 + 0.01 * (0:(n - 1)) + rnorm(n, 0, 0.5),
             start.time = t0 + start, obs.length = 600)
}

test_that("windows.from.table() crops with crop.meas() and keeps the recorded closure times", {
  ow <- list(make_window("A"), make_window("B"))
  win <- data.frame(UniqueID = c("A", "B"),
                    start = c("2024-06-01 10:04:00", "2024-06-01 10:05:00"),
                    end = c("2024-06-01 10:08:00", "2024-06-01 10:05:30"))
  expect_warning(m <- windows.from.table(ow, win), "only 31")
  expect_equal(sort(unique(m$UniqueID)), c("A", "B"))
  a <- m[m$UniqueID == "A", ]
  expect_equal(sum(a$flag), 241)
  expect_equal(unique(a$start.time), t0 + 240)            # crop.meas: window start
  expect_equal(unique(a$end.time), t0 + 480)
  expect_equal(unique(a$cham.close), t0 + 180)            # recorded closure start
  expect_equal(unique(a$cham.open), t0 + 780)
  expect_equal(unique(a$obs.length), 240)
  expect_equal(min(a$Etime[a$flag == 1]), 0)
  # a window that starts before the recorded start (e.g. after a clock correction)
  m2 <- windows.from.table(ow[1], data.frame(UniqueID = "A", start = t0 + 60, end = t0 + 300))
  expect_equal(sum(m2$flag), 241)
  expect_equal(unique(m2$cham.close), t0 + 180)
  expect_equal(unique(m2$start.time), t0 + 60)
  # missing UniqueID is dropped with a message
  expect_message(m3 <- windows.from.table(ow, win[1, ]), "No window for B")
  expect_equal(unique(m3$UniqueID), "A")
  expect_error(windows.from.table(ow, data.frame(UniqueID = "A", start = t0 + 10, end = t0)), "after")
  expect_error(windows.from.table(ow, rbind(win, win)), "unique")
})

test_that("find.rise() locates the sustained CO2 rise", {
  d <- make_window()
  r <- find.rise(d$POSIX.time, d$CO2dry_ppm)
  expect_false(is.null(r))
  expect_lt(abs(as.numeric(r$start - (t0 + 200), units = "secs")), 40)
  expect_gt(r$dur, 300)
  expect_gt(r$dconc, 100)
  expect_null(find.rise(d$POSIX.time, rep(410, nrow(d)) + rnorm(nrow(d), 0, 0.3)))
})

test_that("auto.id.rise() returns the click.peak2 structure and a log", {
  ow <- list(make_window("A"), make_window("B"))
  ow[[2]]$CO2dry_ppm <- 410 + rnorm(900, 0, 0.3)   # no rise -> fallback
  ends <- data.frame(UniqueID = "B", end.time = "2024-06-01 10:07:00")
  dir <- tempfile("autoid_")
  m <- auto.id.rise(ow, end.time = ends, plot.dir = dir)
  lg <- attr(m, "log")
  expect_equal(lg$method, c("auto_rise", "notebook_fallback"))
  expect_true(all(lg$n.flag > 60))
  expect_equal(length(list.files(dir, pattern = "png$")), 2)
  expect_true(all(c("flag", "Etime", "start.time_corr", "end.time_corr", "obs.length_corr") %in% names(m)))
  # fallback window = start + 30 .. notebook end
  expect_equal(unique(m$start.time_corr[m$UniqueID == "B"]), t0 + 180 + 30)
  expect_equal(unique(m$end.time_corr[m$UniqueID == "B"]), t0 + 420)
  # rise window is capped at max.secs
  expect_lte(unique(m$obs.length_corr[m$UniqueID == "A"]), 600)
  expect_equal(min(m$Etime[m$flag == 1 & m$UniqueID == "A"]), 0)
})

test_that("find.rise() follows the logging interval (10 s) and auto.id.rise() warns on non-CO2 defaults", {
  set.seed(2)
  tt <- t0 + 10 * (0:119)
  co2 <- c(rep(410, 20), 410 + 3 * (0:59), rep(590, 40)) + rnorm(120, 0, 0.3)
  r <- find.rise(tt, co2, min.n = 60)
  expect_false(is.null(r))
  expect_gt(r$dur, 400)
  expect_null(find.rise(tt, co2, min.n = 60, gap.secs = 5))   # the earlier fixed 5 s limit
  ow <- list(make_window("A"))
  expect_warning(auto.id.rise(ow, gastype = "CH4dry_ppb"), "CO2 in ppm")
})
