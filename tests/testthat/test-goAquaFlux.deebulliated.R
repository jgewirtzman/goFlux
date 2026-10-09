# Synthetic incubation: known diffusive slope plus two bubbles with
# re-equilibration transients. The de-ebulliated diffusive flux should recover
# the slope; the pre-bubble window sees only the observations before 60 s.
make_trace <- function(slope = 2, steps = c(400, 250), over = c(150, 80), tau = 12,
                       t.b = c(60, 200), n = 360, sd = 1.5, seed = 7) {
  set.seed(seed)
  t <- 0:(n - 1)
  conc <- 2000 + slope * t + rnorm(n, 0, sd)
  for (i in seq_along(t.b)) {
    a <- t >= t.b[i]
    conc[a] <- conc[a] + steps[i] + over[i] * exp(-(t[a] - t.b[i]) / tau)
  }
  data.frame(UniqueID = "syn", POSIX.time = as.POSIXct("2024-06-01 12:00:00", tz = "UTC") + t,
             start.time = as.POSIXct("2024-06-01 12:00:00", tz = "UTC"), obs.length = n,
             Etime = t, flag = 1, CH4dry_ppb = conc, CH4_prec = 1, H2O_ppm = 10000,
             Vtot = 10, Area = 1000, Pcham = 101.325, Tcham = 20)
}
ft <- goFlux:::flux.term(10, 101.325, 1000, 20, 10000 / 1e6)
run <- function(d, ...) suppressWarnings(goAquaFlux(d, "CH4dry_ppb", Vtot = 10, Area = 1000,
                                                    Pcham = 101.325, Tcham = 20, ...))

test_that("deebulliated window recovers a known diffusive slope with two bubbles", {
  d <- make_trace()
  pre <- run(d); deb <- run(d, diffusion.window = "deebulliated")
  expect_equal(pre$flux_summary$diffusive_window, "pre_bubble")
  expect_equal(deb$flux_summary$diffusive_window, "deebulliated")
  expect_gte(nrow(deb$bubbles), 2)
  expect_equal(deb$flux_summary$flux_diffusive, 2 * ft, tolerance = 0.05)
  expect_gt(deb$flux_summary$n_obs.diffusion, 300)
  expect_lt(pre$flux_summary$n_obs.diffusion, 65)
  # ebullition is unchanged by the diffusive window choice
  expect_equal(deb$flux_summary$flux_ebullition, pre$flux_summary$flux_ebullition)
})

test_that("deebulliated gives a diffusive flux when the bubble comes too early", {
  d <- make_trace(t.b = c(20, 200))
  pre <- run(d); deb <- run(d, diffusion.window = "deebulliated")
  expect_true(is.na(pre$flux_summary$flux_diffusive))
  expect_false(is.na(deb$flux_summary$flux_diffusive))
  expect_equal(deb$flux_summary$flux_diffusive, 2 * ft, tolerance = 0.05)
})

test_that("the de-ebulliated trace is returned with the rise samples flagged out", {
  d <- make_trace()
  pre <- run(d); deb <- run(d, diffusion.window = "deebulliated")
  expect_null(pre$deebulliated)
  expect_s3_class(deb$deebulliated, "data.frame")
  expect_equal(nrow(deb$deebulliated), nrow(d))
  expect_true(all(c("UniqueID", "Etime", "flag", "CH4dry_ppb") %in% names(deb$deebulliated)))
  expect_gt(sum(deb$deebulliated$flag == 0), 0)
  expect_equal(sum(deb$deebulliated$flag == 1), deb$flux_summary$n_obs.diffusion)
  # the excluded samples are those of the physical rises, around the bubble times
  ex <- deb$deebulliated$Etime[deb$deebulliated$flag == 0]
  expect_true(all(ex >= 40 & ex <= 230))
  # the remaining trace is a straight line of the known slope
  keep <- deb$deebulliated[deb$deebulliated$flag == 1, ]
  expect_equal(unname(coef(lm(CH4dry_ppb ~ Etime, keep))[2]), 2, tolerance = 0.02)
})

test_that("find.bubbles(second.pass = TRUE) resolves a 12 s pair that the single pass merges", {
  set.seed(3); t <- 0:359; conc <- 2000 + 2 * t + rnorm(360, 0, 1.5)
  for (tb in c(100, 112)) { a <- t >= tb; conc[a] <- conc[a] + 300 + 100 * exp(-(t[a] - tb) / 12) }
  df <- data.frame(Etime = t, CH4dry_ppb = conc)
  b1 <- goFlux:::find.bubbles(df, "CH4dry_ppb", window.size = 15, method = "diff")
  b2 <- goFlux:::find.bubbles(df, "CH4dry_ppb", window.size = 15, method = "diff", second.pass = TRUE)
  expect_equal(nrow(b1), 1); expect_equal(nrow(b2), 2)
  expect_equal(b2$t.step, c(100, 112)); expect_equal(b2$magnitude[1], 300, tolerance = 0.05)
  # default goAquaFlux window is now 15, and second.pass is reachable through bubble.args
  expect_equal(formals(goAquaFlux)$bubble.window.size, 15)
})

test_that("deebulliated diffusive flux carries SE and n, short windows do not fail, and flux.plot.aqua draws the trace", {
  d <- make_trace()
  deb <- run(d, diffusion.window = "deebulliated")
  expect_false(is.na(deb$flux_summary$SE_diffusive)); expect_gt(deb$flux_summary$n_obs.diffusion, 100)
  short <- d[d$Etime < 25, ]   # 25 observations: too few for find.bubbles()
  invisible(capture.output(expect_warning(
    r <- goAquaFlux(short, "CH4dry_ppb", Vtot = 10, Area = 1000, Pcham = 101.325, Tcham = 20,
                    diffusion.minimum_window = 20), "bubble detection skipped")))
  expect_false(is.na(r$flux_summary$flux_diffusive)); expect_null(r$bubbles)
  invisible(capture.output(expect_warning(
    goAquaFlux(d, "CH4dry_ppb", Vtot = 10, Area = 1000, Pcham = 101.325, Tcham = 20,
               bubble_gas = "CO2dry_ppm"), "developed for CH4")))
  skip_if_not_installed("ggplot2")
  p <- suppressWarnings(flux.plot.aqua(deb, d, "CH4dry_ppb"))
  expect_s3_class(p[[1]], "ggplot")
  expect_true(any(vapply(p[[1]]$layers, function(l) "colour" %in% names(l$mapping) &&
                           identical(rlang::eval_tidy(l$mapping$colour), "de-ebulliated"), logical(1))))
})
