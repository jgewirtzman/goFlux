# ported from fluxqc 0.2.4 tests/testthat/test-mdf.R (flag_detection)
d <- example.closures("emission_stem_semirigid")$data
best <- suppressWarnings(best.flux(suppressWarnings(goFlux(d, "CH4dry_ppb"))))
ft <- best$flux.term
fl <- d[d$flag == 1, ]
t.s <- closure.time(as.numeric(fl$POSIX.time))

test_that("default: group second-difference precision, z = 1, t in seconds; only appends", {
  out <- flux.class(best, dataframe = d, gastype = "CH4dry_ppb")
  expect_setequal(setdiff(names(out), names(best)),
                  c("det.prec", "det.t", "det.MDF", "det.class", "det.prec.closure",
                    "det.reg.half", "det.reg.class", "det.reg.ac1"))
  for (nm in names(best)) expect_identical(out[[nm]], as.data.frame(best)[[nm]], info = nm)
  s.i <- hadamard.closure(fl$CH4dry_ppb, as.numeric(fl$POSIX.time))$prec
  expect_equal(out$det.prec.closure, s.i)
  expect_equal(out$det.prec, s.i)                       # one closure: the median is itself
  expect_equal(out$det.t, t.s)
  expect_equal(out$det.MDF, s.i / t.s * ft)
  o95 <- flux.class(best, dataframe = d, gastype = "CH4dry_ppb", conf = 0.95)
  expect_equal(o95$det.MDF, qnorm(0.975) * s.i / t.s * ft)
  expect_equal(out$det.class, "emission")
  # the same duration as goFlux's own MDF
  expect_equal(best$MDF, best$prec / out$det.t * ft)
})

test_that("group precision is the median over the closures of each group", {
  set.seed(4)
  d2 <- d; d2$UniqueID <- "second"; d2$CH4dry_ppb <- d2$CH4dry_ppb + rnorm(nrow(d2), 0, 5)
  both <- rbind(d, d2)
  b2 <- suppressWarnings(best.flux(suppressWarnings(goFlux(both, "CH4dry_ppb"))))
  b2$camp <- c("x", "y")[match(b2$UniqueID, c("emission_stem_semirigid", "second"))]
  h <- function(x) { x <- x[x$flag == 1, ]; hadamard.closure(x$CH4dry_ppb, as.numeric(x$POSIX.time))$prec }
  g <- flux.class(b2, dataframe = both, gastype = "CH4dry_ppb", by = "camp")
  expect_equal(g$det.prec[g$UniqueID == "emission_stem_semirigid"], h(d))
  expect_equal(g$det.prec[g$UniqueID == "second"], h(d2))
  u <- flux.class(b2, dataframe = both, gastype = "CH4dry_ppb")
  expect_equal(u$det.prec, rep(median(c(h(d), h(d2))), 2))
  expect_error(flux.class(b2, dataframe = both, gastype = "CH4dry_ppb", by = "nope"), "by")
})

test_that("closures with very different slopes do not inflate the precision (fluxqc issue #1)", {
  set.seed(5)
  t0 <- as.POSIXct("2024-06-01", tz = "UTC")
  mk <- function(id, slope, k) data.frame(UniqueID = id, POSIX.time = t0 + 3600 * k + 10 * (0:59),
                                          flag = 1, CH4dry_ppb = 2000 + slope * (0:59) + rnorm(60, 0, 1))
  tr <- rbind(mk("slow", 0, 1), mk("fast", 10, 2))
  fr <- data.frame(UniqueID = c("slow", "fast"), flux.term = 0.7, nb.obs = 60, best.flux = 0.1)
  h <- flux.class(fr, dataframe = tr, gastype = "CH4dry_ppb")
  expect_equal(h$det.prec[1], 1, tolerance = 0.2)
  expect_equal(h$det.t, c(600, 600))                     # 60 rows x 10 s, not nb.obs
  expect_equal(h$det.MDF, h$det.prec / 600 * 0.7)
})

test_that("supplied precision: number, vector, or data.frame by UniqueID", {
  fr <- data.frame(UniqueID = c("a", "b", "c", "d"), flux.term = 1, best.flux = c(0.2, -0.2, 0.01, NA))
  z <- qnorm(0.975)
  o <- flux.class(fr, prec = data.frame(UniqueID = c("d", "c", "b", "a"), prec = c(4.4, 4.4, 1.5, 1.5)), t = 100, conf = 0.95)
  expect_equal(o$det.MDF, z * c(1.5, 1.5, 4.4, 4.4) / 100)
  expect_equal(o$det.class, c("emission", "uptake", "below MDF", NA))
  expect_equal(flux.class(fr, prec = 2, t = 100, conf = 0.95)$det.MDF, rep(z * 2 / 100, 4))
  expect_equal(flux.class(fr, prec = 2, t = 100, conf = NULL)$det.MDF, rep(2 / 100, 4))
  expect_equal(flux.class(fr, prec = 1:4, t = 100, conf = 0.95)$det.prec, 1:4)
  expect_error(flux.class(fr, prec = 2), "'t'")
  expect_error(flux.class(fr, prec = 1:3, t = 100, conf = 0.95), "prec")
})

test_that("without prec or dataframe, goFlux's own MDF is used", {
  o <- flux.class(best)
  expect_equal(o$det.MDF, best$MDF)
  expect_equal(o$det.prec, best$prec)
  fake <- best; fake$best.flux <- -best$MDF / 2
  expect_equal(flux.class(fake)$det.class, "below MDF")
  fake$best.flux <- -1
  expect_equal(flux.class(fake)$det.class, "uptake")
})

test_that("the duration follows the logging interval", {
  d5 <- d; d5$POSIX.time <- d5$POSIX.time[1] + 5 * (seq_len(nrow(d5)) - 1)   # same samples every 5 s
  o5 <- flux.class(best, dataframe = d5, gastype = "CH4dry_ppb", prec = 1)
  expect_equal(o5$det.t, 5 * sum(d5$flag == 1))
})


test_that("regression interval: Cowan half-width with the MDF precision; residual lag-1 autocorrelation", {
  out <- flux.class(best, dataframe = d, gastype = "CH4dry_ppb", conf = 0.95)
  tm <- as.numeric(fl$POSIX.time); tm <- tm - tm[1]
  n <- length(tm)
  expect_equal(out$det.reg.half, qt(0.975, n - 2) * out$det.prec / sqrt(sum((tm - mean(tm))^2)) * ft)
  r <- residuals(lm(fl$CH4dry_ppb ~ tm))
  expect_equal(out$det.reg.ac1, cor(r[-1], r[-n]))
  expect_equal(out$det.reg.class, ifelse(out$best.flux > out$det.reg.half, "emission",
                                     ifelse(out$best.flux < -out$det.reg.half, "uptake", "not significant")))
  o90 <- flux.class(best, dataframe = d, gastype = "CH4dry_ppb", reg.conf = 0.90)
  expect_equal(o90$det.reg.half / out$det.reg.half, qt(0.95, n - 2) / qt(0.975, n - 2))
  # supplied precision is used for the interval too
  p2 <- flux.class(best, dataframe = d, gastype = "CH4dry_ppb", prec = 2)
  expect_equal(p2$det.reg.half, qt(0.975, n - 2) * 2 / sqrt(sum((tm - mean(tm))^2)) * ft)
  # no readings, no interval
  expect_false(any(c("det.reg.half", "det.reg.class", "det.reg.ac1") %in% names(flux.class(best, prec = 2, t = 180))))
  expect_error(flux.class(best, dataframe = d, gastype = "CH4dry_ppb", reg.conf = 1.2), "reg.conf")
})

test_that("white noise around a line: interval about sqrt(12/n) of the MDF, residual ac1 near 0", {
  set.seed(11)
  t0 <- as.POSIXct("2024-06-01", tz = "UTC"); n <- 600
  tr <- data.frame(UniqueID = "w", POSIX.time = t0 + 0:(n - 1), flag = 1,
                   CH4dry_ppb = 2000 + 0.01 * (0:(n - 1)) + rnorm(n, 0, 1))
  fr <- data.frame(UniqueID = "w", flux.term = 1, best.flux = 0.01)
  o <- flux.class(fr, dataframe = tr, gastype = "CH4dry_ppb", prec = 1, conf = 0.95)
  expect_equal(o$det.reg.half / o$det.MDF, sqrt(12 / n), tolerance = 0.02)
  expect_lt(abs(o$det.reg.ac1), 0.15)
})
