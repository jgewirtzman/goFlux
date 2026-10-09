# ported from fluxqc 0.2.4 tests/testthat/test-examples.R
keys <- c("emission_stem_semirigid", "emission_upland_stem", "below_detection_stem", "failed_closure_stem",
          "ambiguous_stem", "uptake_stem", "uptake_soil", "high_flux_wetland_stem", "ebullition_floating",
          "emission_li7810_stem")

test_that("example closures load, anonymized, with real geometry", {
  ex <- example.closures()
  expect_named(ex, c("data", "aux"))
  expect_equal(ex$aux$UniqueID, keys)
  expect_equal(sort(unique(ex$data$UniqueID)), sort(keys))
  need <- c("UniqueID", "POSIX.time", "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm", "start.time", "obs.length",
            "Area", "Vtot", "Tcham", "Pcham", "flag", "Etime", "start.time_corr", "end.time_corr",
            "obs.length_corr", "CH4_prec", "CO2_prec")
  expect_true(all(need %in% names(ex$data)))
  expect_true(all(format(ex$data$POSIX.time, "%Y") == "2000"))
  expect_equal(ex$aux$Area[1], 4254.121, tolerance = 1e-6)
  expect_equal(ex$aux$Vtot[1], 7.082698, tolerance = 1e-6)
  dg <- read.csv(system.file("extdata", "example_closures", "example_closures_diagnostics.csv",
                             package = "goFlux"), stringsAsFactors = FALSE)
  expect_equal(dg$key, keys)
  fl <- tapply(ex$data$flag, ex$data$UniqueID, sum)[keys]; nn <- table(ex$data$UniqueID)[keys]
  expect_equal(as.vector(fl), dg$n_flagged); expect_equal(as.vector(nn), dg$n_segment)
  a <- ex$data[ex$data$UniqueID == "emission_stem_semirigid", ]
  expect_equal(nrow(a), 1207); expect_equal(sum(a$flag), 770)
  expect_equal(as.numeric(difftime(a$start.time_corr[1], a$start.time[1], units = "secs")), 91)
  expect_equal(min(a$Etime[a$flag == 1]), 0)
  # continuity of every selected window
  cont <- t(sapply(split(ex$data, ex$data$UniqueID)[keys], function(x) {
    w <- x[x$flag == 1, ]; w <- w[order(w$POSIX.time), ]
    dt <- as.numeric(diff(w$POSIX.time), units = "secs"); dc <- abs(diff(w$CH4dry_ppb))
    c(gap = max(dt) / median(dt), step = max(dc) / mad(dc)) }))
  expect_true(all(cont[, "gap"] <= 2))
  expect_true(all(cont[keys != "ebullition_floating", "step"] <= 6))
  expect_gt(cont["ebullition_floating", "step"], 20)
  re <- example.closures("emission_stem_semirigid", dead.band = 60)$data
  expect_equal(sum(re$flag), sum(re$t_s >= 60))
  expect_equal(nrow(example.closures(c("uptake_soil", "uptake_stem"))$aux), 2)
  expect_error(example.closures("z"), "no example")
})

test_that("ebullition_floating carries one real bubble (goAquaFlux)", {
  e <- example.closures("ebullition_floating")$data
  invisible(capture.output(r <- suppressWarnings(goAquaFlux(e, "CH4dry_ppb"))))
  expect_equal(nrow(r$bubbles), 1L)
  expect_gt(r$flux_summary$flux_ebullition, 15)
})
