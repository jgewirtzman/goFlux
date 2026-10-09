# ported from fluxqc 0.2.4 tests/testthat/test-process.R and test-examples.R
d <- example.closures("emission_stem_semirigid")$data
best.ref <- suppressWarnings(best.flux(suppressWarnings(goFlux(d, "CH4dry_ppb"))))

test_that("process.fluxes(): appends only, best.flux columns untouched", {
  res <- suppressWarnings(process.fluxes(d, "CH4dry_ppb", qc = list(min.obs = 120)))
  expect_named(res, c("fluxes", "dataframe", "gastype", "settings"))
  fx <- res$fluxes
  expect_equal(nrow(fx), 1)
  for (nm in names(best.ref)) expect_identical(fx[[nm]], as.data.frame(best.ref)[[nm]], info = nm)
  added <- setdiff(names(fx), names(best.ref))
  expect_true(all(grepl("^det\\.|^qc\\.", added)))
  expect_true(all(c("det.MDF", "det.class", "qc.c0", "qc.convex", "qc.any") %in% names(fx)))
  expect_equal(fx$det.class, "emission")
  expect_true(is.na(res$settings$conf)); expect_equal(res$settings$z, 1)
  r2 <- suppressWarnings(process.fluxes(d, "CH4dry_ppb", qc = FALSE, prec = 7))
  expect_false("qc.any" %in% names(r2$fluxes))
  expect_equal(r2$fluxes$det.prec, 7)
  expect_null(r2$settings$qc)
})

test_that("flux.plot() runs on process.fluxes()$fluxes", {
  res <- suppressWarnings(process.fluxes(d, "CH4dry_ppb", qc = FALSE))
  pl <- suppressWarnings(flux.plot(res$fluxes, res$dataframe, "CH4dry_ppb"))
  expect_length(pl, 1)
  expect_true(inherits(pl[[1]], "ggplot"))
})

test_that("write.outputs() writes the documented files", {
  res <- suppressWarnings(process.fluxes(d, "CH4dry_ppb"))
  dir <- tempfile("out_")
  p <- suppressWarnings(write.outputs(res, dir, plots = TRUE))
  expect_true(all(c("fluxes", "settings", "plots") %in% names(p)))
  expect_true(all(file.exists(p)))
  js <- jsonlite::read_json(file.path(dir, "settings.json"))
  expect_equal(js$gastype, "CH4dry_ppb")
  expect_equal(nrow(read.csv(file.path(dir, "fluxes.csv"))), 1)
  expect_error(write.outputs(list(), dir), "process.fluxes")
})

test_that("the teaching closures separate (guidelines Figure 4 logic)", {
  ex <- example.closures(c("emission_stem_semirigid", "below_detection_stem", "failed_closure_stem",
                           "ambiguous_stem", "uptake_stem", "uptake_soil", "emission_upland_stem"))
  co2 <- suppressWarnings(best.flux(suppressWarnings(goFlux(ex$data, "CO2dry_ppm"))))
  res <- suppressWarnings(process.fluxes(ex$data, "CH4dry_ppb", co2.flux.result = co2, conf = 0.95))
  fx <- res$fluxes; rownames(fx) <- fx$UniqueID
  expect_equal(fx["emission_stem_semirigid", "det.class"], "emission")
  expect_equal(fx["emission_stem_semirigid", "best.flux"], 0.254, tolerance = 0.02)
  expect_true(fx["emission_stem_semirigid", "co2.tracer"])
  expect_equal(fx["emission_upland_stem", "det.class"], "emission")
  expect_equal(fx["below_detection_stem", "det.class"], "below MDF")
  expect_true(fx["below_detection_stem", "co2.tracer"])
  expect_equal(fx["failed_closure_stem", "det.class"], "below MDF")
  expect_false(fx["failed_closure_stem", "co2.tracer"])
  expect_false(fx["ambiguous_stem", "co2.tracer"])
  expect_equal(fx["uptake_stem", "det.class"], "uptake")
  expect_equal(fx["uptake_soil", "det.class"], "uptake")
  expect_lt(fx["uptake_soil", "best.flux"], -1)
  expect_true(all(!is.na(fx$qc.noisy)))
})

test_that("the two other analyzer formats run through goFlux", {
  ex <- example.closures(c("high_flux_wetland_stem", "emission_li7810_stem"))
  res <- suppressWarnings(process.fluxes(ex$data, "CH4dry_ppb", qc = FALSE))
  fx <- res$fluxes; rownames(fx) <- fx$UniqueID
  expect_gt(fx["high_flux_wetland_stem", "best.flux"], 15)
  expect_equal(fx["emission_li7810_stem", "best.flux"], 0.653, tolerance = 0.05)
})
