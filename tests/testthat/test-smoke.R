test_that("goFlux() and best.flux() run on the bundled example data", {
  data(manID.UGGA)
  flux <- suppressWarnings(goFlux(manID.UGGA, "CH4dry_ppb"))
  expect_s3_class(flux, "data.frame")
  expect_true(all(c("UniqueID", "LM.flux", "HM.flux", "MDF", "prec", "flux.term") %in% names(flux)))
  best <- suppressWarnings(best.flux(flux))
  expect_equal(nrow(best), length(unique(manID.UGGA$UniqueID)))
  expect_true(all(is.finite(best$best.flux)))
})
