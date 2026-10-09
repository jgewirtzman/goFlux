# ported from fluxqc 0.2.4 tests/testthat/test-flux_term.R
test_that("flux.term() is (V/A)(P/RT)(1 - H2O) in mol m-2", {
  expect_equal(flux.term(1, 101.3, 100, 20), 0.1 * 101300 / (8.314 * 293.15), tolerance = 1e-6)
  expect_equal(flux.term(1, 101.3, 100, 20), flux.term(1, 101.3, 100, 20, 0))
  expect_lt(flux.term(1, 101.3, 100, 20, H2O_mol = 0.02), flux.term(1, 101.3, 100, 20))
})

test_that("MDF() is z * p / t * flux.term", {
  expect_equal(MDF(1.1, 540, 2), 1.1 / 540 * 2)
  expect_equal(MDF(1.1, 540, 2, conf = 0.95), qnorm(0.975) * 1.1 / 540 * 2)
  expect_error(MDF(1, 1, 1, conf = 1.2), "conf")
})
