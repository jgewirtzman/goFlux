test_that("flux.asinh scales by the MDF, keeps sign and zero, and inverts", {
  f <- c(-0.5, -0.01, 0, 0.01, 0.5, 50)
  y <- flux.asinh(f, scale = 0.1)
  expect_equal(y, asinh(f / 0.1))
  expect_equal(flux.asinh(y, scale = 0.1, inverse = TRUE), f)
  expect_equal(sign(y), sign(f))
  expect_equal(flux.asinh(0.001, 0.1), 0.01, tolerance = 1e-3)     # nearly linear well below the MDF
  expect_equal(flux.asinh(1000, 0.1), log(2 * 1e4), tolerance = 1e-6)  # logarithmic far above it
  expect_error(flux.asinh(f), "scale")
  expect_error(flux.asinh(f, scale = c(1, 2)), "single positive")
  expect_error(flux.asinh(f, scale = -1), "single positive")
  expect_error(flux.asinh("a", scale = 1), "numeric")
})
