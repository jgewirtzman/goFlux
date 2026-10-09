# ported from fluxqc 0.2.4 tests/testthat/test-import.R
li7810 <- system.file("extdata", "LI7810", "LI7810.data", package = "goFlux")

test_that("import2RData.flat() imports a folder without writing RData", {
  skip_if(!nzchar(li7810), "example LI7810.data not installed")
  wd <- tempfile("wd_"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  expect_message(imp <- import2RData.flat(dirname(li7810), "LI7810", pattern = "\\.data$",
                                          timezone = "UTC"), "imported")
  expect_false(dir.exists(file.path(wd, "RData")))
  expect_equal(unique(imp$source_file), "LI7810.data")
  expect_false(is.unsorted(imp$POSIX.time))
  lst <- suppressMessages(import2RData.flat(dirname(li7810), "LI7810", pattern = "\\.data$",
                                            timezone = "UTC", merge = FALSE))
  expect_type(lst, "list")
  expect_error(import2RData.flat(dirname(li7810), "NOPE"), "no import function")
  expect_error(import2RData.flat(tempdir(), "LI7810", pattern = "zzz$"), "no files")
})

test_that("import2RData.flat() reports a failing file and continues", {
  skip_if(!nzchar(li7810), "example LI7810.data not installed")
  dir <- tempfile("mixed_"); dir.create(dir)
  file.copy(li7810, dir)
  writeLines("this is not an analyzer file", file.path(dir, "bad.data"))
  expect_message(imp <- import2RData.flat(dir, "LI7810", pattern = "\\.data$", timezone = "UTC"),
                 "1 failed")
  expect_equal(unique(imp$source_file), "LI7810.data")
  expect_equal(basename(attr(imp, "failed")), "bad.data")
})

test_that("import.LI7810(dates = ) subsets the bundled file", {
  skip_if(!nzchar(li7810), "example LI7810.data not installed")
  all <- import.LI7810(li7810, timezone = "UTC")
  sub <- import.LI7810(li7810, timezone = "UTC", dates = "2022-12-05")
  expect_equal(nrow(sub), nrow(all[format(all$POSIX.time, "%Y-%m-%d") == "2022-12-05", ]))
  expect_error(import.LI7810(li7810, timezone = "UTC", dates = "1999-01-01"))
})
