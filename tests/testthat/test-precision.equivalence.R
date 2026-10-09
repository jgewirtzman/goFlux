# empirical.prec() against the implementations it replaces, on the bundled
# closures. Reference values are frozen in inst/validation/precision_equivalence/
# (build_reference.R): fluxqc 0.2.4 (ecc7710) precision_hadamard() and
# whole_tree_flux closure_noise().
ref <- read.csv(system.file("validation", "precision_equivalence", "reference.csv", package = "goFlux"),
                stringsAsFactors = FALSE)

fork_closures <- function() {
  ex <- example.closures()$data; ex$source <- "example.closures"
  data(manID.UGGA, manID.G2201i, envir = environment())
  pick <- function(d, src) data.frame(source = src, UniqueID = d$UniqueID, POSIX.time = d$POSIX.time,
                                      flag = d$flag, CO2dry_ppm = d$CO2dry_ppm, CH4dry_ppb = d$CH4dry_ppb)
  all <- rbind(pick(ex, "example.closures"), pick(manID.UGGA, "manID.UGGA"),
               pick(manID.G2201i, "manID.G2201i"))
  all$key <- paste(all$source, all$UniqueID)
  out <- lapply(c("CH4dry_ppb", "CO2dry_ppm"), function(gas) {
    cl <- attr(suppressWarnings(empirical.prec(transform(all, UniqueID = key), gas, by = "UniqueID",
                                               warn = FALSE)), "closures")
    data.frame(key = cl$UniqueID, gastype = gas, prec = cl$prec, dt_s = cl$dt_s,
               n.diff = cl$n.diff, ac1 = cl$ac1, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

test_that("empirical.prec() per closure equals fluxqc 0.2.4 precision_hadamard()", {
  f <- fork_closures()
  i <- match(paste(ref$source, ref$UniqueID, ref$gastype), paste(f$key, f$gastype))
  expect_false(anyNA(i))
  expect_equal(f$prec[i], ref$fluxqc.prec, tolerance = 1e-12)
  expect_equal(f$dt_s[i], ref$fluxqc.dt_s)
  expect_equal(f$n.diff[i], ref$fluxqc.n.diff)
  expect_equal(f$ac1[i], ref$fluxqc.ac1, tolerance = 1e-12)
})

test_that("empirical.prec() per closure equals whole_tree_flux closure_noise()", {
  f <- fork_closures()
  i <- match(paste(ref$source, ref$UniqueID, ref$gastype), paste(f$key, f$gastype))
  expect_equal(f$prec[i], ref$wt.prec, tolerance = 1e-12)
  expect_equal(f$dt_s[i], ref$wt.dt_s)
  expect_equal(f$n.diff[i], ref$wt.n.diff)
  expect_equal(f$ac1[i], ref$wt.ac1, tolerance = 1e-12)
})

test_that("closure.time() equals whole_tree_flux's span + interval", {
  ex <- example.closures()$data; ex <- ex[ex$flag == 1, ]
  r <- ref[ref$source == "example.closures" & ref$gastype == "CH4dry_ppb", ]
  tt <- vapply(split(as.numeric(ex$POSIX.time), ex$UniqueID), closure.time, numeric(1))
  expect_equal(unname(tt[r$UniqueID]), r$wt.t_sec)
})
