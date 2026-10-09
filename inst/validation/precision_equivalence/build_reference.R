# =============================================================================
# build_reference.R -- frozen reference values for the precision equivalence
# test (tests/testthat/test-precision.equivalence.R).
#
# NOT run at install time or by the tests. Computes the per-closure
# second-difference precision of the bundled closures with the two
# implementations that empirical.prec() replaces, read from their source folders
# (environment variables FLUXQC_DIR and WHOLE_TREE_FLUX_DIR) without loading
# either as a package:
#   fluxqc 0.2.4 (commit ecc7710)  R/precision.R    precision_hadamard()
#   whole_tree_flux                precision_helpers.R  closure_noise()
# and writes reference.csv next to this file. Run from the root of the fork:
#   Rscript inst/validation/precision_equivalence/build_reference.R
# =============================================================================

FLUXQC <- file.path(path.expand(Sys.getenv("FLUXQC_DIR")), "R/precision.R")
WT     <- file.path(path.expand(Sys.getenv("WHOLE_TREE_FLUX_DIR")),
                    "data processing/goFlux_reprocessing/precision_helpers.R")
if (!file.exists(FLUXQC) || !file.exists(WT)) stop("set FLUXQC_DIR and WHOLE_TREE_FLUX_DIR")
OUT    <- "inst/validation/precision_equivalence/reference.csv"

fq <- new.env(); sys.source(FLUXQC, envir = fq)
wt <- new.env(); sys.source(WT, envir = wt)

# bundled closures: the ten example closures and the two manID data sets
ex <- read.csv("inst/extdata/example_closures/example_closures.csv", stringsAsFactors = FALSE)
ex$POSIX.time <- as.POSIXct(ex$POSIX.time, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
ex$source <- "example.closures"
load("data/manID.UGGA.RData"); load("data/manID.G2201i.RData")
pick <- function(d, src) data.frame(source = src, UniqueID = d$UniqueID, POSIX.time = d$POSIX.time,
                                    flag = d$flag, CO2dry_ppm = d$CO2dry_ppm, CH4dry_ppb = d$CH4dry_ppb)
all <- rbind(pick(ex, "example.closures"), pick(manID.UGGA, "manID.UGGA"), pick(manID.G2201i, "manID.G2201i"))
all <- all[!is.na(all$flag) & all$flag == 1, ]

rows <- list()
for (src in unique(all$source)) for (id in unique(all$UniqueID[all$source == src]))
  for (gas in c("CH4dry_ppb", "CO2dry_ppm")) {
    d <- all[all$source == src & all$UniqueID == id, ]
    d <- d[order(d$POSIX.time), ]
    h <- fq$hadamard_closure(d[[gas]], as.numeric(d$POSIX.time))
    w <- wt$closure_noise(d[[gas]], d$POSIX.time)
    rows[[length(rows) + 1]] <- data.frame(
      source = src, UniqueID = id, gastype = gas, n = nrow(d),
      fluxqc.prec = h$sigma, fluxqc.dt_s = h$dt_s, fluxqc.n.diff = h$n_diff2, fluxqc.ac1 = h$ac1,
      wt.prec = w$sigma, wt.dt_s = w$dt_s, wt.n.diff = w$n_diff2, wt.ac1 = w$ac1, wt.t_sec = w$t_sec,
      stringsAsFactors = FALSE)
  }
ref <- do.call(rbind, rows)
write.csv(ref, OUT, row.names = FALSE)
cat("wrote", nrow(ref), "rows to", OUT, "\n")
