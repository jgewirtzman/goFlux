# Validation of goAquaFlux(diffusion.window = "deebulliated")
#
# Part A: synthetic traces with known diffusive slope and bubble events
#         (random factorial design, fixed seed).
# Part B: BlueFlux floating-chamber placements (the campaign's own
#         segmentation, geometry and fixed-threshold ebullition flags).
#
# Run from the root of the goFlux source tree on the branch
# feat/aqua-diffusive-deebulliated:
#   Rscript inst/validation/deebulliated_validation.R [N_SYNTH] [BLUEFLUX_ROOT]
# Outputs (in inst/validation/): results_synthetic.csv.gz, summary_synthetic_*.csv,
# summary_sensitivity.csv,
# results_blueflux.csv, events_blueflux.csv, fig_*.png.

args <- commandArgs(trailingOnly = TRUE)
N_SYNTH <- if (length(args) >= 1) as.integer(args[1]) else 2400L
BLUEFLUX <- if (length(args) >= 2) args[2] else
  "~/My Drive/Research/Blueflux/blueflux-ground"
OUT <- "inst/validation"
N_CORES <- max(1L, min(6L, parallel::detectCores() - 1L))

if (file.exists("DESCRIPTION") && any(grepl("^Package: goFlux", readLines("DESCRIPTION")))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else library(goFlux)
if (!"diffusion.window" %in% names(formals(goAquaFlux)))
  stop("this goFlux has no diffusion.window argument: run on the feature branch")

quiet <- function(expr) suppressWarnings(suppressMessages({
  invisible(utils::capture.output(res <- try(expr, silent = TRUE))); res }))

# Chamber used for every synthetic trace (10 L over 0.1 m2, 20 C, 1 atm)
VTOT <- 10; AREA <- 1000; PCH <- 101.325; TCH <- 20; H2O <- 10000
FT <- goFlux:::flux.term(VTOT, PCH, AREA, TCH, H2O / 1e6)  # nmol m-2 s-1 per ppb/s
K_HM <- 0.002                                              # curvature variant, s-1

# ============================================================================
# Part A: synthetic design
# ============================================================================
set.seed(20260924)
design <- data.frame(
  id      = seq_len(N_SYNTH),
  slope   = sample(c(-0.5, 0, 0.05, 0.5, 2, 10), N_SYNTH, TRUE),
  sigma   = sample(c(0.5, 2, 5), N_SYNTH, TRUE),
  noise   = sample(c("gaussian", "ar1"), N_SYNTH, TRUE, prob = c(0.65, 0.35)),
  L       = sample(c(120, 300, 600), N_SYNTH, TRUE),
  dt      = sample(c(1, 3, 5), N_SYNTH, TRUE),
  nb      = sample(0:4, N_SYNTH, TRUE, prob = c(0.25, 0.3, 0.25, 0, 0.2)),
  over    = sample(c(0, 0.2, 0.5), N_SYNTH, TRUE),
  tau     = sample(c(5, 12, 30), N_SYNTH, TRUE),
  ramp    = sample(c(1, 5, 15), N_SYNTH, TRUE),
  curv    = sample(c(FALSE, TRUE), N_SYNTH, TRUE, prob = c(0.75, 0.25)),
  stringsAsFactors = FALSE)
design$nb[design$nb == 3] <- 4L  # levels 0, 1, 2, 4

make_bubbles <- function(p) {
  if (p$nb == 0) return(NULL)
  tb <- runif(p$nb, 5, p$L - 15)
  if (runif(1) < 0.25) tb[1] <- runif(1, 5, 30)                # early bubble
  if (p$nb >= 2 && runif(1) < 0.30) tb[2] <- tb[1] + runif(1, 8, 19)  # close pair
  tb <- sort(pmin(tb, p$L - 10))
  data.frame(tb = tb, step = exp(runif(p$nb, log(5), log(2000))),
             r = p$ramp, over = p$over, tau = p$tau)
}

gen_trace <- function(p, bub) {
  t <- seq(0, p$L, by = p$dt); n <- length(t)
  base <- if (p$curv) 2000 + (p$slope / K_HM) * (1 - exp(-K_HM * t)) else 2000 + p$slope * t
  noise <- if (p$noise == "ar1") {
    phi <- 0.45
    as.numeric(stats::filter(rnorm(n, 0, p$sigma * sqrt(1 - phi^2)), phi, method = "recursive"))
  } else rnorm(n, 0, p$sigma)
  conc <- base + noise
  if (!is.null(bub)) for (j in seq_len(nrow(bub))) {
    b <- bub[j, ]; peak <- b$step * (1 + b$over)
    rise <- t >= b$tb & t < b$tb + b$r
    conc[rise] <- conc[rise] + peak * (t[rise] - b$tb) / b$r
    after <- t >= b$tb + b$r
    conc[after] <- conc[after] + b$step + b$step * b$over * exp(-(t[after] - b$tb - b$r) / b$tau)
  }
  data.frame(UniqueID = sprintf("syn%05d", p$id),
             POSIX.time = as.POSIXct("2024-06-01 12:00:00", tz = "UTC") + t,
             start.time = as.POSIXct("2024-06-01 12:00:00", tz = "UTC"), obs.length = p$L,
             Etime = t, flag = 1, CH4dry_ppb = conc, CH4_prec = 1, H2O_ppm = H2O,
             Vtot = VTOT, Area = AREA, Pcham = PCH, Tcham = TCH)
}

BWS <- NULL   # bubble.window.size; NULL = the package default (15 since this branch; 30 before)
bws_arg <- function() if (is.null(BWS)) list() else list(bubble.window.size = BWS)
run_aqua <- function(d, window) quiet(do.call(goAquaFlux, c(list(
  d, "CH4dry_ppb", Vtot = VTOT, Area = AREA, Pcham = PCH, Tcham = TCH,
  diffusion.window = window), bws_arg())))

match_events <- function(true_t, det, dt) {
  # true_t: peak times of the true bubbles; det: goFlux bubbles table
  tol <- max(10, 3 * dt)
  hit <- rep(FALSE, length(true_t)); used <- rep(FALSE, NROW(det)); idx <- rep(NA_integer_, length(true_t))
  if (NROW(det) > 0) for (j in order(true_t)) {
    tp <- if (!is.null(det$t.peak)) det$t.peak else det$start
    cand <- which(!used & (abs(tp - true_t[j]) <= tol | (det$start - 5 <= true_t[j] & det$end + 5 >= true_t[j])))
    if (length(cand)) { k <- cand[which.min(abs(tp[cand] - true_t[j]))]; hit[j] <- TRUE; used[k] <- TRUE; idx[j] <- k }
  }
  list(hit = hit, idx = idx, fp = sum(!used))
}

one_synth <- function(i) {
  p <- design[i, ]
  set.seed(1e6 + i)
  bub <- make_bubbles(p)
  d <- gen_trace(p, bub)
  true_diff <- p$slope * FT
  true_eb <- if (is.null(bub)) 0 else sum(bub$step) / p$L * FT
  naive <- unname(coef(lm(CH4dry_ppb ~ Etime, d))[2]) * FT
  r1 <- run_aqua(d, "pre_bubble"); r2 <- run_aqua(d, "deebulliated")
  out <- data.frame(p, n_true = NROW(bub), true_diff = true_diff, true_eb = true_eb,
                    first_tb = if (is.null(bub)) NA else bub$tb[1],
                    min_gap = if (NROW(bub) >= 2) min(diff(bub$tb)) else NA,
                    min_step = if (is.null(bub)) NA else min(bub$step),
                    naive = naive, pre = NA, deb = NA, total_pre = NA, total_deb = NA,
                    eb_est = NA, n_det = NA, hits = NA, fp = NA, mag_bias = NA, err_pre = NA, err_deb = NA)
  out$error <- inherits(r1, "try-error") || inherits(r2, "try-error")
  out$error_msg <- if (inherits(r1, "try-error")) as.character(r1) else if (inherits(r2, "try-error")) as.character(r2) else ""
  if (out$error) return(out)
  fs1 <- r1$flux_summary; fs2 <- r2$flux_summary
  out$pre <- fs1$flux_diffusive; out$deb <- fs2$flux_diffusive
  out$total_pre <- fs1$flux_total; out$total_deb <- fs2$flux_total
  out$eb_est <- fs1$flux_ebullition
  det <- r1$bubbles; out$n_det <- NROW(det)
  if (!is.null(bub)) {
    m <- match_events(bub$tb + bub$r, det, p$dt)
    out$hits <- sum(m$hit); out$fp <- m$fp
    if (any(m$hit)) out$mag_bias <- median((det$magnitude[m$idx[m$hit]] - bub$step[m$hit]) / bub$step[m$hit] * 100)
  } else { out$hits <- 0; out$fp <- NROW(det) }
  out
}

run_partA <- function(bws) {
  BWS <<- bws
  cat("Part A: ", N_SYNTH, "synthetic traces on", N_CORES, "cores, bubble.window.size =", if (is.null(bws)) "default" else bws, "\n")
  t0 <- Sys.time()
  syn <- do.call(rbind, parallel::mclapply(seq_len(N_SYNTH), one_synth, mc.cores = N_CORES))
  cat("  done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
  derive(syn)
}
derive <- function(syn) {
syn$err_naive <- syn$naive - syn$true_diff
syn$err_pre <- syn$pre - syn$true_diff
syn$err_deb <- syn$deb - syn$true_diff
rel <- function(err, tr) ifelse(tr != 0, err / abs(tr) * 100, NA)
syn$rel_naive <- rel(syn$err_naive, syn$true_diff)
syn$rel_pre <- rel(syn$err_pre, syn$true_diff)
syn$rel_deb <- rel(syn$err_deb, syn$true_diff)
syn$true_total <- syn$true_diff + syn$true_eb
syn$rel_total_pre <- rel(syn$total_pre - syn$true_total, syn$true_total)
syn$rel_total_deb <- rel(syn$total_deb - syn$true_total, syn$true_total)
syn$step_sigma <- syn$min_step / syn$sigma
syn$pre_na <- is.na(syn$pre)
syn$deb_na <- is.na(syn$deb)
syn$timing <- cut(syn$first_tb, c(0, 30, 60, 120, Inf), c("<30 s", "30-60 s", "60-120 s", ">120 s"))
syn$step_bin <- cut(syn$step_sigma, c(0, 3, 10, 30, 100, Inf), c("<3", "3-10", "10-30", "30-100", ">100"))
syn
}
syn30 <- run_partA(30)   # the old default, for the before/after table
syn <- run_partA(NULL)   # the current default (15)
con <- gzfile(file.path(OUT, "results_synthetic.csv.gz"), "w"); write.csv(syn, con, row.names = FALSE); close(con)

summ <- function(x, by) {
  g <- if (is.null(by)) factor(rep("all", nrow(x))) else interaction(x[by], drop = TRUE, sep = " | ")
  f <- function(v, fun) tapply(v, g, fun)
  med <- function(v) round(median(v, na.rm = TRUE), 1); mabs <- function(v) round(median(abs(v), na.rm = TRUE), 3)
  data.frame(group = levels(g), n = as.vector(table(g)),
             hit_rate = round(f(x$hits, sum) / pmax(f(x$n_true, sum), 1), 3),
             fp_per_trace = round(f(x$fp, mean), 3),
             mag_bias_pct = f(x$mag_bias, med),
             pre_NA = round(f(x$pre_na, mean), 3),
             rel_naive = f(x$rel_naive, med), rel_pre = f(x$rel_pre, med), rel_deb = f(x$rel_deb, med),
             abs_naive = f(x$err_naive, mabs), abs_pre = f(x$err_pre, mabs), abs_deb = f(x$err_deb, mabs),
             rel_total_pre = f(x$rel_total_pre, med), rel_total_deb = f(x$rel_total_deb, med),
             row.names = NULL)
}
ok <- syn[!syn$error, ]
with_b <- ok[ok$n_true > 0, ]
tabs <- list(
  overall = summ(ok, NULL), by_nb = summ(ok, "nb"), by_step = summ(with_b, "step_bin"),
  by_sigma_noise = summ(ok, c("sigma", "noise")), by_timing = summ(with_b, "timing"),
  by_dt = summ(with_b, "dt"), by_L = summ(with_b, "L"), by_ramp = summ(with_b, "ramp"),
  by_over_tau = summ(with_b, c("over", "tau")), by_slope = summ(ok, "slope"),
  by_curv = summ(ok, "curv"), by_noise = summ(ok, "noise"),
  close_pairs = summ(with_b[!is.na(with_b$min_gap), ], NULL),
  close_pairs_lt20 = summ(with_b[!is.na(with_b$min_gap) & with_b$min_gap < 20, ], NULL),
  bubble_free = summ(ok[ok$n_true == 0, ], c("sigma", "noise")))
for (nm in names(tabs)) write.csv(tabs[[nm]], file.path(OUT, paste0("summary_synthetic_", nm, ".csv")), row.names = FALSE)

# ---- figures A -------------------------------------------------------------
png(file.path(OUT, "fig1_step_sigma.png"), 1400, 600, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hb <- tabs$by_step
barplot(hb$hit_rate, names.arg = hb$group, ylim = c(0, 1), col = "grey70",
        xlab = "smallest step / sigma", ylab = "hit rate", main = "event detection")
boxplot(mag_bias ~ step_bin, with_b, outline = FALSE, xlab = "smallest step / sigma",
        ylab = "settled-step bias (%)", main = "magnitude", col = "grey90"); abline(h = 0, lty = 2)
dev.off()

png(file.path(OUT, "fig2_timing.png"), 1400, 600, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
tb <- with_b[with_b$slope != 0, ]
x <- data.frame(bias = c(tb$rel_naive, tb$rel_pre, tb$rel_deb),
                method = rep(c("naive LM", "pre_bubble", "deebulliated"), each = nrow(tb)),
                timing = rep(tb$timing, 3))
x$method <- factor(x$method, c("naive LM", "pre_bubble", "deebulliated"))
boxplot(bias ~ method + timing, x, outline = FALSE, las = 2, col = c("grey50", "grey75", "white"),
        ylab = "diffusive flux bias (%)", xlab = "", main = "bias vs first-bubble time", ylim = c(-100, 200),
        names = rep(c("naive", "pre", "deb"), 4)); abline(h = 0, lty = 2)
mtext(levels(x$timing), side = 3, at = c(2, 5, 8, 11), line = -1.2, cex = 0.8)
na <- tapply(with_b$pre_na, with_b$timing, mean)
barplot(na, ylim = c(0, 1), col = "grey70", ylab = "fraction NA", xlab = "first-bubble time",
        main = "pre_bubble returns NA")
dev.off()

png(file.path(OUT, "fig3_false_positives.png"), 1000, 600, res = 130)
bf <- ok[ok$n_true == 0, ]
fp <- tapply(bf$n_det > 0, list(bf$noise, bf$sigma), mean)
barplot(fp, beside = TRUE, ylim = c(0, max(0.3, max(fp, na.rm = TRUE) * 1.2)), col = c("grey40", "grey85"),
        xlab = "sigma (ppb)", ylab = "fraction of bubble-free traces with >= 1 event",
        legend.text = rownames(fp), main = "false positives on bubble-free traces")
dev.off()

png(file.path(OUT, "fig4_ar1_curvature.png"), 1400, 600, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (case in c("noise", "curv")) {
  s <- with_b[with_b$slope != 0, ]
  lab <- if (case == "noise") ifelse(s$noise == "ar1", "AR(1)", "Gaussian") else ifelse(s$curv, "HM curvature", "linear")
  x <- data.frame(bias = c(s$rel_pre, s$rel_deb), method = rep(c("pre", "deb"), each = nrow(s)), case = rep(lab, 2))
  boxplot(bias ~ method + case, x, outline = FALSE, col = c("grey75", "white"), ylim = c(-100, 150),
          ylab = "diffusive flux bias (%)", xlab = "", main = if (case == "noise") "noise model" else "diffusive shape")
  abline(h = 0, lty = 2)
}
dev.off()

# ---- before/after table (synthetic part) -----------------------------------
headline <- function(x) {
  x <- x[!x$error, ]; w <- x[x$n_true > 0, ]; b0 <- x[x$n_true == 0, ]; sl <- w[w$slope != 0, ]
  hr <- function(z) round(sum(z$hits) / sum(z$n_true), 3)
  medp <- function(v) round(median(v, na.rm = TRUE), 1)
  out <- c(hit_single = hr(w[w$nb == 1, ]), hit_multi = hr(w[w$nb >= 2, ]),
           hit_close_pairs_lt20s = hr(w[!is.na(w$min_gap) & w$min_gap < 20, ]),
           fp_bubble_free_gaussian = round(mean(b0$n_det[b0$noise == "gaussian"] > 0), 3),
           fp_bubble_free_ar1 = round(mean(b0$n_det[b0$noise == "ar1"] > 0), 3),
           step_bias_pct = medp(w$mag_bias))
  for (tm in levels(sl$timing)) out[paste0("pre_bias_", gsub(" ", "", tm))] <- medp(sl$rel_pre[sl$timing == tm])
  for (tm in levels(sl$timing)) out[paste0("deb_bias_", gsub(" ", "", tm))] <- medp(sl$rel_deb[sl$timing == tm])
  for (k in c(1, 2, 4)) out[paste0("pre_bias_nb", k)] <- medp(sl$rel_pre[sl$nb == k])
  for (k in c(1, 2, 4)) out[paste0("deb_bias_nb", k)] <- medp(sl$rel_deb[sl$nb == k])
  out["pre_NA_frac"] <- round(mean(w$pre_na), 3); out["deb_NA_frac"] <- round(mean(w$deb_na), 3)
  out
}
ba <- data.frame(metric = names(headline(syn30)), old_30 = unname(headline(syn30)), new_15 = unname(headline(syn)))
write.csv(ba, file.path(OUT, "summary_before_after.csv"), row.names = FALSE)

# ============================================================================
# Part B: BlueFlux floating-chamber placements
# ============================================================================
bf_file <- file.path(path.expand(BLUEFLUX), "output/ebullition/all_traces.rds")
if (file.exists(bf_file)) {
  cat("Part B: BlueFlux placements\n")
  tr <- readRDS(bf_file)
  AREA_B <- 324.3; VTOT_B <- 4.318; PCH_B <- 101.325; TCH_B <- 26   # campaign floating chamber (soil_water_dims.csv, auxfile)
  FT_B <- goFlux:::flux.term(VTOT_B, PCH_B, AREA_B, TCH_B, 0.01)

  cluster_jumps <- function(tj, gap = 10) {   # campaign point flags -> events
    if (!length(tj)) return(numeric(0)); tj <- sort(tj); ev <- tj[1]
    for (t in tj[-1]) if (t - ev[length(ev)] > gap) ev <- c(ev, t)
    ev
  }
  one_bf <- function(e) {
    s <- e$summary; d0 <- e$trace
    d0 <- d0[!is.na(d0$CH4_ppm), ]
    if (nrow(d0) < 40) return(NULL)
    d <- data.frame(UniqueID = s$placement_id, POSIX.time = d0$datetime,
                    start.time = d0$datetime[1], obs.length = max(d0$elapsed_sec),
                    Etime = d0$elapsed_sec, flag = 1, CH4dry_ppb = d0$CH4_ppm * 1000,
                    CH4_prec = 1, H2O_ppm = pmax(d0$H2O_ppm, 0),
                    Vtot = VTOT_B, Area = AREA_B, Pcham = PCH_B, Tcham = TCH_B)
    r1 <- quiet(do.call(goAquaFlux, c(list(d, "CH4dry_ppb", Vtot = VTOT_B, Area = AREA_B, Pcham = PCH_B, Tcham = TCH_B), bws_arg())))
    r2 <- quiet(do.call(goAquaFlux, c(list(d, "CH4dry_ppb", Vtot = VTOT_B, Area = AREA_B, Pcham = PCH_B, Tcham = TCH_B,
                           diffusion.window = "deebulliated"), bws_arg())))
    if (inherits(r1, "try-error") || inherits(r2, "try-error")) return(NULL)
    camp_ev <- cluster_jumps(d0$elapsed_sec[d0$is_jump %in% TRUE])
    det <- r1$bubbles
    m <- match_events(camp_ev, det, median(diff(d0$elapsed_sec)))
    naive <- unname(coef(lm(CH4dry_ppb ~ Etime, d))[2])
    list(row = data.frame(
      placement_id = s$placement_id, analyzer = s$analyzer, date = as.character(s$date), site = s$site,
      n_points = nrow(d0), duration = max(d0$elapsed_sec), dt = round(median(diff(d0$elapsed_sec)), 2),
      camp_n_jumps = sum(d0$is_jump %in% TRUE), camp_n_events = length(camp_ev),
      camp_total_ppb = s$total_ebullitive_ppm * 1000, camp_diff_ppb_s = s$diffusive_rate_clean_ppm_s * 1000,
      gof_n_events = NROW(det), matched = sum(m$hit), gof_unmatched = m$fp,
      gof_total_ppb = if (NROW(det)) sum(det$magnitude, na.rm = TRUE) else 0,
      gof_max_mag = if (NROW(det)) max(det$magnitude, na.rm = TRUE) else NA,
      slope_naive = naive,
      slope_pre = r1$flux_summary$flux_diffusive / FT_B, slope_deb = r2$flux_summary$flux_diffusive / FT_B,
      n_pre = r1$flux_summary$n_obs.diffusion, n_deb = r2$flux_summary$n_obs.diffusion,
      window_deb = r2$flux_summary$diffusive_window,
      flux_eb = r1$flux_summary$flux_ebullition, stringsAsFactors = FALSE),
      events = if (NROW(det)) cbind(placement_id = s$placement_id, det[, c("start", "end", "t.step", "t.peak", "magnitude", "overshoot", "tau", "SE", "reequil.complete")]) else NULL,
      trace = d, deb = r2$deebulliated, det = det, camp_ev = camp_ev)
  }
  run_partB <- function(bws) {
    BWS <<- bws
    bres <- lapply(tr, one_bf); bres <- bres[!sapply(bres, is.null)]
    bf <- do.call(rbind, lapply(bres, `[[`, "row"))
    bf$flagged <- bf$camp_n_jumps > 0
    # matched sample: every unflagged placement from the analyzer-days that have a flagged one
    fd <- unique(paste(bf$analyzer, bf$date)[bf$flagged])
    bf$matched_sample <- bf$flagged | paste(bf$analyzer, bf$date) %in% fd
    bf$slope_diff <- bf$slope_deb - bf$slope_pre
    list(bf = bf, bres = bres)
  }
  b30 <- run_partB(30); bf30 <- b30$bf
  bB <- run_partB(NULL); bf <- bB$bf; bres <- bB$bres
  write.csv(bf, file.path(OUT, "results_blueflux.csv"), row.names = FALSE)
  ev <- do.call(rbind, lapply(bres, `[[`, "events"))
  if (!is.null(ev)) write.csv(ev, file.path(OUT, "events_blueflux.csv"), row.names = FALSE)

  # figures B
  png(file.path(OUT, "fig5_blueflux.png"), 1400, 600, res = 130)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  rng <- range(c(bf$slope_pre, bf$slope_deb, bf$slope_naive), na.rm = TRUE)
  plot(bf$slope_pre, bf$slope_deb, xlim = rng, ylim = rng, pch = ifelse(bf$flagged, 16, 1),
       col = ifelse(bf$flagged, "red", "grey40"), xlab = "pre_bubble diffusive slope (ppb/s)",
       ylab = "deebulliated diffusive slope (ppb/s)", main = "BlueFlux placements", log = "")
  abline(0, 1, lty = 2); legend("topleft", c("campaign-flagged", "unflagged"), pch = c(16, 1), col = c("red", "grey40"), bty = "n")
  plot(bf$camp_n_events, bf$gof_n_events, pch = 16, col = adjustcolor("black", 0.4),
       xlab = "campaign events (0.10 ppm threshold, clustered)", ylab = "goFlux events", main = "detection agreement")
  abline(0, 1, lty = 2)
  dev.off()

  # page of traces: every campaign-flagged placement plus the unflagged ones
  # with the largest goFlux events (the candidate false positives)
  worst <- c(which(bf$flagged)[order(-bf$camp_n_events[bf$flagged])],
             head(which(!bf$flagged & bf$gof_n_events > 0)[order(-bf$gof_max_mag[!bf$flagged & bf$gof_n_events > 0])], 12 - sum(bf$flagged)))
  worst <- head(worst, 12)
  png(file.path(OUT, "fig6_worst_traces.png"), 1800, 2000, res = 120)
  par(mfrow = c(4, 3), mar = c(3.5, 4, 2.5, 1))
  for (i in worst) {
    b <- bres[[which(sapply(bres, function(z) z$row$placement_id) == bf$placement_id[i])]]
    d <- b$trace; plot(d$Etime, d$CH4dry_ppb, type = "l", col = "grey30", xlab = "Etime (s)", ylab = "CH4 (ppb)",
                       main = sprintf("%s\ncampaign %d / goFlux %d events; pre %.2f, deb %.2f ppb/s", bf$placement_id[i],
                                      bf$camp_n_events[i], bf$gof_n_events[i], bf$slope_pre[i], bf$slope_deb[i]), cex.main = 0.8)
    if (length(b$camp_ev)) abline(v = b$camp_ev, col = "red", lty = 3)
    if (NROW(b$det)) { rect(b$det$start, par("usr")[3], b$det$end, par("usr")[4], col = adjustcolor("orange", 0.2), border = NA)
      abline(v = b$det$t.peak, col = "orange") }
    if (!is.null(b$deb)) { points(b$deb$Etime, b$deb$CH4dry_ppb, col = ifelse(b$deb$flag == 1, "steelblue", "grey70"), pch = ".", cex = 2) }
  }
  dev.off()
  cat("  placements run:", nrow(bf), "of", length(tr), "\n")
  ba_blueflux <- function(b) {
    f <- b[b$flagged, ]; u <- b[!b$flagged, ]
    c(bf_events_on_flagged = sum(f$gof_n_events), bf_matched_campaign = sum(f$matched),
      bf_unflagged_with_events = sum(u$gof_n_events > 0), bf_unflagged_events_gt100ppb = sum(u$gof_max_mag > 100, na.rm = TRUE),
      bf_pre_NA = sum(is.na(b$slope_pre)), bf_deb_NA = sum(is.na(b$slope_deb)))
  }
  ba <- rbind(ba, data.frame(metric = names(ba_blueflux(bf30)), old_30 = unname(ba_blueflux(bf30)), new_15 = unname(ba_blueflux(bf))))
  write.csv(ba, file.path(OUT, "summary_before_after.csv"), row.names = FALSE)
  cat("\nBefore/after (bubble.window.size 30 vs 15):\n"); print(ba, row.names = FALSE)
} else cat("Part B skipped: ", bf_file, " not found\n")

# ============================================================================
# Part C: detector settings on the hard synthetic cases
# ============================================================================
cat("Part C: sensitivity of find.bubbles settings\n")
hard <- ok$id[ok$n_true >= 2 | ok$n_true == 0 | (ok$n_true == 1 & ok$ramp == 15)]
close_ids <- ok$id[!is.na(ok$min_gap) & ok$min_gap < 20]
BWS <- NULL
settings <- list(
  default_15         = list(),
  window30_old       = list(bubble.window.size = 30),
  min_gap5           = list(bubble.args = list(min_gap = 5)),
  k3                 = list(bubble.args = list(k = 3)),
  second_pass        = list(bubble.args = list(second.pass = TRUE)),
  settle1            = list(bubble.args = list(settle.mult = 1)),
  second_pass_settle1 = list(bubble.args = list(second.pass = TRUE, settle.mult = 1)))
one_sens <- function(i, st) {
  p <- design[i, ]; set.seed(1e6 + i); bub <- make_bubbles(p); d <- gen_trace(p, bub)
  r <- quiet(do.call(goAquaFlux, c(list(d, "CH4dry_ppb", Vtot = VTOT, Area = AREA, Pcham = PCH, Tcham = TCH,
                                        diffusion.window = "deebulliated"), st)))
  if (inherits(r, "try-error")) return(NULL)
  det <- r$bubbles; m <- if (is.null(bub)) list(hit = logical(0), idx = integer(0), fp = NROW(det)) else match_events(bub$tb + bub$r, det, p$dt)
  data.frame(id = i, n_true = NROW(bub), n_det = NROW(det), hits = sum(m$hit), fp = m$fp,
             mag_bias = if (any(m$hit)) median((det$magnitude[m$idx[m$hit]] - bub$step[m$hit]) / bub$step[m$hit] * 100) else NA,
             rel_deb = if (p$slope != 0) (r$flux_summary$flux_diffusive - p$slope * FT) / abs(p$slope * FT) * 100 else NA,
             eb_rel = if (!is.null(bub)) (r$flux_summary$flux_ebullition - sum(bub$step) / p$L * FT) / (sum(bub$step) / p$L * FT) * 100 else NA)
}
sens <- do.call(rbind, lapply(names(settings), function(nm) {
  r <- do.call(rbind, parallel::mclapply(hard, one_sens, st = settings[[nm]], mc.cores = N_CORES))
  wb <- r[r$n_true > 0, ]; bfree <- r[r$n_true == 0, ]
  data.frame(setting = nm, n = nrow(r),
             hit_rate = round(sum(wb$hits) / sum(wb$n_true), 3),
             hit_rate_multi = round(sum(wb$hits[wb$n_true >= 2]) / sum(wb$n_true[wb$n_true >= 2]), 3),
             hit_rate_close_lt20 = round(sum(wb$hits[wb$id %in% close_ids]) / sum(wb$n_true[wb$id %in% close_ids]), 3),
             fp_bubble_free = round(mean(bfree$n_det > 0), 3), fp_per_bubble_trace = round(mean(wb$fp), 3),
             mag_bias_pct = round(median(wb$mag_bias, na.rm = TRUE), 1),
             mag_abs_pct_q90 = round(quantile(abs(wb$mag_bias), 0.9, na.rm = TRUE), 1),
             eb_rel_pct = round(median(wb$eb_rel, na.rm = TRUE), 1),
             rel_deb_pct = round(median(wb$rel_deb, na.rm = TRUE), 1),
             rel_deb_abs_q90 = round(quantile(abs(wb$rel_deb), 0.9, na.rm = TRUE), 1))
}))
print(sens)
write.csv(sens, file.path(OUT, "summary_sensitivity.csv"), row.names = FALSE)
cat("done\n")
