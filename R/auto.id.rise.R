#' Longest sustained concentration rise in a trace
#'
#' Scans a concentration series for the longest stretch during which a
#' chamber was evidently closed: the concentration rises by more than
#' \code{rise} over \code{rise.secs}, then keeps climbing (never dropping more
#' than \code{drop} below its running maximum) with no gap longer than
#' \code{gap.secs}. Used by \code{\link[goFlux]{auto.id.rise}}. The defaults
#' of \code{rise}, \code{drop} and \code{conc.range} are for CO2 in ppm;
#' set them for any other gas or unit.
#'
#' @param time POSIXct (or numeric, seconds) time stamps.
#' @param conc numerical; concentrations (ppm for CO2).
#' @param rise numerical; minimum increase over \code{rise.secs} that starts a
#'             candidate rise, units of \code{conc}. Default 6 (ppm CO2).
#' @param rise.secs numerical; look-ahead interval (s) for the start test.
#'                  Default 60.
#' @param drop numerical; tolerated drop below the running maximum before the
#'             rise ends, units of \code{conc}. Default 8 (ppm CO2).
#' @param min.dur numerical; minimum duration (s) of an accepted rise.
#'                Default 90.
#' @param gap.secs numerical; a gap between consecutive samples longer than
#'                 this (s) ends the rise. Default \code{NULL}: three times the
#'                 logging interval (median spacing of \code{time}), at least
#'                 5 s.
#' @param conc.range numerical vector of length 2; values outside this range
#'                   are discarded before scanning (e.g. analyzer error
#'                   values). Default \code{c(300, 20000)} (ppm CO2).
#' @param min.n numerical; minimum number of valid samples needed to scan.
#'              Default 120.
#'
#' @returns \code{NULL} if no rise is found; otherwise a list with
#'          \code{start}, \code{end} (same class as \code{time}), \code{dur}
#'          (seconds) and \code{dconc} (rise, units of \code{conc}).
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{auto.id.rise}}
#'
#' @examples
#' set.seed(1)
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' co2 <- c(rep(410, 200), 410 + 0.3 * (0:400), rep(530, 300)) + rnorm(901, 0, 0.3)
#' find.rise(t0 + 0:900, co2)
#' @export
find.rise <- function(time, conc, rise = 6, rise.secs = 60, drop = 8,
                      min.dur = 90, gap.secs = NULL, conc.range = c(300, 20000),
                      min.n = 120) {
  if(length(time) != length(conc)) stop("'time' and 'conc' must have the same length")
  ok <- !is.na(conc) & !is.na(time) & conc > conc.range[1] & conc < conc.range[2]
  time <- time[ok]; conc <- conc[ok]
  n <- length(time)
  if(n < min.n) return(NULL)
  ts <- as.numeric(time)
  if(is.null(gap.secs)){
    dt <- diff(ts); dt <- dt[dt > 0]
    gap.secs <- max(5, 3 * if(length(dt)) stats::median(dt) else 1)
  }
  best <- NULL
  i <- 1L
  while(i <= n){
    k <- findInterval(ts[i] + rise.secs, ts)   # last index with time <= t_i + rise.secs
    if(k <= i) k <- min(i + 1L, n)
    if(k >= n && ts[k] - ts[i] < rise.secs) break
    if(conc[k] - conc[i] > rise){
      j <- i; peak <- conc[i]
      while(j < n && conc[j + 1] >= peak - drop && (ts[j + 1] - ts[j]) < gap.secs){
        j <- j + 1L; peak <- max(peak, conc[j])
      }
      dur <- ts[j] - ts[i]
      if(dur >= min.dur && (is.null(best) || dur > best$dur)){
        best <- list(start = time[i], end = time[j], dur = dur, dconc = conc[j] - conc[i])}
      i <- j + 10L
    } else i <- i + 5L
  }
  best
}

#' Automatic observation windows from the sustained CO2 rise
#'
#' A non-interactive alternative to \code{\link[goFlux]{click.peak2}} that can
#' run under \code{Rscript}. For each observation window it finds the longest
#' sustained rise of \code{gastype} with \code{\link[goFlux]{find.rise}} (a
#' sealed chamber on a respiring surface accumulates CO2), removes the first
#' \code{trim.secs} (mixing), and caps the window at the notebook end time plus
#' \code{slack.secs} or at \code{max.secs}, whichever comes first. When no rise
#' is found it falls back to notebook start + \code{trim.secs} to the notebook
#' end (or \code{fallback.secs}). Optionally one diagnostic PNG per
#' measurement, with a stacked panel per gas, is written for review; windows
#' to redo can then be clicked with \code{\link[goFlux]{click.peak2}}.
#'
#' @param ow.list list of data.frames; output of \code{\link[goFlux]{obs.win}}
#'                (each with \code{UniqueID}, \code{POSIX.time} and
#'                \code{start.time}). A single data.frame is split by
#'                \code{UniqueID}.
#' @param gastype character string; the column used to find the rise.
#'                Default \code{"CO2dry_ppm"}.
#' @param gases character vector; the gas columns drawn in the diagnostic
#'              plots, one panel each. Default: \code{gastype} and, when
#'              present, \code{"CH4dry_ppb"} and \code{"N2Odry_ppb"}.
#' @param end.time optional notebook end times: a data.frame with the columns
#'                 \code{UniqueID} and \code{end.time} (POSIXct, or character
#'                 in \code{tz}), or a named POSIXct vector. Missing entries are
#'                 allowed.
#' @param trim.secs numerical; seconds removed from the start of the rise.
#'                  Default 30.
#' @param max.secs numerical; maximum window length (s). Default 600.
#' @param min.secs numerical; shorter windows are extended to
#'                 \code{min.secs + 30} s and reported as \code{"short"}.
#'                 Default 60.
#' @param slack.secs numerical; seconds allowed past the notebook end time.
#'                   Default 60.
#' @param fallback.secs numerical; window length used when neither a rise nor
#'                      a notebook end time is available. Default 300.
#' @param plot.dir character string; folder for the diagnostic PNGs.
#'                 Default \code{NULL}: no plots.
#' @param tz character string; time zone of character \code{end.time}.
#'           Default \code{"UTC"}.
#' @param warn.length numerical; windows with fewer observations give a
#'                    warning. Default 60.
#' @param ... further arguments to \code{\link[goFlux]{find.rise}}; set
#'            \code{rise}, \code{drop} and \code{conc.range} when
#'            \code{gastype} is not CO2 in ppm.
#'
#' @returns A data.frame in the format of \code{\link[goFlux]{click.peak2}}
#'          (all windows bound together, with \code{flag}, \code{Etime},
#'          \code{start.time_corr}, \code{end.time_corr} and
#'          \code{obs.length_corr}), with an attribute \code{"log"}: a
#'          data.frame (\code{UniqueID}, \code{method}, \code{auto.start},
#'          \code{auto.end}, \code{n.flag}, \code{dconc}) describing each
#'          choice.
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{find.rise}}, \code{\link[goFlux]{click.peak2}},
#'          \code{\link[goFlux]{windows.from.table}}
#'
#' @examples
#' set.seed(1)
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' tt <- t0 + 0:900
#' ow <- list(data.frame(UniqueID = "tree1", POSIX.time = tt,
#'   CO2dry_ppm = c(rep(410, 200), 410 + 0.3 * (0:400), rep(530, 300)) + rnorm(901, 0, 0.3),
#'   CH4dry_ppb = 2000 + 0.01 * (0:900) + rnorm(901, 0, 0.5),
#'   start.time = t0 + 180, obs.length = 600))
#' manID <- auto.id.rise(ow)
#' attr(manID, "log")
#' @export
auto.id.rise <- function(ow.list, gastype = "CO2dry_ppm", gases = NULL,
                         end.time = NULL, trim.secs = 30, max.secs = 600,
                         min.secs = 60, slack.secs = 60, fallback.secs = 300,
                         plot.dir = NULL, tz = "UTC", warn.length = 60, ...) {
  if(is.data.frame(ow.list)) ow.list <- split(ow.list, ow.list$UniqueID)
  if(!is.list(ow.list) || length(ow.list) == 0) stop("'ow.list' must be a non-empty list of data.frames")
  en <- NULL
  if(!is.null(end.time)){
    if(is.data.frame(end.time)){
      if(!all(c("UniqueID", "end.time") %in% names(end.time))){
        stop("a data.frame 'end.time' must contain 'UniqueID' and 'end.time'")}
      e <- end.time$end.time
      if(!inherits(e, "POSIXct")) e <- as.POSIXct(e, tz = tz)
      en <- e; names(en) <- end.time$UniqueID
    } else {
      en <- end.time
      if(!inherits(en, "POSIXct")) en <- as.POSIXct(en, tz = tz)
    }
  }
  if(!is.null(plot.dir)) dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
  if(!grepl("^CO2.*_ppm$", gastype) && !all(c("rise", "drop", "conc.range") %in% names(list(...)))){
    warning("the defaults of 'rise', 'drop' and 'conc.range' are for CO2 in ppm; ",
            "set them for '", gastype, "'", call. = FALSE)}

  out <- vector("list", length(ow.list)); log <- vector("list", length(ow.list))
  for(k in seq_along(ow.list)){
    d <- as.data.frame(ow.list[[k]]); uid <- unique(d$UniqueID)[1]
    if(!"start.time" %in% names(d)) stop(uid, ": no 'start.time' column")
    st <- unique(d$start.time)[1]
    if(is.na(st)) stop(uid, ": 'start.time' is NA")
    if(!gastype %in% names(d)) stop("column '", gastype, "' not found in 'ow.list'")
    en.uid <- if(!is.null(en) && uid %in% names(en) && !is.na(en[[uid]])) en[[uid]] else NA
    r <- find.rise(d$POSIX.time, d[[gastype]], ...)
    method <- "auto_rise"
    if(is.null(r)){
      s0 <- st + trim.secs
      e0 <- if(!is.na(en.uid)) en.uid else st + fallback.secs
      method <- "notebook_fallback"
    } else {
      s0 <- r$start + trim.secs
      e0 <- r$end
      if(!is.na(en.uid)) e0 <- min(e0, en.uid + slack.secs)
      e0 <- min(e0, s0 + max.secs)
    }
    if(as.numeric(e0 - s0, units = "secs") < min.secs){
      e0 <- s0 + min.secs + 30; method <- paste(method, "short")}
    d2 <- flag.window(d, s0, e0)
    n.flag <- sum(d2$flag)
    if(n.flag < warn.length){
      warning(uid, ": only ", n.flag, " observations in window (< warn.length = ",
              warn.length, ")", call. = FALSE)}
    out[[k]] <- d2
    log[[k]] <- data.frame(UniqueID = uid, method = method,
                           auto.start = format(s0, "%Y-%m-%d %H:%M:%S"),
                           auto.end = format(e0, "%Y-%m-%d %H:%M:%S"),
                           n.flag = n.flag,
                           dconc = if(is.null(r)) NA_real_ else round(r$dconc, 1),
                           stringsAsFactors = FALSE)
    if(!is.null(plot.dir)){
      fl <- d2$flag == 1
      g <- if(is.null(gases)) unique(c(gastype, intersect(c("CH4dry_ppb", "N2Odry_ppb"), names(d))))
           else intersect(gases, names(d))
      grDevices::png(file.path(plot.dir, paste0(sprintf("%03d", k), "_", uid, ".png")),
                     900, 350 * length(g))
      op <- par(mfrow = c(length(g), 1), mar = c(3, 4.5, 2.5, 1))
      for(i in seq_along(g)){
        graphics::plot(d$POSIX.time, d[[g[i]]], col = ifelse(fl, "red", "grey40"), pch = 16,
                       cex = 0.5, ylab = g[i], xlab = "",
                       main = if(i == 1) paste(uid, "-", method) else "")
        graphics::abline(v = st, col = "blue")
        if(!is.na(en.uid)) graphics::abline(v = en.uid, col = "blue", lty = 2)
      }
      par(op); grDevices::dev.off()
    }
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  attr(res, "log") <- do.call(rbind, log)
  res
}

# Flag one observation window the way click.peak2() does: flag (1 inside
# [start, end]), Etime (s since start), start.time_corr, end.time_corr and
# obs.length_corr.
flag.window <- function(d, start, end) {
  if(end <= start) stop("the window end must be after its start")
  d$flag <- as.numeric(d$POSIX.time >= start & d$POSIX.time <= end)
  d$Etime <- as.numeric(d$POSIX.time - start, units = "secs")
  d$start.time_corr <- start
  d$end.time_corr <- end
  d$obs.length_corr <- as.numeric(end - start, units = "secs")
  d
}
